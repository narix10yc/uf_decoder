#include "UF/UFDecoder3D.h"
#include "Utils/PrintSpan.h"
#include <iostream>
#include <list>

using namespace uf;

static bool isInValidEdgeRange(int idx, int L, int T) {
  return 0 <= idx && idx < (3 * L * L * T - L * L);
}

static bool isInValidVertexRange(int idx, int L, int T) {
  return 0 <= idx && idx < (L * L * T);
}

EdgeCoordinate3D::EdgeCoordinate3D(int idx, int L, int T) {
  assert(isInValidEdgeRange(idx, L, T));
  t = (idx / (3 * L * L));
  idx %= (3 * L * L);
  if (idx < L * L) {
    group = Horizontal;
  } else if (idx < 2 * L * L) {
    group = Vertical;
    idx -= L * L;
  } else {
    group = Temporal;
    idx -= 2 * L * L;
  }
  r = idx / L;
  c = idx % L;
}

std::ostream& uf::operator<<(std::ostream& os, const EdgeCoordinate3D& coord) {
  os << "(";
  switch (coord.group) {
  case EdgeCoordinate3D::Horizontal:
    os << "H";
    break;
  case EdgeCoordinate3D::Vertical:
    os << "V";
    break;
  case EdgeCoordinate3D::Temporal:
    os << "T";
    break;
  default:
    assert(false && "Unknown edge group");
    os << "Unknown";
  }
  os << ", " << coord.t << ", " << coord.r << ", " << coord.c << ")";
  return os;
}

std::ostream& uf::operator<<(std::ostream& os, const BitArray& arr) {
  int b = 0;
  while (b < arr.nbits()) {
    os << (arr.test(b) ? '1' : '0');
    ++b;
    if (b % 8 == 0)
      os << ' ';
    if (b % 64 == 0)
      os << '\n';
  }
  os << '\n';
  return os;
}

namespace {
// (t,r,c)
struct VertexCoord {
  int t, r, c;
  int idx;
};

std::ostream& operator<<(std::ostream& os, const VertexCoord& coord) {
  return os << "(" << coord.t << "," << coord.r << "," << coord.c << ")";
}

// A sparse-dense set data structure to store cluster boundary vertices
// This structure is move-only
class ClusterBoundaryVertices {
  // These two arrays will never to resized
  std::unique_ptr<VertexCoord[]> values_;
  std::unique_ptr<int[]> indices_;
  unsigned valueSize_;
  unsigned size_;

  static constexpr int INVALID = -1;

public:
  ClusterBoundaryVertices(unsigned size)
      : values_(std::make_unique<VertexCoord[]>(size)),
        indices_(std::make_unique<int[]>(size)), valueSize_(0), size_(size) {
    std::fill_n(indices_.get(), size, -1);
  }

  bool contains(const VertexCoord& coord) const {
    assert(0 <= coord.idx && coord.idx < size_);
    return indices_[coord.idx] != INVALID;
  }

  // returns true of inserted; false if already exists
  bool insert(const VertexCoord& coord) {
    if (contains(coord))
      return false; // already exists
    values_[valueSize_] = coord;
    indices_[coord.idx] = valueSize_;
    ++valueSize_;
    return true;
  }

  // returns true if removed; false if not found
  bool erase(const VertexCoord& coord) {
    if (!contains(coord))
      return false; // not found
    const int idx = indices_[coord.idx];
    assert(idx != INVALID);
    assert(idx < valueSize_);
    // Move the last element to the position of the removed element
    --valueSize_;
    values_[idx] = values_[valueSize_];
    indices_[values_[idx].idx] = idx;
    indices_[coord.idx] = INVALID;
    return true;
  }

  // returns the number of values present
  unsigned size() const { return valueSize_; }

  // returns true if popped; false if empty. If popped, set is updated.
  bool pop(VertexCoord& coord) {
    if (valueSize_ == 0)
      return false;
    --valueSize_;
    coord = values_[valueSize_];
    indices_[coord.idx] = INVALID;
    return true;
  }
};

class ClusterDSU {
  std::unique_ptr<int[]> parent;
  std::unique_ptr<int[]> size;

public:
  ClusterDSU(int n)
      : parent(std::make_unique<int[]>(n)), size(std::make_unique<int[]>(n)) {
    for (int i = 0; i < n; i++)
      parent[i] = i;
    std::fill_n(size.get(), n, 1);
  }

  int find(int x) { return parent[x] == x ? x : parent[x] = find(parent[x]); }

  // returns the root of the merged tree. If they are already in the same set,
  // returns -1;
  int unite(int a, int b) {
    a = find(a);
    b = find(b);
    if (a == b)
      return -1;

    if (size[a] < size[b])
      std::swap(a, b); // ensure a is bigger
    parent[b] = a;
    size[a] += size[b];
    return a; // return the root of the merged tree
  }

  int set_size(int x) { return size[find(x)]; }
};

struct Int3 {
  int t, r, c;
};

// route is given by
// - horizontal edges (t0, r0, c0) to (t0, r0, c1), followed by
// - vertical edges (t0, r0, c1) to (t0, r1, c1), followed by
// - temporal edges (t0, r1, c1) to (t1, r1, c1).
void findRouteAndSetErasure(BitArray& erasure,
                            const VertexCoord& A,
                            const VertexCoord& B,
                            int L) {
  Int3 route{A.t, A.r, A.c};
  int x1; // the target coordinate component
// horizontal edges
label_horizontal:
  assert(route.c == A.c);
  if (A.c == B.c) {
    // nothing to do on horizontal edges
    goto label_vertical;
  } else if (A.c < B.c) {
    route.c = A.c;
    x1 = B.c;
  } else {
    route.c = B.c;
    x1 = A.c;
  }
  assert(route.c < x1);
  // we will go from route.c to x1
  if (x1 - route.c < L / 2) {
    // go right
    assert(route.c < x1);
    for (; route.c < x1; ++route.c) {
      erasure.set(edgeIdx_H(route.t, route.r, route.c, L));
    }
  } else {
    // go left
    std::swap(route.c, x1);
    assert(route.c > x1);
    for (; route.c > x1; --route.c) {
      erasure.set(edgeIdx_H(route.t, route.r, route.c, L));
    }
  }
label_vertical:
  // we might have gone from B to A
  route.c = B.c;
  assert(route.r == A.r);
  if (A.r == B.r) {
    // nothing to do on vertical edges
    goto label_temporal;
  } else if (A.r < B.r) {
    route.r = A.r;
    x1 = B.r;
  } else {
    route.r = B.r;
    x1 = A.r;
  }
  // we will go from route.r to x1
  if (x1 - route.r < L / 2) {
    // go down
    assert(route.r < x1);
    for (; route.r < x1; ++route.r) {
      erasure.set(edgeIdx_V(route.t, route.r, route.c, L));
    }
  } else {
    // go up
    std::swap(route.r, x1);
    assert(route.r > x1);
    for (; route.r > x1; --route.r) {
      erasure.set(edgeIdx_V(route.t, route.r, route.c, L));
    }
  }
label_temporal:
  // we might have gone from B to A
  route.r = B.r;
  assert(route.t == A.t);
  if (A.t == B.t) {
    // nothing to do on temporal edges
    return;
  } else if (A.t < B.t) {
    route.t = A.t;
    x1 = B.t;
  } else {
    route.t = B.t;
    x1 = A.t;
  }
  // there is no wrap around in the temporal axis.
  // always go from route.t to x1
  for (; route.t < x1; ++route.t) {
    erasure.set(edgeIdx_T(route.t, route.r, route.c, L));
  }
}

} // anonymous namespace

void uf::clustering(BitArray& erasure,
                    const BitArray& syndromes,
                    int L,
                    int T) {
  erasure.unset_all();

  // LLT equals to the number of vertices
  const int LLT = L * L * T;
  assert(syndromes.nbits() == LLT);

  ClusterDSU dsu(LLT);

  // initialize erasure
  erasure = syndromes;

  enum ClusterStatus { Even, Odd, Invalid };

  struct Cluster {
    VertexCoord root;
    ClusterBoundaryVertices boundaryVertices;
    ClusterStatus status;

    inline void toggleParity() {
      assert(status != Invalid);
      status = (status == Even) ? Odd : Even;
    }

    // reserveSize should be set to the total number of vertices
    Cluster(const VertexCoord& root, unsigned reserveSize)
        : root(root), boundaryVertices(reserveSize), status(Odd) {
      boundaryVertices.insert(root);
    }
  };

  // create a circular reference:
  // - cluster to boundary vertices via the clusters variable
  // - vertex to cluster via clusterLookup

  // clusters with odd parity
  std::list<Cluster> clusters;
  // indexed by vertex indices. only valid for root vertices
  std::unique_ptr<Cluster*[]> lookup_(std::make_unique<Cluster*[]>(LLT));
  const auto lookupCluster = [&](int vertex) -> Cluster* {
    assert(isInValidVertexRange(vertex, L, T));
    return lookup_[dsu.find(vertex)];
  };

  std::fill_n(lookup_.get(), LLT, nullptr);

  const auto mergeClusters = [&](Cluster* A, Cluster* B) -> void {
    std::cerr << "Merging clusters with roots " << A->root << " and " << B->root
              << "\n";
    auto mergedRoot = dsu.unite(A->root.idx, B->root.idx);
    std::cerr << "Merged root is " << mergedRoot << "\n";
    assert(mergedRoot != -1);
    // invalidate the smaller cluster
    // ensure B is the larger one
    if (mergedRoot != B->root.idx)
      std::swap(A, B);
    assert(mergedRoot == A->root.idx);
    // - move the boundary vertices to the bigger cluster
    VertexCoord movedVertex;
    while (B->boundaryVertices.pop(movedVertex)) {
      A->boundaryVertices.insert(movedVertex);
    }
    // - invalidate the smaller cluster
    A->status = Invalid;
    B->toggleParity();
    // - find a route between the two roots and add the route to the erasure
    findRouteAndSetErasure(erasure, A->root, B->root, L);
  };

  for (int idx : syndromes.ones()) {
    // idx is a vertex index
    // pre-compute all (t,r,c) coordinates for future uses
    const int t = idx / (L * L);
    const int r = (idx / L) % L;
    const int c = idx % L;

    clusters.emplace_back(VertexCoord{t, r, c, idx}, LLT);
    lookup_[idx] = &clusters.back();
    // ensure the circular reference
    assert(lookup_[idx]->root.idx == idx);
  }

  while (!clusters.empty()) {
    // find the smallest-size cluster to grow
    auto minIt = clusters.begin();
    auto minSize = minIt->boundaryVertices.size();
    for (auto it = minIt, end = clusters.end(); it != end; ++it) {
      // skip invalid cluster. skip even parity clusters.
      if (it->status != Odd)
        continue;
      if (it->boundaryVertices.size() < minSize) {
        minIt = it;
        minSize = it->boundaryVertices.size();
      }
    }

    // no more cluster to grow
    if (minIt->status != Odd)
      break;

    // grow the cluster
    Cluster* cluster = &(*minIt);
    Cluster* growCluster = nullptr;
    assert(cluster != nullptr);
    VertexCoord vertex;
    VertexCoord growVertex;
    while (cluster->boundaryVertices.pop(vertex)) {
      // up vertex (t,r-1,c)
      growVertex.t = vertex.t;
      growVertex.r = (vertex.r == 0) ? (L - 1) : (vertex.r - 1);
      growVertex.c = vertex.c;
      growVertex.idx = trc(growVertex.t, growVertex.r, growVertex.c, L);
      growCluster = lookupCluster(growVertex.idx);
      if (growCluster == nullptr) {
        // growVertex is not in any cluster. add it to the boundary
        cluster->boundaryVertices.insert(growVertex);
      } else if (growCluster == cluster) {
        // Grown into the same cluster's boundary vertices. This could happen
        // when the cluster has grown a whole loop in this direction (up-down or
        // left-right). Just do nothing
      } else {
        assert(growCluster->status != Invalid);
        // Grown into another cluster's boundary vertices. Merge them.
        mergeClusters(cluster, growCluster);
        // re-evaluate which cluster to grow
        break;
      }
      // down vertex (t,r+1,c)
      growVertex.t = vertex.t;
      growVertex.r = (vertex.r == L - 1) ? 0 : (vertex.r + 1);
      growVertex.c = vertex.c;
      growVertex.idx = trc(growVertex.t, growVertex.r, growVertex.c, L);
      growCluster = lookupCluster(growVertex.idx);
      if (growCluster == nullptr) {
        cluster->boundaryVertices.insert(growVertex);
      } else if (growCluster == cluster) {
      } else {
        assert(growCluster->status != Invalid);
        mergeClusters(cluster, growCluster);
        break;
      }
      // left vertex (t,r,c-1)
      growVertex.t = vertex.t;
      growVertex.r = vertex.r;
      growVertex.c = (vertex.c == 0) ? (L - 1) : (vertex.c - 1);
      growVertex.idx = trc(growVertex.t, growVertex.r, growVertex.c, L);
      growCluster = lookupCluster(growVertex.idx);
      if (growCluster == nullptr) {
        cluster->boundaryVertices.insert(growVertex);
      } else if (growCluster == cluster) {
      } else {
        assert(growCluster->status != Invalid);
        mergeClusters(cluster, growCluster);
        break;
      }
      // right vertex (t,r,c+1)
      growVertex.t = vertex.t;
      growVertex.r = vertex.r;
      growVertex.c = (vertex.c == L - 1) ? 0 : (vertex.c + 1);
      growVertex.idx = trc(growVertex.t, growVertex.r, growVertex.c, L);
      growCluster = lookupCluster(growVertex.idx);
      if (growCluster == nullptr) {
        cluster->boundaryVertices.insert(growVertex);
      } else if (growCluster == cluster) {
      } else {
        assert(growCluster->status != Invalid);
        mergeClusters(cluster, growCluster);
        break;
      }
      // forward vertex (t+1,r,c)
      if (vertex.t < T - 1) {
        growVertex.t = vertex.t + 1;
        growVertex.r = vertex.r;
        growVertex.c = vertex.c;
        growVertex.idx = trc(growVertex.t, growVertex.r, growVertex.c, L);
        growCluster = lookupCluster(growVertex.idx);
        if (growCluster == nullptr) {
          cluster->boundaryVertices.insert(growVertex);
        } else if (growCluster == cluster) {
          assert(false && "Should never happen in the temporal axis");
        } else {
          assert(growCluster->status != Invalid);
          mergeClusters(cluster, growCluster);
          break;
        }
      }
      // backward vertex (t-1,r,c)
      if (vertex.t > 0) {
        growVertex.t = vertex.t - 1;
        growVertex.r = vertex.r;
        growVertex.c = vertex.c;
        growVertex.idx = trc(growVertex.t, growVertex.r, growVertex.c, L);
        growCluster = lookupCluster(growVertex.idx);
        if (growCluster == nullptr) {
          cluster->boundaryVertices.insert(growVertex);
        } else if (growCluster == cluster) {
          assert(false && "Should never happen in the temporal axis");
        } else {
          assert(growCluster->status != Invalid);
          mergeClusters(cluster, growCluster);
          break;
        }
      }
    } // while pop a vertex from the boundary
  } // while clusters not empty
}

void uf::peelingDecode(BitArray& correction,
                       const BitArray& erasure,
                       BitArray syndromes,
                       int L,
                       int T) {
  const int LL = L * L;
  const int LL2 = 2 * LL;
  const int LL3 = 3 * LL;
  const int LLT = L * L * T;

  assert(correction.nbits() == 3 * LLT - LL);
  assert(erasure.nbits() == 3 * LLT - LL);
  assert(syndromes.nbits() == LLT);

  // (LL * t) + (L * r) + (c). No range check
  const auto trc = [L](int t, int r, int c) -> int {
    return t * L * L + r * L + c;
  };

  // dsu and degrees are to be indexed by syndromes (vertices)
  // total number of vertices is L * L * T
  DisjointSetUnion dsu(LLT);
  auto degrees = std::make_unique<uint8_t[]>(LLT);
  std::fill_n(degrees.get(), LLT, static_cast<uint8_t>(0));

  // forest is to be indexed by edges
  BitArray forest(erasure.nbits());
  correction.unset_all();

  /* ---- Step 1: build spanning forest ---- */

  // edge (H,t,r,c) connects vertex (t,r,c) and vertex (t,r,c+1)
  for (int t = 0; t < T; ++t) {
    for (int r = 0; r < L; ++r) {
      for (int c = 0; c < L - 1; ++c) {
        const int edge = edgeIdx_H(t, r, c, L); // edge (H,t,r,c)
        const int vertexA = trc(t, r, c);       // vertex (t,r,c)
        const int vertexB = vertexA + 1;        // vertex (t,r,c+1)

        assert(isInValidEdgeRange(edge, L, T));
        assert(isInValidVertexRange(vertexA, L, T));
        assert(isInValidVertexRange(vertexB, L, T));

        if (erasure.test(edge) && dsu.unite(vertexA, vertexB)) {
          forest.flip(edge);
          degrees[vertexA] += 1;
          degrees[vertexB] += 1;
        }
      } // for c = 0 .. L-2
      // c = L - 1 case
      const int edge = edgeIdx_H(t, r, L - 1, L); // edge (H,t,r,L-1)
      const int vertexA = trc(t, r, L - 1);       // vertex (t,r,L-1)
      const int vertexB = trc(t, r, 0);           // vertex (t,r,0)

      assert(isInValidEdgeRange(edge, L, T));
      assert(isInValidVertexRange(vertexA, L, T));
      assert(isInValidVertexRange(vertexB, L, T));

      if (erasure.test(edge) && dsu.unite(vertexA, vertexB)) {
        forest.flip(edge);
        degrees[vertexA] += 1;
        degrees[vertexB] += 1;
      }
    } // for r = 0 .. L-1
  } // for t = 0 .. T-1

  // edge (V,t,r,c) connects vertex (t,r,c) and vertex (t,r+1,c)
  for (int t = 0; t < T; ++t) {
    for (int c = 0; c < L; ++c) {
      for (int r = 0; r < L - 1; ++r) {
        const int edge = edgeIdx_V(t, r, c, L); // edge (V,t,r,c)
        const int vertexA = trc(t, r, c);       // vertex (t,r,c)
        const int vertexB = vertexA + L;        // vertex (t,r+1,c)

        assert(isInValidEdgeRange(edge, L, T));
        assert(isInValidVertexRange(vertexA, L, T));
        assert(isInValidVertexRange(vertexB, L, T));

        if (erasure.test(edge) && dsu.unite(vertexA, vertexB)) {
          forest.flip(edge);
          degrees[vertexA] += 1;
          degrees[vertexB] += 1;
        }
      } // for r = 0 .. L-2
      // r = L - 1 case
      const int edge = edgeIdx_V(t, L - 1, c, L); // edge (V,t,L-1,c)
      const int vertexA = trc(t, L - 1, c);       // vertex (t,L-1,c)
      const int vertexB = trc(t, 0, c);           // vertex (t,0,c)

      assert(isInValidEdgeRange(edge, L, T));
      assert(isInValidVertexRange(vertexA, L, T));
      assert(isInValidVertexRange(vertexB, L, T));

      if (erasure.test(edge) && dsu.unite(vertexA, vertexB)) {
        forest.flip(edge);
        degrees[vertexA] += 1;
        degrees[vertexB] += 1;
      }
    } // for c = 0 .. L-1
  } // for t = 0 .. T-1

  // edge (T,t,r,c) connects vertex (t,r,c) and vertex (t+1,r,c)
  for (int t = 0; t < T - 1; ++t) {
    for (int r = 0; r < L; ++r) {
      for (int c = 0; c < L; ++c) {
        const int edge = edgeIdx_T(t, r, c, L); // edge (T,t,r,c)
        const int vertexA = trc(t, r, c);       // vertex (t,r,c)
        const int vertexB = vertexA + LL;       // vertex (t+1,r,c)

        assert(isInValidEdgeRange(edge, L, T));
        assert(isInValidVertexRange(vertexA, L, T));
        assert(isInValidVertexRange(vertexB, L, T));

        if (erasure.test(edge) && dsu.unite(vertexA, vertexB)) {
          forest.flip(edge);
          degrees[vertexA] += 1;
          degrees[vertexB] += 1;
        }
      } // for c = 0 .. L-1
    } // for r = 0 .. L-1
  } // for t = 0 .. T-2

  /* ---- Step 2: identify leaves ---- */

  struct Leaf {
    int edge;
    int pendantNode;
    // We will need the trc information to the non-pendant node
    int nonPendantNode_t;
    int nonPendantNode_r;
    int nonPendantNode_c;
  };

  // Given a pendant node (t,r,c), update the leaf information
  // by finding the edge and the non-pendant node.
  const auto updateLeaf = [=](int t, int r, int c, Leaf& leaf) {
    assert(leaf.pendantNode == trc(t, r, c));
    int edge;
    // up edge (V,t,r-1,c)
    edge = edgeIdx_V(t, (r == 0) ? (L - 1) : (r - 1), c, L);
    if (forest.test(edge)) {
      leaf.edge = edge;
      // up star (t,r-1,c)
      leaf.nonPendantNode_t = t;
      leaf.nonPendantNode_r = (r == 0) ? (L - 1) : (r - 1);
      leaf.nonPendantNode_c = c;
      return;
    }
    // down edge (V,t,r,c)
    edge = edgeIdx_V(t, r, c, L);
    if (forest.test(edge)) {
      leaf.edge = edge;
      // down star (t,r+1,c)
      leaf.nonPendantNode_t = t;
      leaf.nonPendantNode_r = (r == L - 1) ? 0 : (r + 1);
      leaf.nonPendantNode_c = c;
      return;
    }
    // left edge (H,t,r,c-1)
    edge = edgeIdx_H(t, r, (c == 0) ? (L - 1) : (c - 1), L);
    if (forest.test(edge)) {
      leaf.edge = edge;
      // left star (t,r,c-1)
      leaf.nonPendantNode_t = t;
      leaf.nonPendantNode_r = r;
      leaf.nonPendantNode_c = (c == 0) ? (L - 1) : (c - 1);
      return;
    }
    // right edge (H,t,r,c)
    edge = edgeIdx_H(t, r, c, L);
    if (forest.test(edge)) {
      leaf.edge = edge;
      // right star (t,r,c+1)
      leaf.nonPendantNode_t = t;
      leaf.nonPendantNode_r = r;
      leaf.nonPendantNode_c = (c == L - 1) ? 0 : (c + 1);
      return;
    }
    // forward edge (T,t,r,c)
    if (t < T - 1) {
      edge = edgeIdx_T(t, r, c, L);
      if (forest.test(edge)) {
        leaf.edge = edge;
        // forward star (t+1,r,c)
        leaf.nonPendantNode_t = t + 1;
        leaf.nonPendantNode_r = r;
        leaf.nonPendantNode_c = c;
        return;
      }
    }
    // backward edge (T,t-1,r,c)
    assert(t > 0);
    edge = edgeIdx_T(t - 1, r, c, L);
    assert(forest.test(edge));
    leaf.edge = edge;
    // backward star (t-1,r,c)
    leaf.nonPendantNode_t = t - 1;
    leaf.nonPendantNode_r = r;
    leaf.nonPendantNode_c = c;
  }; // updateLeaf

  // leaves are to be indexed by stars (vertices)
  auto leaves = std::make_unique<Leaf[]>(LLT);
  int leafIdx = 0;
  for (int t = 0; t < T; ++t) {
    for (int r = 0; r < L; ++r) {
      for (int c = 0; c < L; ++c) {
        const int star = trc(t, r, c);
        if (degrees[star] != 1)
          continue;
        // star is now a pendant node
        // search for the edge
        auto& leaf = leaves[leafIdx];
        leaf.pendantNode = star;
        updateLeaf(t, r, c, leaf);

        // Avoid adding repeated leaves
        const int otherStar = trc(leaf.nonPendantNode_t,
                                  leaf.nonPendantNode_r,
                                  leaf.nonPendantNode_c);
        assert(degrees[otherStar] > 0);
        if (degrees[otherStar] > 1) {
          // accept this leaf
          ++leafIdx;
          continue;
        }
        // We have a singleton tree
        // mark the other pendant node as having no degree to skip
        // subsequent search
        degrees[otherStar] = 0;
        if (syndromes.test(star)) {
          assert(syndromes.test(otherStar));
          assert(!correction.test(leaf.edge));
          correction.set(leaf.edge);

          // no need to update syndromes at these two stars
          // as we won't use them anymore
          // syndromes.unset(star);
          // syndromes.unset(otherStar);
        } else {
          assert(!syndromes.test(otherStar));
        }
      } // for c = 0 .. L-1
    } // for r = 0 .. L-1
  } // for t = 0 .. T-1

  // Step 3: Peeling
  std::cerr << "There are " << leafIdx << " pendant nodes in the forest.\n";

  while (leafIdx > 0) {
    const auto [edge, pStar, npStar_t, npStar_r, npStar_c] = leaves[--leafIdx];
    const int npStar = trc(npStar_t, npStar_r, npStar_c);

    assert(forest.test(edge));
    assert(degrees[pStar] == 1);
    assert(degrees[npStar] >= 1);

    // Fast path: if edge is a singleton leaf (the last leaf of this tree)
    if (degrees[npStar] == 1) {
      // std::cerr << "Fast path: edge is a singleton tree.\n";
      if (syndromes.test(pStar)) {
        assert(syndromes.test(npStar));
        assert(!correction.test(edge));
        correction.set(edge);
        // no need to update syndromes at these two stars
        // as we won't use them anymore
        // syndromes.unset(pendantStar);
        // syndromes.unset(nonPendantStar);
      }
      forest.unset(edge);
      degrees[pStar] = 0;
      degrees[npStar] = 0;
      continue;
    } // end of fast path

    // Substep 1: Remove edge from the forest and update leaves.
    // DSU is not designed for removal. We won't update dsu.
    forest.unset(edge);
    // no need to update pendantStar's degree as we won't use it anymore
    // degrees[pendantStar] = 0;

    auto& leaf = leaves[leafIdx];
    if (degrees[npStar] == 1) {
      // nonPendantStar is now a pendant node
      // search for the leaf edge. update the leaf in-place
      leaf.pendantNode = npStar;
      updateLeaf(npStar_t, npStar_r, npStar_c, leaf);

      // avoid repeating the same leaf
      // only accept the new leaf if the new pendant node has >1 degree
      const auto new_npStar = trc(
          leaf.nonPendantNode_t, leaf.nonPendantNode_r, leaf.nonPendantNode_c);
      if (degrees[new_npStar] > 1)
        ++leafIdx;
    }

    // Substep 2: Update correction and syndromes.
    if (syndromes.test(pStar)) {
      correction.flip(leaf.edge);
      syndromes.flip(npStar);
    }
  } // while leafIdx > 0
}
