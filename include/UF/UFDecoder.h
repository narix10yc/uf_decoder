#ifndef UF_UFDECODER_H
#define UF_UFDECODER_H

#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <memory>
#include <queue>
#include <vector>

#include "ToricCodeVisualizer.h"

// #define UF_ENABLE_PEELING_VISUALIZATION

// clang-format off
#ifdef UF_ENABLE_PEELING_VISUALIZATION
#define DEBUG_PEELING_DECODER(X) do { X } while (0)
#else
#define DEBUG_PEELING_DECODER(X) {}
#endif
// clang-format on
namespace uf {

// Disjoint Set Union
template <int L> struct DSU {
  static constexpr int LL = L * L;
  std::unique_ptr<int[]> parents;
  std::unique_ptr<uint8_t[]> ranks;

  DSU() {
    parents = std::make_unique<int[]>(LL);
    ranks = std::make_unique<uint8_t[]>(LL);
    for (int i = 0; i < LL; ++i) {
      parents[i] = i;
      ranks[i] = 0;
    }
  }

  int find(int x) {
    while (parents[x] != x) {
      parents[x] = parents[parents[x]]; // path compression
      x = parents[x];
    }
    return x;
  }

  // Unite sets of a and b; return true if merged
  bool unite(int a, int b) {
    a = find(a);
    b = find(b);
    if (a == b) {
      // already in the same tree
      return false;
    }
    // ensure a has larger rank
    if (ranks[a] < ranks[b])
      std::swap(a, b);
    parents[b] = a;
    if (ranks[a] == ranks[b])
      ++ranks[a];
    return true;
  }
};

struct StarNeighbour {
  int upEdge, downEdge, leftEdge, rightEdge;
  int upStar, downStar, leftStar, rightStar;
};

// returns r * L + c
template <int L> inline constexpr int rc_(int r, int c) {
  assert(0 <= r && r < L);
  assert(0 <= c && c < L);
  return r * L + c; // row-major index
}

// The indices for adjacent edges and stars for a star
template <int L> std::unique_ptr<StarNeighbour[]> build_star_adjacencies() {
  constexpr int LL = L * L;
  auto adjacencies = std::make_unique<StarNeighbour[]>(LL);
  constexpr auto rc = [](int r, int c) { return rc_<L>(r, c); };

  // Star (vertex) (r, c) is connected to:
  // - edges:
  // up edge    (V, r - 1, c    ),
  // down edge  (V, r,     c    ),
  // left edge  (H, r,     c - 1),
  // right edge (H, r,     c    ),
  // - stars:
  // up star    (r - 1, c    ),
  // down star  (r + 1, c    ),
  // left star  (r,     c - 1),
  // right star (r,     c + 1).

  // clang-format off

  // First Row
  // (r,c) = (0,0)
  adjacencies[0].upEdge    = rc(L-1, 0) + LL;
  adjacencies[0].downEdge  = rc(0, 0) + LL;
  adjacencies[0].leftEdge  = rc(0, L-1);
  adjacencies[0].rightEdge = rc(0, 0);
  adjacencies[0].upStar    = rc(0, 0);
  adjacencies[0].downStar  = rc(1, 0);
  adjacencies[0].leftStar  = rc(0, L-1);
  adjacencies[0].rightStar = rc(0, 1);

  // (r,c) = (0, 1 .. L-2)
  for (int c = 1; c < L-1; ++c) {
    adjacencies[rc(0, c)].upEdge    = rc(L-1, c) + LL;
    adjacencies[rc(0, c)].downEdge  = rc(0, c) + LL;
    adjacencies[rc(0, c)].leftEdge  = rc(0, c-1);
    adjacencies[rc(0, c)].rightEdge = rc(0, c);
    adjacencies[rc(0, c)].upStar    = rc(L-1, c);
    adjacencies[rc(0, c)].downStar  = rc(1, c);
    adjacencies[rc(0, c)].leftStar  = rc(0, c-1);
    adjacencies[rc(0, c)].rightStar = rc(0, c+1);
  }

  // (r,c) = (0, L-1)
  adjacencies[rc(0, L-1)].upEdge    = rc(L-1, L-1) + LL;
  adjacencies[rc(0, L-1)].downEdge  = rc(0, L-1) + LL;
  adjacencies[rc(0, L-1)].leftEdge  = rc(0, L-2);
  adjacencies[rc(0, L-1)].rightEdge = rc(0, L-1);
  adjacencies[rc(0, L-1)].upStar    = rc(L-1, L-1);
  adjacencies[rc(0, L-1)].downStar  = rc(1, L-1);
  adjacencies[rc(0, L-1)].leftStar  = rc(0, L-2);
  adjacencies[rc(0, L-1)].rightStar = rc(0, 0);

  // Middle Rows
  for (int r = 1; r < L-1; ++r) {
    // (r,c) = (1 .. L-2, 0)
    adjacencies[rc(r, 0)].upEdge    = rc(r-1, 0) + LL;
    adjacencies[rc(r, 0)].downEdge  = rc(r, 0) + LL;
    adjacencies[rc(r, 0)].leftEdge  = rc(r, L-1);
    adjacencies[rc(r, 0)].rightEdge = rc(r, 0);
    adjacencies[rc(r, 0)].upStar    = rc(r-1, 0);
    adjacencies[rc(r, 0)].downStar  = rc(r+1, 0);
    adjacencies[rc(r, 0)].leftStar  = rc(r, L-1);
    adjacencies[rc(r, 0)].rightStar = rc(r, 1);

    // (r,c) = (1 .. L-2, 1 .. L-2)
    for (int c = 1; c < L-1; ++c) {
      adjacencies[rc(r, c)].upEdge    = rc(r-1, c) + LL;
      adjacencies[rc(r, c)].downEdge  = rc(r, c) + LL;
      adjacencies[rc(r, c)].leftEdge  = rc(r, c-1);
      adjacencies[rc(r, c)].rightEdge = rc(r, c);
      adjacencies[rc(r, c)].upStar    = rc(r-1, c);
      adjacencies[rc(r, c)].downStar  = rc(r+1, c);
      adjacencies[rc(r, c)].leftStar  = rc(r, c-1);
      adjacencies[rc(r, c)].rightStar = rc(r, c+1); 
    }

    // (r,c) = (1 .. L-2, L-1)
    adjacencies[rc(r, L-1)].upEdge    = rc(r-1, L-1) + LL;
    adjacencies[rc(r, L-1)].downEdge  = rc(r, L-1) + LL;
    adjacencies[rc(r, L-1)].leftEdge  = rc(r, L-2);
    adjacencies[rc(r, L-1)].rightEdge = rc(r, L-1);
    adjacencies[rc(r, L-1)].upStar    = rc(r-1, L-1);
    adjacencies[rc(r, L-1)].downStar  = rc(r+1, L-1);
    adjacencies[rc(r, L-1)].leftStar  = rc(r, L-2);
    adjacencies[rc(r, L-1)].rightStar = rc(r, 0);
  }

  // Last Row
  // (r,c) = (L-1,0)
  adjacencies[rc(L-1, 0)].upEdge    = rc(L-2, 0) + LL;
  adjacencies[rc(L-1, 0)].downEdge  = rc(L-1, 0) + LL;
  adjacencies[rc(L-1, 0)].leftEdge  = rc(L-1, L-1);
  adjacencies[rc(L-1, 0)].rightEdge = rc(L-1, 0);
  adjacencies[rc(L-1, 0)].upStar    = rc(L-1, 0);
  adjacencies[rc(L-1, 0)].downStar  = rc(0, 0);
  adjacencies[rc(L-1, 0)].leftStar  = rc(L-1, L-1);
  adjacencies[rc(L-1, 0)].rightStar = rc(L-1, 1);

  // (r,c) = (L-1, 1 .. L-2)
  for (int c = 1; c < L-1; ++c) {
    adjacencies[rc(L-1, c)].upEdge    = rc(L-2, c) + LL;
    adjacencies[rc(L-1, c)].downEdge  = rc(L-1, c) + LL;
    adjacencies[rc(L-1, c)].leftEdge  = rc(L-1, c-1);
    adjacencies[rc(L-1, c)].rightEdge = rc(L-1, c);
    adjacencies[rc(L-1, c)].upStar    = rc(L-2, c);
    adjacencies[rc(L-1, c)].downStar  = rc(0, c);
    adjacencies[rc(L-1, c)].leftStar  = rc(L-1, c-1);
    adjacencies[rc(L-1, c)].rightStar = rc(L-1, c+1);
  }

  // (r,c) = (L-1, L-1)
  adjacencies[rc(L-1, L-1)].upEdge    = rc(L-2, L-1) + LL;
  adjacencies[rc(L-1, L-1)].downEdge  = rc(L-1, L-1) + LL;
  adjacencies[rc(L-1, L-1)].leftEdge  = rc(L-1, L-2);
  adjacencies[rc(L-1, L-1)].rightEdge = rc(L-1, L-1);
  adjacencies[rc(L-1, L-1)].upStar    = rc(L-2, L-1);
  adjacencies[rc(L-1, L-1)].downStar  = rc(0, L-1);
  adjacencies[rc(L-1, L-1)].leftStar  = rc(L-1, L-2);
  adjacencies[rc(L-1, L-1)].rightStar = rc(L-1, 0);

  // clang-format on

  return adjacencies;
} // build_adjacencies

// template <int L>
// constexpr std::array<ToricCodeError, 2 * L * L> build_star_qubit_adjacency()
// {
//   // Star (vertex) (r, c) is connected to:
//   // up qubit    (V, r - 1, c    ),
//   // down qubit  (V, r,     c    ),
//   // left qubit  (H, r,     c - 1),
//   // right qubit (H, r,     c    ),

//   using Code = ToricCode<L>;

//   std::array<ToricCodeError, 2 * L * L> adjacency;
//   for (int r = 0; r < L; ++r) {
//     for (int c = 0; c < L; ++c) {
//       const int starIdx = rc_<L>(r, c);
//       // up qubit
//       adjacency[starIdx].set(Code::qubitIdx_VG(r - 1, c));
//       // down qubit
//       adjacency[starIdx].set(Code::qubitIdx_VG(r, c));
//       // left qubit
//       adjacency[starIdx].set(Code::qubitIdx_HG(r, c - 1));
//       // right qubit
//       adjacency[starIdx].set(Code::qubitIdx_HG(r, c));
//     }
//   }
//   return adjacency;
// }

// template<int L>
// constexpr std::array<ToricCodeSyndrome, L*L> build_star_star_adjacency() {
//   // Star (vertex) (r, c) is connected to:
//   // up star    (r - 1, c    ),
//   // down star  (r + 1, c    ),
//   // left star  (r,     c - 1),
//   // right star (r,     c + 1).
//   constexpr auto rc = [](int r, int c) {

//   }
//   std::array<ToricCodeSyndrome, L * L> adjacency;
//   for (int r = 0; r < L; ++r) {
//     for (int c = 0; c < L; ++c) {
//       const int starIdx = rc_<L>(r, c);
// }

// Return the erasure set
template <int L>
ToricCode<L>::Error cluster(const typename ToricCode<L>::Syndrome xSyndrome,
                            const StarNeighbour* adjacencies) {
  assert(xSyndrome.parity() == 0 && "Syndrome must have even parity");
  return {};
}

template <int L>
ToricCode<L>::Error
peeling_decode_Z(const typename ToricCode<L>::Error& erasure,
                 typename ToricCode<L>::Syndrome xSyndrome,
                 const StarNeighbour* adjacencies) {
  using Error = typename ToricCode<L>::Error;
  using Syndrome = typename ToricCode<L>::Syndrome;
  constexpr auto rc = [](int r, int c) { return rc_<L>(r, c); };
  static constexpr int LL = L * L;

  ToricCodeVisualizer<L> visualizer;

  // dsu stores the spanning forest of the syndrome graph (a graph with LL
  // vertices). dsu is to be indexed by stars
  DSU<L> dsu;
  auto degrees = std::make_unique<uint8_t[]>(LL);

  // zero-initialized. forest is to be indexed by qubits
  Error forest{};

  // Step 1: Build the spanning forest

  // edge (H,r,c) connects vertex (r,c) and vertex (r,c+1)
  for (int r = 0; r < L; ++r) {
    for (int c = 0; c < L - 1; ++c) {
      const int edge = rc(r, c);       // edge (H,r,c)
      const int vertexA = rc(r, c);    // vertex (r,c)
      const int vertexB = vertexA + 1; // vertex (r,c+1)
      if (erasure.test(edge) && dsu.unite(vertexA, vertexB)) {
        forest.flip(edge);
        degrees[vertexA] += 1;
        degrees[vertexB] += 1;
      }
    } // for c = 0 .. L-2
    // c = L - 1 case
    const int edge = rc(r, L - 1);    // edge (H,r,L-1)
    const int vertexA = rc(r, L - 1); // vertex (r,L-1)
    const int vertexB = rc(r, 0);     // vertex (r,0)
    if (erasure.test(edge) && dsu.unite(vertexA, vertexB)) {
      forest.flip(edge);
      degrees[vertexA] += 1;
      degrees[vertexB] += 1;
    }
  }

  // edge (V,r,c) connects vertex (r,c) and vertex (r+1,c)
  for (int c = 0; c < L; ++c) {
    for (int r = 0; r < L - 1; ++r) {
      const int edge = LL + rc(r, c);   // edge (V,r,c)
      const int vertexA = rc(r, c);     // vertex (r,c)
      const int vertexB = rc(r + 1, c); // vertex (r+1,c)
      if (erasure.test(edge) && dsu.unite(vertexA, vertexB)) {
        forest.flip(edge);
        degrees[vertexA] += 1;
        degrees[vertexB] += 1;
      }
    } // for r = 0 .. L-2
    // r = L - 1 case
    const int edge = LL + rc(L - 1, c); // edge (V,L-1,c)
    const int vertexA = rc(L - 1, c);   // vertex (L-1,c)
    const int vertexB = rc(0, c);       // vertex (0,c)
    if (erasure.test(edge) && dsu.unite(vertexA, vertexB)) {
      forest.flip(edge);
      degrees[vertexA] += 1;
      degrees[vertexB] += 1;
    }
  }

  // Step 2: Initialize leaf nodes

  // zero-initialized
  Error correction{};
  struct Leaf {
    int edge;
    int pendantNode;
    int nonPendantNode;
  };

  auto leaves = std::make_unique<Leaf[]>(LL);
  int leafIdx = 0;
  for (int star = 0; star < LL; ++star) {
    if (degrees[star] != 1)
      continue;
    // star is now a pendant node
    // search for the edge
    int edge;
    int otherStar;
    const auto& adj = adjacencies[star];
    if (forest.test(adj.upEdge)) {
      edge = adj.upEdge;
      otherStar = adj.upStar;
    } else if (forest.test(adj.downEdge)) {
      edge = adj.downEdge;
      otherStar = adj.downStar;
    } else if (forest.test(adj.leftEdge)) {
      edge = adj.leftEdge;
      otherStar = adj.leftStar;
    } else {
      assert(forest.test(adj.rightEdge));
      edge = adj.rightEdge;
      otherStar = adj.rightStar;
    }

    // Avoid adding repeated leaves
    assert(degrees[otherStar] > 0);
    if (degrees[otherStar] > 1) {
      leaves[leafIdx++] = {edge, star, otherStar};
      continue;
    }
    // We have a singleton tree
    // mark the other pendant node as having no degree to skip
    // subsequent search
    degrees[otherStar] = 0;
    if (xSyndrome.test(star)) {
      assert(xSyndrome.test(otherStar));
      assert(!correction.test(edge));
      correction.set(edge);

      // no need to update syndromes at these two stars
      // as we won't use them anymore
      // xSyndrome.unset(star);
      // xSyndrome.unset(otherStar);
    } else {
      assert(!xSyndrome.test(otherStar));
    }
  }

  std::ofstream ofs;
  int fileIdx = 0;
  DEBUG_PEELING_DECODER(
      std::cerr << "There are " << leafIdx << " pendant nodes in the forest.\n";
      visualizer.clearDocument();
      visualizer.drawLattice();
      visualizer.drawError(forest, "#4A0");
      visualizer.drawSyndromeX(xSyndrome);
      ofs.open("peeling." + std::to_string(fileIdx++) + ".svg");
      visualizer.render(ofs);
      ofs.close(););

  // Step 3: Peeling

  while (leafIdx > 0) {
    const auto [edge, pendantStar, nonPendantStar] = leaves[--leafIdx];
    DEBUG_PEELING_DECODER(
        std::cerr << "Step " << fileIdx << ": [edge, pStar, npStar] = [" << edge
                  << ", " << pendantStar << ", " << nonPendantStar << "]\n";);
    assert(forest.test(edge));
    assert(degrees[pendantStar] == 1);
    assert(degrees[nonPendantStar] >= 1);

    // Fast path: if edge is a singleton leaf (the last leaf of this tree)
    if (degrees[nonPendantStar] == 1) {
      DEBUG_PEELING_DECODER(std::cerr
                                << "Fast path: edge is a singleton tree.\n";);
      if (xSyndrome.test(pendantStar)) {
        assert(xSyndrome.test(nonPendantStar));
        assert(!correction.test(edge));
        correction.set(edge);
        // no need to update syndromes at these two stars
        // as we won't use them anymore
        // xSyndrome.unset(pendantStar);
        // xSyndrome.unset(nonPendantStar);
      }
      forest.unset(edge);
      degrees[pendantStar] = 0;
      degrees[nonPendantStar] = 0;
      goto label_continue;
    }

    // Substep 1: Remove edge from the forest and update leaves.
    // DSU is not designed for removal. We won't update dsu.
    forest.unset(edge);
    // no need to update pendantStar's degree as we won't use it anymore
    // degrees[pendantStar] = 0;

    if (degrees[nonPendantStar] == 1) {
      // nonPendantStar is now a pendant node
      // search for the leaf edge. update the leaf in-place
      const auto& adj = adjacencies[nonPendantStar];
      if (forest.test(adj.upEdge))
        leaves[leafIdx] = {adj.upEdge, nonPendantStar, adj.upStar};
      else if (forest.test(adj.downEdge))
        leaves[leafIdx] = {adj.downEdge, nonPendantStar, adj.downStar};
      else if (forest.test(adj.leftEdge))
        leaves[leafIdx] = {adj.leftEdge, nonPendantStar, adj.leftStar};
      else {
        assert(forest.test(adj.rightEdge));
        leaves[leafIdx] = {adj.rightEdge, nonPendantStar, adj.rightStar};
      }

      // avoid repeating the same leaf
      if (degrees[leaves[leafIdx].nonPendantNode] > 1)
        ++leafIdx;
    }

    // Substep 2: Update correction and syndromes.
    if (xSyndrome.test(pendantStar)) {
      correction.flip(edge);
      xSyndrome.flip(nonPendantStar);
    }

  label_continue:
    DEBUG_PEELING_DECODER(
        visualizer.clearDocument(); visualizer.drawLattice();
        visualizer.drawError(forest, "#4A0");
        visualizer.drawSyndromeX(xSyndrome);
        ofs.open("peeling." + std::to_string(fileIdx++) + ".svg");
        visualizer.render(ofs);
        ofs.close(););
    continue;

  } // while leafIdx > 0

  return correction;

} // peeling_decode_Z

// Horizontal edges: [0 .. L*L-1], edge (r,c) connects (r,c)->(r,c+1)
// Vertical edges:   [L*L .. 2*L*L-1], edge (r,c) connects (r,c)->(r+1,c)
//
// Given a Z operator on data edges (Error), return two bits
// (b_h, b_v) where:
//   b_h = 1 iff it anticommutes with the horizontal X logical (a loop of X
//         on the horizontal edges of some fixed row r0)
//   b_v = 1 iff it anticommutes with the vertical X logical (a loop of X
//         on the vertical edges of some fixed column c0)
//
// Any choice of row r0 and column c0 works (homology invariance).
// @return bit 0 gives the horizontal parity, bit 1 gives the vertical parity.
template <int L>
int compute_logical_op(typename ToricCode<L>::Error const& zOp,
                       int r0 = 0,
                       int c0 = 0) {
  int b_h = 0; // parity vs X_h (horizontal loop on row r0)
  for (int c = 0; c < L; ++c) {
    int hIdx = r0 * L + c; // horizontal edge (r0, c)
    b_h ^= zOp.test(hIdx) ? 1 : 0;
  }

  int b_v = 0; // parity vs X_v (vertical loop on col c0)
  for (int r = 0; r < L; ++r) {
    int vIdx = L * L + r * L + c0; // vertical edge (r, c0)
    b_v ^= zOp.test(vIdx) ? 1 : 0;
  }

  return b_h | (b_v << 1); // return as a single int
} // compute_logical_op

} // namespace uf

#endif // UF_UFDECODER_H