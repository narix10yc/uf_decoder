#pragma once
#include <array>
#include <cassert>
#include <cstdint>
#include <queue>
#include <vector>

// This decoder depends on your cleaned-up Toric code header.
#include "ToricCode.h"

// ==============================================================
// Union-Find 3D (space–time) decoder for Toric code (periodic in space)
//
// Decodes Z errors from a stream of X-syndrome measurements over time.
// Time is OPEN by default (first/last measurements act as boundaries).
// IMPORTANT: When all inter-round XORs are zero (identical snapshots),
//            there are no 3D detection events. In that case we **fall back**
//            to a 2D union-find decode on a single snapshot (by default the
//            last snapshot). This handles the "perfect measurement" case.
// ==============================================================

namespace uf {

template <int L> class UF3D_ToricDecoder {
  // ---- Aliases to the Toric code types ----
  using Code = ToricCode<L>;

public:
  using ErrBits = typename Code::ErrBitArray; // 2 * L * L
  using SynBits = typename Code::SynBitArray; // L * L

private:
  // ---- Space–time lattice sizes ----
  static constexpr int V = L * L; // checks per spatial layer

  // ---- Indexing helpers for nodes (t, r, c) ----
  static inline int modL(int x) {
    int m = x % L;
    return m < 0 ? m + L : m;
  }
  static inline int node_index(int t, int r, int c, int Tlayers) {
    (void)Tlayers;
    return t * V + (r * L + c);
  }
  static inline void node_coords(int idx, int Tlayers, int& t, int& r, int& c) {
    (void)Tlayers;
    t = idx / V;
    int rc = idx % V;
    r = rc / L;
    c = rc % L;
  }

  // Map two adjacent star checks to a data-qubit *global* index (0..2LL-1).
  static int data_qubit_from_adjacent_stars(int r1, int c1, int r2, int c2) {
    r1 = modL(r1);
    c1 = modL(c1);
    r2 = modL(r2);
    c2 = modL(c2);
    // Horizontal neighbors
    if (r1 == r2) {
      int dc = (c2 - c1 + L) % L;
      if (dc == 1)
        return /*H(r1,c1)*/ (r1 * L + c1);
      if (dc == L - 1)
        return /*H(r2,c2)*/ (r2 * L + c2);
    }
    // Vertical neighbors
    if (c1 == c2) {
      int dr = (r2 - r1 + L) % L;
      if (dr == 1)
        return /*V(r1,c1)*/ (V + r1 * L + c1);
      if (dr == L - 1)
        return /*V(r2,c2)*/ (V + r2 * L + c2);
    }
    assert(false && "stars not adjacent on torus");
    return -1;
  }

  // ---- State for a single decode call ----
  int T_{};       // number of measurement rounds
  int Tlayers_{}; // = T_ - 1 detection-time layers
  int N_{};       // total nodes = Tlayers_ * V

  // detection events (1 = odd vertex)
  std::vector<uint8_t> det_; // size N_

  // BFS ownership/parent pointers
  std::vector<int>
      owner_; // size N_, -1 if unclaimed; else DSU representative id
  std::vector<int> parent_; // size N_, predecessor node index
  enum Step : uint8_t { NONE = 0, LEFT, RIGHT, UP, DOWN, TIME_BACK, TIME_FWD };
  std::vector<Step> parent_step_; // size N_

  // Disjoint set union over cluster representatives + two time boundaries
  std::vector<int> dsu_parent_; // size N_ + 2
  std::vector<uint8_t> dsu_rank_;
  std::vector<uint8_t> dsu_odd_;    // parity: 1 if cluster currently odd
  std::vector<uint8_t> dsu_active_; // cache: grows only while odd

  // Pseudo DSU ids for time boundaries
  int B_LOWER_{}; // below t=0
  int B_UPPER_{}; // above t=Tlayers_-1

  struct Meet {
    int a;
    int b;
    Step via;
    bool boundary;
    int which; /*0=lower,1=upper*/
  };
  std::vector<Meet> meets_;

  struct QItem {
    int node;
  };
  std::queue<QItem> q_;

  // ---- DSU helpers ----
  int find(int x) {
    if (dsu_parent_[x] == x)
      return x;
    return dsu_parent_[x] = find(dsu_parent_[x]);
  }
  int unite(int a, int b) {
    a = find(a);
    b = find(b);
    if (a == b)
      return a;
    if (dsu_rank_[a] < dsu_rank_[b])
      std::swap(a, b);
    dsu_parent_[b] = a;
    if (dsu_rank_[a] == dsu_rank_[b])
      ++dsu_rank_[a];
    dsu_odd_[a] ^= dsu_odd_[b];
    dsu_active_[a] = dsu_odd_[a];
    return a;
  }

  // ==========================================================
  // (1) Build detection events from a time series of syndromes
  // ==========================================================
  void build_detection_events(const std::vector<SynBits>& syndromes) {
    T_ = static_cast<int>(syndromes.size());
    Tlayers_ = T_ - 1;
    N_ = Tlayers_ * V;
    det_.assign(N_, 0);
    for (int t = 0; t < Tlayers_; ++t) {
      for (int v = 0; v < V; ++v) {
        const bool bit = (syndromes[t + 1].test(v) ^ syndromes[t].test(v));
        det_[node_index(t, v / L, v % L, Tlayers_)] = static_cast<uint8_t>(bit);
      }
    }
  }

  bool any_detection() const {
    for (auto b : det_)
      if (b)
        return true;
    return false;
  }

  // ==============================================
  // (2) Initialize clusters and union-find forest
  // ==============================================
  void initialize_clusters() {
    owner_.assign(N_, -1);
    parent_.assign(N_, -1);
    parent_step_.assign(N_, NONE);

    dsu_parent_.resize(N_ + 2);
    dsu_rank_.assign(N_ + 2, 0);
    dsu_odd_.assign(N_ + 2, 0);
    dsu_active_.assign(N_ + 2, 0);

    B_LOWER_ = N_ + 0; // time below first layer
    B_UPPER_ = N_ + 1; // time above last layer
    dsu_parent_[B_LOWER_] = B_LOWER_;
    dsu_parent_[B_UPPER_] = B_UPPER_;

    // Initialize all nodes and enqueue detection-event seeds
    for (int idx = 0; idx < N_; ++idx) {
      dsu_parent_[idx] = idx;
      if (det_[idx]) {
        owner_[idx] = idx; // seed owns itself
        dsu_odd_[idx] = 1; // odd cluster
        dsu_active_[idx] = 1;
        q_.push({idx});
      }
    }
  }

  // ==============================================
  // (3) Growth: multi-source BFS + DSU merges
  // ==============================================
  void record_boundary_meet(int u, int which_boundary) {
    // Record that the frontier from node u touched a time boundary; this
    // toggles the parity of its cluster.
    meets_.push_back({u,
                      -1,
                      which_boundary == 0 ? TIME_BACK : TIME_FWD,
                      true,
                      which_boundary});
    int ru = find(owner_[u]);
    if (!dsu_active_[ru])
      return;          // already even
    dsu_odd_[ru] ^= 1; // absorb one defect into the boundary
    dsu_active_[ru] = dsu_odd_[ru];
    // Optional: merge with boundary representative for bookkeeping
    unite(ru, which_boundary == 0 ? B_LOWER_ : B_UPPER_);
  }

  void try_edge(int u, int v, Step step) {
    if (v < 0 || v >= N_)
      return; // defensive
    int ru = find(owner_[u]);
    if (!dsu_active_[ru])
      return; // cluster already satisfied

    if (owner_[v] == -1) {
      // Claim fresh node for this cluster
      owner_[v] = ru;
      parent_[v] = u;
      parent_step_[v] = step;
      q_.push({v});
    } else {
      // Frontier meets another (possibly same) cluster
      int rv = find(owner_[v]);
      if (ru != rv) {
        meets_.push_back({u, v, step, false, -1});
        unite(ru, rv);
      }
    }
  }

  void grow_until_even() {
    while (!q_.empty()) {
      int u = q_.front().node;
      q_.pop();
      // Seeds are guaranteed to have owners; skip if cluster already even
      int ru = find(owner_[u]);
      if (!dsu_active_[ru])
        continue;

      int t, r, c;
      node_coords(u, Tlayers_, t, r, c);

      // --- Temporal neighbors ---
      if (t == 0)
        record_boundary_meet(u, /*lower*/ 0);
      else {
        int v = node_index(t - 1, r, c, Tlayers_);
        if (owner_[v] == -1 || find(owner_[v]) != ru)
          try_edge(u, v, TIME_BACK);
      }

      if (t == Tlayers_ - 1)
        record_boundary_meet(u, /*upper*/ 1);
      else {
        int v = node_index(t + 1, r, c, Tlayers_);
        if (owner_[v] == -1 || find(owner_[v]) != ru)
          try_edge(u, v, TIME_FWD);
      }

      // --- Spatial neighbors (periodic wrap) ---
      {
        int v = node_index(t, r, (c + L - 1) % L, Tlayers_);
        if (owner_[v] == -1 || find(owner_[v]) != ru)
          try_edge(u, v, LEFT);
      }
      {
        int v = node_index(t, r, (c + 1) % L, Tlayers_);
        if (owner_[v] == -1 || find(owner_[v]) != ru)
          try_edge(u, v, RIGHT);
      }
      {
        int v = node_index(t, (r + L - 1) % L, c, Tlayers_);
        if (owner_[v] == -1 || find(owner_[v]) != ru)
          try_edge(u, v, UP);
      }
      {
        int v = node_index(t, (r + 1) % L, c, Tlayers_);
        if (owner_[v] == -1 || find(owner_[v]) != ru)
          try_edge(u, v, DOWN);
      }
    }
  }

  // ==============================================
  // (4) Peeling: emit a spatial Z-correction
  // ==============================================
  static inline void
  toggle_spatial_edge(ErrBits& corr, int r1, int c1, int r2, int c2) {
    int q = data_qubit_from_adjacent_stars(r1, c1, r2, c2);
    corr.flip(q);
  }

  // Follow parent pointers from u back towards a seed (detection event) and
  // toggle spatial edges along the way. Returns the seed node index on success,
  // or -1 if no seed was encountered (shouldn’t happen for properly rooted
  // trees).
  int peel_path_to_source(int u, ErrBits& corr) {
    while (u != -1 && !det_[u]) {
      int p = parent_[u];
      Step s = parent_step_[u];
      if (p == -1)
        break; // reached a root that is not a detection (defensive)
      int t0, r0, c0;
      node_coords(p, Tlayers_, t0, r0, c0);
      int t1, r1, c1;
      node_coords(u, Tlayers_, t1, r1, c1);
      if (s == LEFT || s == RIGHT || s == UP || s == DOWN) {
        toggle_spatial_edge(corr, r0, c0, r1, c1);
      }
      u = p;
    }
    return (u != -1 && det_[u]) ? u : -1;
  }

  void peel_and_emit_correction(ErrBits& corr) {
    for (const auto& m : meets_) {
      if (m.boundary) {
        (void)peel_path_to_source(m.a,
                                  corr); // boundary contributes no spatial edge
      } else {
        (void)peel_path_to_source(m.a, corr);
        (void)peel_path_to_source(m.b, corr);
        if (m.via == LEFT || m.via == RIGHT || m.via == UP || m.via == DOWN) {
          int ta, ra, ca;
          node_coords(m.a, Tlayers_, ta, ra, ca);
          int tb, rb, cb;
          node_coords(m.b, Tlayers_, tb, rb, cb);
          toggle_spatial_edge(corr, ra, ca, rb, cb);
        }
      }
    }
  }

  // =============================================================
  // 2D fallback: union-find on a single snapshot (perfect measurement)
  // =============================================================
  ErrBits decode2D_from_snapshot(const SynBits& syn) {
    // Data structures local to 2D
    std::vector<uint8_t> det2(V, 0);
    std::vector<int> owner2(V, -1), parent2(V, -1);
    enum Step2 : uint8_t { N2 = 0, L2, R2, U2, D2 };
    std::vector<uint8_t> pstep2(V, N2);
    std::vector<int> dsu_p(V), dsu_r(V, 0), dsu_parity(V, 0), dsu_active(V, 0);

    struct Meet2 {
      int a;
      int b;
      Step2 via;
    };
    std::vector<Meet2> meets2;
    std::queue<int> q2;

    auto fnd = [&](auto self, int x) -> int {
      return dsu_p[x] == x ? x : dsu_p[x] = self(self, dsu_p[x]);
    };
    auto uft = [&](int a, int b) {
      a = fnd(fnd, a);
      b = fnd(fnd, b);
      if (a == b)
        return a;
      if (dsu_r[a] < dsu_r[b])
        std::swap(a, b);
      dsu_p[b] = a;
      if (dsu_r[a] == dsu_r[b])
        ++dsu_r[a];
      dsu_parity[a] ^= dsu_parity[b];
      dsu_active[a] = dsu_parity[a];
      return a;
    };

    // Seeds
    for (int v = 0; v < V; ++v) {
      dsu_p[v] = v;
      if (syn.test(v)) {
        det2[v] = 1;
        owner2[v] = v;
        dsu_parity[v] = 1;
        dsu_active[v] = 1;
        q2.push(v);
      }
    }

    auto rc_of = [&](int idx) { return std::pair<int, int>{idx / L, idx % L}; };
    auto idx_of = [&](int r, int c) { return r * L + c; };

    auto try_edge2 = [&](int u, int v, Step2 via) {
      int ru = fnd(fnd, u);
      if (!dsu_active[ru])
        return; // satisfied cluster
      if (owner2[v] == -1) {
        owner2[v] = ru;
        parent2[v] = u;
        pstep2[v] = via;
        q2.push(v);
      } else {
        int rv = fnd(fnd, owner2[v]);
        if (ru != rv) {
          meets2.push_back({u, v, via});
          uft(ru, rv);
        }
      }
    };

    // BFS growth on 4-neighborhood with toric wrap
    while (!q2.empty()) {
      int u = q2.front();
      q2.pop();
      int ru = fnd(fnd, u);
      if (!dsu_active[ru])
        continue;
      auto [r, c] = rc_of(u);
      try_edge2(u, idx_of(r, (c + L - 1) % L), L2);
      try_edge2(u, idx_of(r, (c + 1) % L), R2);
      try_edge2(u, idx_of((r + L - 1) % L, c), U2);
      try_edge2(u, idx_of((r + 1) % L, c), D2);
    }

    // Peeling
    ErrBits corr; // zeroed
    auto toggle_spatial = [&](int u, int v, Step2 s) {
      auto [r0, c0] = rc_of(u);
      auto [r1, c1] = rc_of(v);
      int q = data_qubit_from_adjacent_stars(r0, c0, r1, c1);
      corr.flip(q);
    };

    auto peel_to_source = [&](int u) {
      while (u != -1 && !det2[u]) {
        int p = parent2[u];
        Step2 s = (Step2)pstep2[u];
        if (p == -1)
          break;
        toggle_spatial(p, u, s);
        u = p;
      }
      return (u != -1 && det2[u]) ? u : -1;
    };

    for (const auto& m : meets2) {
      (void)peel_to_source(m.a);
      (void)peel_to_source(m.b);
      toggle_spatial(m.a, m.b, m.via);
    }
    return corr;
  }

public:
  UF3D_ToricDecoder() = default;

  // Main entry: decode Z errors from a stream of X-syndrome rounds.
  // Input: syndromes[t] is the X-syndrome (stars) at round t, t=0..T-1.
  // Output: ErrBits mask of Z corrections over data qubits (2*L*L bits).
  // If all inter-round XORs are zero, falls back to 2D UF on the last snapshot.
  ErrBits decode(const std::vector<SynBits>& syndromes) {
    ErrBits corr; // zero-initialized
    const int T = static_cast<int>(syndromes.size());
    if (T <= 0)
      return corr;
    if (T == 1) {
      // Single round → pure 2D decode
      return decode2D_from_snapshot(syndromes.back());
    }

    // Phases
    build_detection_events(syndromes);
    initialize_clusters();
    grow_until_even();
    peel_and_emit_correction(corr);
    return corr;
  }
};

template <int L>
ToricCode<L>::ErrBitArray
decode(const std::vector<typename ToricCode<L>::SynBitArray>& syndromes) {
  UF3D_ToricDecoder<L> decoder;
  return decoder.decode(syndromes);
}
} // namespace uf
