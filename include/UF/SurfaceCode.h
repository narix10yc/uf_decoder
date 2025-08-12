#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <random>

#include "SVG/Document.h"

namespace uf {

//==============================================================
// StaticBitArray<Nbits>: fixed-size bitset backed by 64-bit words
//==============================================================
template <int Nbits> class StaticBitArray {
  static_assert(Nbits > 0, "Nbits must be positive");
  static constexpr int kWords = (Nbits + 63) / 64;
  std::array<uint64_t, kWords> w_{}; // zero-initialized

  // Mark unused bits
  static constexpr uint64_t tail_mask() {
    if constexpr ((Nbits & 63) == 0)
      return ~0ULL;
    else
      return (~0ULL) >> (64 - (Nbits & 63));
  }

public:
  static constexpr int size_bits() { return Nbits; }
  static constexpr int size_words() { return kWords; }

  void clear() { w_.fill(0ULL); }

  void set_all() {
    w_.fill(~0ULL);
    w_.back() &= tail_mask();
  }

  void set(int i) { w_[i >> 6] |= (1ULL << (i & 63)); }
  void reset_bit(int i) { w_[i >> 6] &= ~(1ULL << (i & 63)); }
  void flip(int i) { w_[i >> 6] ^= (1ULL << (i & 63)); }
  bool test(int i) const { return (w_[i >> 6] >> (i & 63)) & 1ULL; }

  uint64_t word(int i) const { return w_[i]; }
  uint64_t* words() { return w_.data(); }
  const uint64_t* words() const { return w_.data(); }

  void xor_with(const StaticBitArray& other) {
    for (int i = 0; i < kWords; ++i)
      w_[i] ^= other.w_[i];
  }

  static inline int parity64(uint64_t x) {
#if defined(__GNUG__) || defined(__clang__)
    return __builtin_parityll(x);
#else
    x ^= x >> 32;
    x ^= x >> 16;
    x ^= x >> 8;
    x ^= x >> 4;
    x &= 0xF;
    return (0x6996 >> x) & 1;
#endif
  }

  int masked_parity(const StaticBitArray& mask) const {
    int p = 0;
    for (int i = 0; i < kWords; ++i)
      p ^= parity64(w_[i] & mask.w_[i]);
    return p & 1;
  }

  void display(std::ostream& os) const {
    int b = 0;
    while (b < Nbits) {
      os << (test(b) ? '1' : '0');
      ++b;
      if (b % 8 == 0)
        os << ' ';
      if (b % 64 == 0)
        os << '\n';
    }
    os << '\n';
  }
};

//==============================================================
// SurfaceCode<L>: compile-time L, runtime periodic vs planar
//==============================================================
template <int L> class SurfaceCode {
  static_assert(L >= 2, "L must be >= 2");

public:
  // Maximum number of data qubits occurs for the toric code: 2*L*L
  static constexpr int NH_TORIC = L * L; // horizontal edges if periodic
  static constexpr int NV_TORIC = L * L; // vertical edges if periodic
  // Toric codes has 2*L*L qubits
  static constexpr int NQ_MAX = 2 * L * L;
  using Bits = StaticBitArray<NQ_MAX>;

  // private:
public:
  bool periodic_; // true: toric, false: planar

  // Qubits are divided into two groups:
  // horizontal edges (H) and vertical edges (V).

  // horizontal lattice: H_rows_ = L; H_cols_ = periodic? L : (L-1)
  int H_rows_, H_cols_;
  // vertical lattice:   V_rows_ = periodic? L : (L-1); V_cols_ = L
  int V_rows_, V_cols_;
  // number of *usable* H/V edges for this geometry
  int Nh_used_, Nv_used_;
  // Nh_used_ + Nv_used_ = total number of usable qubits
  int Nq_used_{};

  // Error bitfields (zero-initialized)
  Bits xErr_{};
  Bits zErr_{};

  // Stabilizer incidence masks (always sized to L*L "slots").
  static constexpr int N_STARS = L * L;                  // vertices
  static constexpr int N_PLAQ_TORIC = L * L;             // faces in toric
  static constexpr int N_PLAQ_PLANR = (L - 1) * (L - 1); // faces in planar

  std::array<Bits, N_STARS> starMasks_{}; // X-checks: detect Z errors
  std::array<Bits, N_PLAQ_TORIC>
      plaquetteMasks_{}; // Z-checks: detect X errors; last (L*L - (L-1)^2)
                         // entries empty for planar

  // ---------------- Indexing helpers ----------------
  // Horizontal edges: between (r,c) and (r,c+1)
  static constexpr int hIndex_max(int r, int c) {
    return r * L + c;
  } // uses L columns (toric layout)
  // Vertical edges: between (r,c) and (r+1,c)
  static constexpr int vIndex_max(int r, int c) {
    return NH_TORIC + r * L + c;
  } // uses L rows (toric layout)

  static constexpr int modL(int x) {
    int m = x % L;
    return m < 0 ? m + L : m;
  }

  inline void add_h_edge(Bits& m, int r, int c) const {
    if (periodic_) {
      m.set(hIndex_max(modL(r), modL(c)));
    } else {
      if (0 <= r && r < L && 0 <= c && c < (L - 1))
        m.set(hIndex_max(r, c));
    }
  }
  inline void add_v_edge(Bits& m, int r, int c) const {
    if (periodic_) {
      m.set(vIndex_max(modL(r), modL(c)));
    } else {
      if (0 <= r && r < (L - 1) && 0 <= c && c < L)
        m.set(vIndex_max(r, c));
    }
  }

  void build_masks() {
    // --- Plaquettes ---
    int nP = periodic_ ? N_PLAQ_TORIC : N_PLAQ_PLANR;
    int rMax = periodic_ ? L : (L - 1);
    int cMax = periodic_ ? L : (L - 1);
    int idx = 0;
    for (int r = 0; r < rMax; ++r) {
      for (int c = 0; c < cMax; ++c) {
        Bits m; // zero-initialized
        // face (r,c): top H(r,c), right V(r,c+1), bottom H(r+1,c), left V(r,c)
        add_h_edge(m, r, c);
        add_v_edge(m, r, c + 1);
        add_h_edge(m, r + 1, c);
        add_v_edge(m, r, c);
        plaquetteMasks_[idx++] = m;
      }
    }
    // Zero-out any unused trailing mask slots (planar)
    for (; idx < N_PLAQ_TORIC; ++idx)
      plaquetteMasks_[idx].clear();

    // --- Stars ---
    // (always LxL vertices; boundary stars have fewer incident edges)
    for (int r = 0; r < L; ++r) {
      for (int c = 0; c < L; ++c) {
        Bits m;
        // up, down, left, right
        add_v_edge(m, r - 1, c);
        add_v_edge(m, r, c);
        add_h_edge(m, r, c - 1);
        add_h_edge(m, r, c);
        starMasks_[r * L + c] = m;
      }
    }
  }

public:
  explicit SurfaceCode(bool isPeriodic) : periodic_(isPeriodic) {
    H_rows_ = L;
    H_cols_ = periodic_ ? L : (L - 1);
    V_rows_ = periodic_ ? L : (L - 1);
    V_cols_ = L;
    Nh_used_ = H_rows_ * H_cols_;
    Nv_used_ = V_rows_ * V_cols_;
    Nq_used_ = Nh_used_ + Nv_used_;
    build_masks();
  }

  // ------------ Introspection ------------
  static constexpr int gridSize() { return L; }
  bool isPeriodic() const { return periodic_; }
  int numQubits() const { return Nq_used_; }
  int numPlaquettes() const { return periodic_ ? N_PLAQ_TORIC : N_PLAQ_PLANR; }
  int numStars() const { return N_STARS; }

  // ------------ Error manipulation ------------
  void clearAllErr() {
    xErr_.clear();
    zErr_.clear();
  }
  void setXErr(int q) {
    assert(0 <= q && q < NQ_MAX);
    xErr_.set(q);
  }
  void setZErr(int q) {
    assert(0 <= q && q < NQ_MAX);
    zErr_.set(q);
  }
  void flipXErr(int q) {
    assert(0 <= q && q < NQ_MAX);
    xErr_.flip(q);
  }
  void flipZErr(int q) {
    assert(0 <= q && q < NQ_MAX);
    zErr_.flip(q);
  }
  bool hasXErr(int q) const {
    assert(0 <= q && q < NQ_MAX);
    return xErr_.test(q);
  }
  bool hasZErr(int q) const {
    assert(0 <= q && q < NQ_MAX);
    return zErr_.test(q);
  }

  template <class URNG> void injectPauliNoise(double pX, double pZ, URNG& rng) {
    assert(0.0 <= pX && pX <= 1.0);
    assert(0.0 <= pZ && pZ <= 1.0);
    std::bernoulli_distribution bx(pX), bz(pZ);
    for (int q = 0; q < Nq_used_; ++q) {
      if (bx(rng))
        xErr_.flip(q);
      if (bz(rng))
        zErr_.flip(q);
    }
  }
  template <class URNG> void injectXNoise(double p, URNG& rng) {
    std::bernoulli_distribution b(p);
    for (int q = 0; q < Nq_used_; ++q)
      if (b(rng))
        xErr_.flip(q);
  }
  template <class URNG> void injectZNoise(double p, URNG& rng) {
    std::bernoulli_distribution b(p);
    for (int q = 0; q < Nq_used_; ++q)
      if (b(rng))
        zErr_.flip(q);
  }

  // ------------ Syndromes ------------
  StaticBitArray<N_PLAQ_TORIC> computeZSyndrome() const {
    // extra tail entries will remain 0 for planar
    StaticBitArray<N_PLAQ_TORIC> syn;
    const int nP = numPlaquettes();
    for (int i = 0; i < nP; ++i)
      if (xErr_.masked_parity(plaquetteMasks_[i]))
        syn.set(i);
    return syn;
  }
  StaticBitArray<N_STARS> computeXSyndrome() const {
    StaticBitArray<N_STARS> syn;
    for (int i = 0; i < N_STARS; ++i)
      if (zErr_.masked_parity(starMasks_[i]))
        syn.set(i);
    return syn;
  }

  // ------------ Geometry helpers ------------
  int numHRows() const { return L; }
  int numHCols() const { return H_cols_; }
  int numVRows() const { return V_rows_; }
  int numVCols() const { return L; }

  // horizontal coordinate to global index
  int qubitIdx_HG(int r, int c) const {
    assert(0 <= r && r < L);
    assert(0 <= c && c < L);
    return hIndex_max(r, c);
  }
  // vertical coordinate to global index
  int qubitIdx_VG(int r, int c) const {
    assert(0 <= r && r < L);
    assert(0 <= c && c < L);
    return vIndex_max(r, c);
  }
  // global index to horizontal coordinate
  std::pair<int, int> qubitIdx_GH(int hIdx) const {
    return std::make_pair(hIdx / L, hIdx % L);
  }
  // global index to vertical coordinate
  std::pair<int, int> qubitIdx_GV(int vIdx) const {
    assert(vIdx >= NH_TORIC && vIdx < NQ_MAX);
    int idx = vIdx - NH_TORIC;
    return std::make_pair(idx / L, idx % L);
  }

  // Flip by coordinate (silently a no-op in planar if edge is out-of-range)
  void flip_h_X(int r, int c) {
    if (periodic_ || c < L - 1)
      flipXErr(qubitIdx_HG(r, c));
  }
  void flip_h_Z(int r, int c) {
    if (periodic_ || c < L - 1)
      flipZErr(qubitIdx_HG(r, c));
  }
  void flip_v_X(int r, int c) {
    if (periodic_ || r < L - 1)
      flipXErr(qubitIdx_VG(r, c));
  }
  void flip_v_Z(int r, int c) {
    if (periodic_ || r < L - 1)
      flipZErr(qubitIdx_VG(r, c));
  }

  // Generate SVG representation of the surface code
  void generateSVG(svg::Document& doc) const {
    // Generate SVG code for the surface code representation
    constexpr int PADDING = 15;  // Padding around the grid
    constexpr int CELLSIZE = 40; // Size of each cell in the grid

    auto& path = doc.addPath();

    // horizontal lines
    for (int row = 0; row < L; ++row) {
      path.moveTo(PADDING, PADDING + row * CELLSIZE)
          .hTo(PADDING + L * CELLSIZE);
    }

    // vertical lines
    for (int col = 0; col < L; ++col) {
      path.moveTo(PADDING + col * CELLSIZE, PADDING)
          .vTo(PADDING + L * CELLSIZE);
    }

    // bottom and right edges for periodic boundary conditions
    {
      auto& path = doc.addPath()
                       .moveTo(PADDING + L * CELLSIZE, PADDING)
                       .vTo(PADDING + L * CELLSIZE)
                       .hTo(PADDING);
      if (periodic_)
        path.setStrokeDashArray("5,5");
    }

    constexpr float STROKE_WIDTH = 3.0f;
    // X error
    {
      auto& path = doc.addPath().setStrokeWidth(STROKE_WIDTH).setStroke("red");
      // horizontal edges
      for (int gIdx = 0; gIdx < Nh_used_; ++gIdx) {
        if (hasXErr(gIdx) && !hasZErr(gIdx)) {
          auto [r, c] = qubitIdx_GH(gIdx);
          path.moveTo(PADDING + c * CELLSIZE, PADDING + r * CELLSIZE)
              .hTo(PADDING + (c + 1) * CELLSIZE);
        }
      }
      // vertical edges
      for (int gIdx = NH_TORIC; gIdx < Nv_used_ + NH_TORIC; ++gIdx) {
        if (hasXErr(gIdx) && !hasZErr(gIdx)) {
          auto [r, c] = qubitIdx_GV(gIdx);
          path.moveTo(PADDING + c * CELLSIZE, PADDING + r * CELLSIZE)
              .vTo(PADDING + (r + 1) * CELLSIZE);
        }
      }
    }

    // Z error
    {
      auto& path = doc.addPath().setStrokeWidth(STROKE_WIDTH).setStroke("blue");
      // horizontal edges
      for (int gIdx = 0; gIdx < Nh_used_; ++gIdx) {
        if (hasZErr(gIdx) && !hasXErr(gIdx)) {
          auto [r, c] = qubitIdx_GH(gIdx);
          path.moveTo(PADDING + c * CELLSIZE, PADDING + r * CELLSIZE)
              .hTo(PADDING + (c + 1) * CELLSIZE);
        }
      }
      // vertical edges
      for (int gIdx = NH_TORIC; gIdx < Nv_used_ + NH_TORIC; ++gIdx) {
        if (hasZErr(gIdx) && !hasXErr(gIdx)) {
          auto [r, c] = qubitIdx_GV(gIdx);
          path.moveTo(PADDING + c * CELLSIZE, PADDING + r * CELLSIZE)
              .vTo(PADDING + (r + 1) * CELLSIZE);
        }
      }
    }

    // Y error
    {
      auto& path =
          doc.addPath().setStrokeWidth(STROKE_WIDTH).setStroke("green");
      // horizontal edges
      for (int gIdx = 0; gIdx < Nh_used_; ++gIdx) {
        if (hasXErr(gIdx) && hasZErr(gIdx)) {
          auto [r, c] = qubitIdx_GH(gIdx);
          path.moveTo(PADDING + c * CELLSIZE, PADDING + r * CELLSIZE)
              .hTo(PADDING + (c + 1) * CELLSIZE);
        }
      }
      // vertical edges
      for (int gIdx = NH_TORIC; gIdx < Nv_used_ + NH_TORIC; ++gIdx) {
        if (hasXErr(gIdx) && hasZErr(gIdx)) {
          auto [r, c] = qubitIdx_GV(gIdx);
          path.moveTo(PADDING + c * CELLSIZE, PADDING + r * CELLSIZE)
              .vTo(PADDING + (r + 1) * CELLSIZE);
        }
      }
    }

    // Syndrome
    {
      auto zSyndrome = computeZSyndrome();
      auto xSyndrome = computeXSyndrome();

      // X syndrome: red circles on vertices
      for (int i = 0; i < numPlaquettes(); ++i) {
        if (xSyndrome.test(i)) {
          auto [r, c] = qubitIdx_GH(i);
          doc.addCircle(PADDING + c * CELLSIZE,
                        PADDING + r * CELLSIZE,
                        CELLSIZE * 0.1f)
              .setFill("red");
        }
      }

      // Z syndrome: blue rectangles on vertices
      for (int i = 0; i < numStars(); ++i) {
        if (zSyndrome.test(i)) {
          auto [r, c] = qubitIdx_GH(i);
          doc.addRect(PADDING + (c + 0.1f) * CELLSIZE,
                      PADDING + (r + 0.1f) * CELLSIZE,
                      CELLSIZE * 0.8f,
                      CELLSIZE * 0.8f)
              .setFill("blue").setOpacity(0.5f);
        }
      }

    }

  } // generateSVG
};

} // namespace uf
