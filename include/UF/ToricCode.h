#ifndef UF_TORICCODE_H
#define UF_TORICCODE_H

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <random>

#include "SVG/Document.h"

namespace uf {

// StaticBitArray<Nbits>: fixed-size bitset backed by 64-bit words
template <int Nbits> class StaticBitArray {
  static_assert(Nbits > 0, "Nbits must be positive");
  static constexpr int NBitsPerWord = 64;
  static constexpr int NWords = (Nbits + NBitsPerWord - 1) / NBitsPerWord;
  std::array<uint64_t, NWords> w_{}; // zero-initialized

  // Mark unused bits
  static constexpr uint64_t tail_mask() {
    if constexpr ((Nbits & 63) == 0)
      return ~0ULL;
    else
      return (~0ULL) >> (64 - (Nbits & 63));
  }

public:
  static consteval int num_bits() { return Nbits; }
  static consteval int num_words() { return NWords; }

  void clear() { w_.fill(0ULL); }

  void set_all() {
    w_.fill(~0ULL);
    w_.back() &= tail_mask();
  }

  // Set bit i to 1
  void set(int i) { w_[i >> 6] |= (1ULL << (i & 63)); }
  // Set bit i to 0
  void unset(int i) { w_[i >> 6] &= ~(1ULL << (i & 63)); }
  // Toggle bit i
  void flip(int i) { w_[i >> 6] ^= (1ULL << (i & 63)); }
  // Check if bit i is 1
  bool test(int i) const { return (w_[i >> 6] >> (i & 63)) & 1ULL; }

  void set_word(int i, uint64_t value) {
    assert(0 <= i && i < NWords);
    w_[i] = value;
  }

  uint64_t word(int i) const { return w_[i]; }
  uint64_t* words() { return w_.data(); }
  const uint64_t* words() const { return w_.data(); }

  void xor_with(const StaticBitArray& other) {
    for (int i = 0; i < NWords; ++i)
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
    for (int i = 0; i < NWords; ++i)
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

template <int L> using ToricCodeError = StaticBitArray<2 * L * L>;

template <int L> using ToricCodeSyndrome = StaticBitArray<L * L>;

// ToricCode<L>: Toric code with L x L lattice
template <int L> class ToricCode {
  static_assert(L >= 2, "L must be >= 2");

  // L * L
  static constexpr int LL = L * L;

  // Number of qubits: equals to 2 * L * L
  static constexpr int NQubits = 2 * LL;

public:
  // Error bit arrays (length 2*L*L)
  using Error = ToricCodeError<L>;
  // Syndrome bit arrays (length L*L)
  using Syndrome = ToricCodeSyndrome<L>;

  // Error bitfields (zero-initialized)
  Error xErr{};
  Error zErr{};
  Error erasureErr{};

  // ---------------- Indexing helpers ----------------
  // Qubits (edges) are divided into two groups:
  // horizontal edges (H) and vertical edges (V).
  // For toric code, H edges and V edges both have L rows and L columns. We
  // enumerate them row-wise.

  // There is also a global qubit index (G), given by enumerating H edges
  // followed by V edges.

  // horizontal coordinate to global index
  // both r and c must be in [0, L)
  static int qubitIdx_HG(int r, int c) {
    assert(0 <= r && r < L);
    assert(0 <= c && c < L);
    return r * L + c;
  }

  // vertical coordinate to global index
  // both r and c must be in [0, L)
  static int qubitIdx_VG(int r, int c) {
    assert(0 <= r && r < L);
    assert(0 <= c && c < L);
    return LL + r * L + c;
  }

  // global index to horizontal coordinate
  // idx must be in [0, LL)
  static std::pair<int, int> qubitIdx_GH(int idx) {
    assert(0 <= idx && idx < LL);
    return std::make_pair(idx / L, idx % L);
  }

  // global index to vertical coordinate
  // idx must be in [LL, 2*LL)
  static std::pair<int, int> qubitIdx_GV(int idx) {
    assert(LL <= idx && idx < 2 * LL);
    return std::make_pair((idx - LL) / L, (idx - LL) % L);
  }

public:
  ToricCode() {}

  // Grid size is L
  static constexpr int gridSize() { return L; }

  // ------------ Error manipulation ------------

  template <class URNG> void injectPauliNoise(float pX, float pZ, URNG& rng) {
    assert(0.0 <= pX && pX <= 1.0);
    assert(0.0 <= pZ && pZ <= 1.0);
    if (pX > 0.0f) {
      std::bernoulli_distribution b(pX);
      for (int q = 0; q < NQubits; ++q)
        if (b(rng))
          xErr.flip(q);
    }
    if (pZ > 0.0f) {
      std::bernoulli_distribution b(pZ);
      for (int q = 0; q < NQubits; ++q)
        if (b(rng))
          zErr.flip(q);
    }
  }

  template <class URNG> void injectErasureNoise(float p, URNG& rng) {
    assert(0.0 <= p && p <= 1.0);
    if (p > 0.0f) {
      std::uniform_real_distribution<float> d(0.0f, 1.0f);
      for (int q = 0; q < NQubits; ++q) {
        auto r = d(rng);
        if (r < 0.25f * p) {
          // Erasure with I
          erasureErr.flip(q);
        } else if (r < 0.5f * p) {
          // Erasure with X
          xErr.flip(q);
          erasureErr.flip(q);
        } else if (r < 0.75f * p) {
          // Erasure with Z
          zErr.flip(q);
          erasureErr.flip(q);
        } else if (r < p) {
          // Erasure with Y
          xErr.flip(q);
          zErr.flip(q);
          erasureErr.flip(q);
        }
      }
    }
  }

  // ------------ Syndromes ------------
  Syndrome computeZSyndrome() const {
    Syndrome syn;
    // Z syndrome (plaquettes) is computed by checking X errors on surrounding
    // edges

    // Remark: L is a compile-time constant. Compilers will optimize away all
    // constants.

    // Plaquette (face) (r, c) is surrounded by:
    // top edge    (H, r,     c    ),
    // bottom edge (H, r + 1, c    ),
    // left edge   (V, r,     c    ),
    // right edge  (V, r,     c + 1).

    int parity;
    for (int r = 0; r < L - 1; ++r) {
      for (int c = 0; c < L - 1; ++c) {
        parity = 0;
        parity ^= xErr.test(qubitIdx_HG(r, c));     // top (H, r, c)
        parity ^= xErr.test(qubitIdx_HG(r + 1, c)); // bottom (H, r + 1, c)
        parity ^= xErr.test(qubitIdx_VG(r, c));     // left (V, r, c)
        parity ^= xErr.test(qubitIdx_VG(r, c + 1)); // right (V, r, c + 1)
        if (parity)
          syn.set(r * L + c);
      }
      // c = L - 1 case
      parity = 0;
      parity ^= xErr.test(qubitIdx_HG(r, L - 1));     // top (H, r, c)
      parity ^= xErr.test(qubitIdx_HG(r + 1, L - 1)); // bottom (H, r + 1, c)
      parity ^= xErr.test(qubitIdx_VG(r, L - 1));     // left (V, r, c)
      parity ^= xErr.test(qubitIdx_VG(r, 0));         // right (V, r, c + 1)
      if (parity)
        syn.set(r * L + (L - 1));
    }

    // r = L - 1 case
    for (int c = 0; c < L - 1; ++c) {
      parity = 0;
      parity ^= xErr.test(qubitIdx_HG(L - 1, c));     // top (H, r, c)
      parity ^= xErr.test(qubitIdx_HG(L - 1, c + 1)); // bottom (H, r + 1, c)
      parity ^= xErr.test(qubitIdx_VG(L - 1, c));     // left (V, r, c)
      parity ^= xErr.test(qubitIdx_VG(L - 1, c + 1)); // right (V, r, c + 1)
      if (parity)
        syn.set((L - 1) * L + c);
    }
    // r = L - 1, c = L - 1 case
    parity = 0;
    parity ^= xErr.test(qubitIdx_HG(L - 1, L - 1)); // top (H, r, c)
    parity ^= xErr.test(qubitIdx_HG(0, L - 1));     // bottom (H, r + 1, c)
    parity ^= xErr.test(qubitIdx_VG(L - 1, L - 1)); // left (V, r, c)
    parity ^= xErr.test(qubitIdx_VG(L - 1, 0));     // right (V, r, c + 1)
    if (parity)
      syn.set((L - 1) * L + (L - 1)); // LL - 1

    return syn;
  }

  Syndrome computeXSyndrome() const {
    Syndrome syn;
    // X syndrome (stars) is computed by checking Z errors on connected edges

    // Remark: L is a compile-time constant. Compilers will optimize away all
    // constants.

    // Star (vertex) (r, c) is connected to:
    // up edge    (V, r - 1, c    ),
    // down edge  (V, r,     c    ),
    // left edge  (H, r,     c - 1),
    // right edge (H, r,     c    ).

    // r = 0, c = 0 case
    int parity;
    parity = 0;
    parity ^= zErr.test(qubitIdx_VG(L - 1, 0)); // top (V, r - 1, c)
    parity ^= zErr.test(qubitIdx_VG(0, 0));     // bottom (V, r, c)
    parity ^= zErr.test(qubitIdx_HG(0, L - 1)); // left (H, r, c - 1)
    parity ^= zErr.test(qubitIdx_HG(0, 0));     // right (H, r, c)
    if (parity)
      syn.set(0);

    // r = 0 case
    for (int c = 1; c < L; ++c) {
      parity = 0;
      parity ^= zErr.test(qubitIdx_VG(L - 1, c)); // top (V, r - 1, c)
      parity ^= zErr.test(qubitIdx_VG(0, c));     // bottom (V, r, c)
      parity ^= zErr.test(qubitIdx_HG(0, c - 1)); // left (H, r, c - 1)
      parity ^= zErr.test(qubitIdx_HG(0, c));     // right (H, r, c)
      if (parity)
        syn.set(c);
    }

    for (int r = 1; r < L; ++r) {
      for (int c = 1; c < L; ++c) {
        parity = 0;
        parity ^= zErr.test(qubitIdx_VG(r - 1, c)); // top (V, r - 1, c)
        parity ^= zErr.test(qubitIdx_VG(r, c));     // bottom (V, r, c)
        parity ^= zErr.test(qubitIdx_HG(r, c - 1)); // left (H, r, c - 1)
        parity ^= zErr.test(qubitIdx_HG(r, c));     // right (H, r, c)
        if (parity)
          syn.set(r * L + c);
      }
    }
    return syn;
  }

};

} // namespace uf

#endif // UF_TORICCODE_H