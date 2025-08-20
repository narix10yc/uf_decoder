#ifndef UF_UFDECODER3D_H
#define UF_UFDECODER3D_H

#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <memory>
#include <queue>
#include <vector>

namespace uf {

struct EdgeCoordinate3D {
  enum Group {
    Horizontal, // H edges
    Vertical,   // V edges
    Temporal    // T edges
  };
  Group group;
  int t; // time coordinate
  int r; // row coordinate
  int c; // column coordinate

  EdgeCoordinate3D(int idx, int L, int T);
};

std::ostream& operator<<(std::ostream& os, const EdgeCoordinate3D& coord);

static inline int edgeIdx_H(int t, int r, int c, int L) {
  assert(0 <= r && r < L);
  assert(0 <= c && c < L);
  return (3 * L * L) * t + (L * r) + c;
}

static inline int edgeIdx_V(int t, int r, int c, int L) {
  assert(0 <= r && r < L);
  assert(0 <= c && c < L);
  return (3 * L * L) * t + (L * r) + c + (L * L);
}

static inline int edgeIdx_T(int t, int r, int c, int L) {
  assert(0 <= r && r < L);
  assert(0 <= c && c < L);
  return (3 * L * L) * t + (L * r) + c + (2 * L * L);
}

class BitArray {
  size_t nbits_;
  size_t nwords_;
  uint64_t* data_;

  static inline int ctz64(uint64_t x) {
    // x != 0 is guaranteed by callers
    return __builtin_ctzll(x);
  }

  uint64_t last_word_mask() const {
    const unsigned r = unsigned(nbits_ & 63);
    return r ? ((r == 64 ? ~0ULL : ((1ULL << r) - 1ULL))) : ~0ULL;
  }

public:
  // Zero-initialized BitArray
  BitArray(size_t nbits)
      : nbits_(nbits), nwords_((nbits + 63) / 64),
        data_(new uint64_t[nwords_]()) {
    unset_all();
  }

  ~BitArray() { delete[] data_; }

  BitArray(const BitArray& other)
      : nbits_(other.nbits_), nwords_(other.nwords_),
        data_(new uint64_t[nwords_]) {
    std::memcpy(data_, other.data_, sizeInBytes());
  }

  BitArray& operator=(const BitArray& other) {
    if (this == &other)
      return *this;
    delete[] data_;
    nbits_ = other.nbits_;
    nwords_ = other.nwords_;
    data_ = new uint64_t[nwords_];
    std::memcpy(data_, other.data_, sizeInBytes());
    return *this;
  }

  BitArray(BitArray&& other) noexcept
      : nbits_(other.nbits_), nwords_(other.nwords_), data_(other.data_) {
    other.data_ = nullptr;
  }

  BitArray& operator=(BitArray&& other) noexcept {
    if (this == &other)
      return *this;
    delete[] data_;
    nbits_ = other.nbits_;
    nwords_ = other.nwords_;
    data_ = other.data_;
    other.data_ = nullptr;
    return *this;
  }

  void set(size_t index) {
    assert(index < nbits_);
    data_[index / 64] |= (1ULL << (index % 64));
  }

  void unset(size_t index) {
    assert(index < nbits_);
    data_[index / 64] &= ~(1ULL << (index % 64));
  }

  void flip(size_t index) {
    assert(index < nbits_);
    data_[index / 64] ^= (1ULL << (index % 64));
  }

  void set_word(size_t index, uint64_t value) {
    assert(index < nwords_);
    data_[index] = value;
  }

  bool test(size_t index) const {
    assert(index < nbits_);
    return data_[index / 64] & (1ULL << (index % 64));
  }

  void unset_all() { std::memset(data_, 0, sizeInBytes()); }

  size_t nbits() const { return nbits_; }
  size_t nwords() const { return nwords_; }
  size_t sizeInBytes() const { return nwords_ * sizeof(uint64_t); }

  BitArray& operator|=(const BitArray& other) {
    assert(nbits_ == other.nbits_);
    for (size_t i = 0; i < nwords_; ++i)
      data_[i] |= other.data_[i];
    return *this;
  }

  BitArray& operator&=(const BitArray& other) {
    assert(nbits_ == other.nbits_);
    for (size_t i = 0; i < nwords_; ++i)
      data_[i] &= other.data_[i];
    return *this;
  }

  BitArray& operator^=(const BitArray& other) {
    assert(nbits_ == other.nbits_);
    for (size_t i = 0; i < nwords_; ++i)
      data_[i] ^= other.data_[i];
    return *this;
  }

  /* ---- Iterator ---- */

  class OneIterator {
    const BitArray* ba_ = nullptr;
    size_t word_idx_ = 0; // which 64-bit word we are in
    uint64_t word_ = 0;   // remaining set bits in current word (masked)
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = size_t;
    using difference_type = std::ptrdiff_t;
    using reference = size_t; // returns by value
    using pointer = void;

    OneIterator() = default;

    // begin ctor
    explicit OneIterator(const BitArray* ba) : ba_(ba), word_idx_(0), word_(0) {
      if (!ba_)
        return;
      if (ba_->nwords_ == 0) {
        word_idx_ = 0;
        word_ = 0;
        return;
      }
      // load first non-empty (masked) word
      advance_to_next_nonempty_from(word_idx_);
    }

    // end ctor
    static OneIterator make_end(const BitArray* ba) {
      OneIterator it;
      it.ba_ = ba;
      it.word_idx_ = ba ? ba->nwords_ : 0;
      it.word_ = 0;
      return it;
    }

    // dereference -> index of current set bit
    size_t operator*() const {
      // valid only if not at end
      int bit = ctz64(word_);
      return word_idx_ * 64 + size_t(bit);
    }

    OneIterator& operator++() {
      // clear current lowest set bit
      word_ &= (word_ - 1);
      if (word_ == 0) {
        // move to next non-empty word
        ++word_idx_;
        advance_to_next_nonempty_from(word_idx_);
      }
      return *this;
    }

    OneIterator operator++(int) {
      OneIterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const OneIterator& rhs) const {
      // same container and same position
      return ba_ == rhs.ba_ && word_idx_ == rhs.word_idx_ && word_ == rhs.word_;
    }
    bool operator!=(const OneIterator& rhs) const { return !(*this == rhs); }

  private:
    void advance_to_next_nonempty_from(size_t start) {
      // advance word_idx_ to first word with any 1-bits; set word_ accordingly
      const size_t nwords = ba_->nwords_;
      const uint64_t last_mask = ba_->last_word_mask();
      word_idx_ = start;
      while (word_idx_ < nwords) {
        uint64_t w = ba_->data_[word_idx_];
        if (word_idx_ == nwords - 1)
          w &= last_mask; // mask off unused tail bits
        if (w) {
          word_ = w;
          return;
        }
        ++word_idx_;
      }
      // reached end
      word_ = 0;
    }
  };

  // begin/end for one-bits
  OneIterator one_indices() const { return OneIterator(this); }
  OneIterator one_indices_end() const { return OneIterator::make_end(this); }

  // range-for support:
  struct OneRange {
    const BitArray* ba;
    OneIterator begin() const { return ba->one_indices(); }
    OneIterator end() const { return BitArray::OneIterator::make_end(ba); }
  };
  OneRange ones() const { return OneRange{this}; }

}; // class BitArray

std::ostream& operator<<(std::ostream& os, const BitArray& arr);

class DisjointSetUnion {
  std::unique_ptr<unsigned[]> parents;
  std::unique_ptr<uint8_t[]> ranks;
  unsigned size_;

public:
  DisjointSetUnion(unsigned size)
      : parents(std::make_unique<unsigned[]>(size)),
        ranks(std::make_unique<uint8_t[]>(size)), size_(size) {
    for (unsigned i = 0; i < size; ++i)
      parents[i] = i;
    std::memset(ranks.get(), 0, size * sizeof(uint8_t));
  }

  unsigned find(unsigned x) {
    assert(x < size_);
    while (parents[x] != x) {
      parents[x] = parents[parents[x]]; // path compression
      x = parents[x];
    }
    return x;
  }

  bool unite(unsigned a, unsigned b) {
    assert(a < size_ && b < size_);
    a = find(a);
    b = find(b);
    if (a == b)
      return false; // already in the same set
    if (ranks[a] < ranks[b])
      std::swap(a, b);
    parents[b] = a;
    if (ranks[a] == ranks[b])
      ++ranks[a];
    return true;
  }
};

static inline int trc(int t, int r, int c, int L) {
  return t * L * L + r * L + c;
}

void clustering(BitArray& erasure, const BitArray& syndromes, int L, int T);

void peelingDecode(BitArray& correction,
                   const BitArray& erasure,
                   BitArray syndromes,
                   int L,
                   int T);

} // namespace uf

#endif // UF_UFDECODER3D_H