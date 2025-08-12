#ifndef UTILS_SMALLSTRING_H
#define UTILS_SMALLSTRING_H

#include <cassert>
#include <cstddef>
#include <ostream>
#include <string>
#include <string_view>

namespace utils {

// SmallString is designed to hold up to 15 char inline, and larger strings on
// the heap. Larger strings are stored in a std::string on the heap.
class SmallString {
private:
  // 15 bytes for inline storage + 1 byte for flags
  alignas(std::string*) std::byte data_[16];

  static constexpr int NumInline = 15;

  static constexpr std::byte ZERO_TAG{0x00};
  static constexpr std::byte HEAP_TAG{0xFF};

  std::string const* get_heap_string() const {
    return *reinterpret_cast<std::string* const*>(data_);
  }

  std::string* get_heap_string() {
    return *reinterpret_cast<std::string**>(data_);
  }

  void set_heap_string(std::string* str) {
    assert(get_flag() == HEAP_TAG);
    *reinterpret_cast<std::string**>(data_) = str;
  }

  std::byte get_flag() const { return data_[NumInline]; }

  void set_flag(std::byte flag) { data_[NumInline] = flag; }

public:
  SmallString() {
    set_flag(ZERO_TAG); // empty inline storage
  }

  SmallString(const char* str) : SmallString(std::string_view(str)) {}

  SmallString(const std::string& str) : SmallString(std::string_view(str)) {}

  SmallString(std::string_view sv) {
    if (sv.size() <= NumInline) {
      set_flag(static_cast<std::byte>(sv.size()));
      std::memcpy(data_, sv.data(), sv.size());
      return;
    }
    // dynamic storage
    set_flag(HEAP_TAG);
    set_heap_string(new std::string(sv));
  }

  // Copy constructor
  SmallString(const SmallString& other) {
    if (other.get_flag() == HEAP_TAG) {
      // copy heap string
      set_flag(HEAP_TAG);
      set_heap_string(new std::string(*other.get_heap_string()));
    } else {
      set_flag(other.get_flag());
      std::memcpy(data_, other.data_, NumInline);
    }
  }

  // Move constructor
  SmallString(SmallString&& other) noexcept {
    if (other.get_flag() == HEAP_TAG) {
      // steal the heap string
      set_flag(HEAP_TAG);
      set_heap_string(other.get_heap_string());
      other.set_flag(ZERO_TAG); // reset moved-from object
    } else {
      // this will also copy the flag byte
      std::memcpy(data_, other.data_, 16);
    }
  }

  // Copy assignment operator
  SmallString& operator=(const SmallString& other) {
    if (this == &other)
      return *this;
    this->~SmallString();
    new (this) SmallString(other);
    return *this;
  }

  // Move assignment operator
  SmallString& operator=(SmallString&& other) noexcept {
    if (this == &other)
      return *this;
    this->~SmallString();
    new (this) SmallString(std::move(other));
    return *this;
  }

  // Destructor
  ~SmallString() {
    if (get_flag() == HEAP_TAG)
      delete get_heap_string();
  }

  std::string_view view() const {
    if (get_flag() == HEAP_TAG)
      return std::string_view(*get_heap_string());
    return std::string_view(reinterpret_cast<const char*>(data_),
                            static_cast<size_t>(get_flag()));
  }

  bool empty() const {
    return get_flag() == ZERO_TAG ||
           (get_flag() == HEAP_TAG && get_heap_string()->empty());
  }

  size_t size() const {
    if (get_flag() == HEAP_TAG)
      return get_heap_string()->size();
    return static_cast<size_t>(get_flag());
  }

  friend std::ostream& operator<<(std::ostream& os, const SmallString& str) {
    auto view = str.view();
    os.write(view.data(), view.size());
    return os;
  }

}; // SmallString

inline bool operator==(const SmallString& lhs, const SmallString& rhs) {
  return lhs.view() == rhs.view();
}

inline bool operator<(const SmallString& lhs, const SmallString& rhs) {
  return lhs.view() < rhs.view();
}

} // namespace utils

#endif // UTILS_SMALLSTRING_H