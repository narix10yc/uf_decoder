// LLVM-style RTTI helpers
// As this project does not use LLVM, we implement a minimal version of the
// casting utilities.
// This file mimics llvm/Support/Casting.h
// Functions are under the `utils` namespace.

#ifndef UTILS_CASTING_H
#define UTILS_CASTING_H

namespace utils {

template <typename To, typename From> inline To* dyn_cast(From* base) {
  if (To::classof(base)) {
    return static_cast<To*>(base);
  }
  return nullptr;
}

template <typename To, typename From>
inline const To* dyn_cast(const From* base) {
  if (To::classof(base)) {
    return static_cast<const To*>(base);
  }
  return nullptr;
}

template <typename To, typename From> inline bool isa(const From* base) {
  return To::classof(base);
}

} // namespace utils

#endif // UTILS_CASTING_H