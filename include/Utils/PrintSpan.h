#ifndef UTILS_PRINTSPAN_H
#define UTILS_PRINTSPAN_H

#include <iostream>
#include <span>

namespace utils {

template <typename T> void printSpan(std::ostream& os, std::span<T> span) {
  if (span.empty()) {
    os << "[]";
    return;
  }
  os << "[" << span[0];
  for (size_t i = 1; i < span.size(); ++i)
    os << ", " << span[i];
  os << "]";
}

} // namespace utils

#endif // UTILS_PRINTSPAN_H