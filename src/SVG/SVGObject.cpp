#include "SVG/SVGObject.h"

utils::SmallString svg::to_ss(float v) {
  char buf[64];
  auto [p, ec] = std::to_chars(buf,
                               buf + sizeof(buf),
                               v,
                               std::chars_format::general,
                               9); // ~float precision
  if (ec != std::errc{}) {
    // very rare fallback
    int n = std::snprintf(buf, sizeof(buf), "%.9g", static_cast<double>(v));
    p = buf + std::max(0, n);
  }
  // trim trailing zeros and '.' (no locale)
  char* first = buf;
  char* last = p;
  char* dot = nullptr;
  for (char* q = first; q != last; ++q)
    if (*q == '.') {
      dot = q;
      break;
    }
  if (dot) {
    while (last > dot + 1 && *(last - 1) == '0')
      --last;
    if (last > dot && *(last - 1) == '.')
      --last;
  }
  if ((last - first) >= 2 && first[0] == '-' && first[1] == '0') { // -0 -> 0
    bool allZero = true;
    for (const char* q = first + 1; q != last; ++q)
      if (*q != '0' && *q != '.') {
        allZero = false;
        break;
      }
    if (allZero) {
      first[0] = '0';
      last = first + 1;
    }
  }
  return utils::SmallString(
      std::string_view(first, static_cast<size_t>(last - first)));
}