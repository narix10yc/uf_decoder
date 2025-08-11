#ifndef SVG_STYLE_H
#define SVG_STYLE_H

#include "Utils/SmallString.h"
#include <iostream>

namespace svg {

struct Style {
  utils::SmallString stroke = "#000";
  double stroke_width = 1.0;
  utils::SmallString fill = "none";
  utils::SmallString stroke_dasharray; // e.g. "5,3"
  utils::SmallString stroke_linecap;   // "butt", "round", "square"
  utils::SmallString stroke_linejoin;  // "miter", "round", "bevel"
  double opacity = 1.0;
};

} // namespace svg

std::ostream& operator<<(std::ostream& os, const svg::Style& style) {
  os << "style=\"";
  os << "stroke:" << style.stroke << ";";
  os << "stroke-width:" << style.stroke_width << ";";
  os << "fill:" << style.fill << ";";
  if (!style.stroke_dasharray.empty())
    os << "stroke-dasharray:" << style.stroke_dasharray << ";";
  if (!style.stroke_linecap.empty())
    os << "stroke-linecap:" << style.stroke_linecap << ";";
  if (!style.stroke_linejoin.empty())
    os << "stroke-linejoin:" << style.stroke_linejoin << ";";
  if (style.opacity != 1.0)
    os << "opacity:" << style.opacity << ";";
  os << "\"";
  return os;
}

#endif // SVG_STYLE_H