#ifndef SVG_TEXT_H
#define SVG_TEXT_H

#include "SVG/SVGObject.h"

namespace svg {

class Text : public SVGObject<Text> {
  float x_, y_;                // Position of the text
  utils::SmallString content_; // Text content

public:
  Text(float x, float y, utils::SmallString content)
      : SVGObject<Text>(TextKind), x_(x), y_(y), content_(content) {}

  static constexpr const char* element_name = "text";

  static bool classof(const SVGObjectBase* base) {
    return base->getKind() == TextKind;
  }

  Text& setFontFamily(utils::SmallString family) {
    return setStyle("font-family", std::move(family));
  }
  Text& setFontSize(float size) {
    return setStyle("font-size", to_ss(size));
  }
  Text& setFontWeight(utils::SmallString weight) {
    return setStyle("font-weight", std::move(weight));
  }
  Text& setFill(utils::SmallString color) {
    return setStyle("fill", std::move(color));
  }
  Text& setTextAnchor(utils::SmallString anchor) {
    return setStyle("text-anchor", std::move(anchor));
  }
  Text& setDominantBaseline(utils::SmallString baseline) {
    return setStyle("dominant-baseline", std::move(baseline));
  }
  Text& setOpacity(float op) { return setStyle("opacity", to_ss(op)); }

  void render_impl(std::ostream& os) const {
    write_open(os);
    attr(os, "x", x_);
    attr(os, "y", y_);
    write_styles(os);
    os << '>' << content_ << "</" << element_name << ">\n";
  }
};

} // namespace svg

#endif // SVG_TEXT_H