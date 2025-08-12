#ifndef SVG_LINE_H
#define SVG_LINE_H

#include "SVG/SVGObject.h"

namespace svg {

class Line : public SVGObject<Line> {
  float x1_, y1_, x2_, y2_;

public:
  Line(float x1, float y1, float x2, float y2)
      : SVGObject<Line>(LineKind), x1_(x1), y1_(y1), x2_(x2), y2_(y2) {}

  static constexpr const char* element_name = "line";

  static bool classof(const SVGObjectBase* base) {
    return base->getKind() == LineKind;
  }

  Line& setStroke(utils::SmallString color) {
    return setStyle("stroke", color);
  }
  Line& setStrokeWidth(float w) {
    return setStyle("stroke-width", to_ss(w));
  }
  Line& setDashArray(utils::SmallString pattern) {
    return setStyle("stroke-dasharray", pattern);
  }
  Line& setLineCap(utils::SmallString cap) {
    return setStyle("stroke-linecap", cap);
  }
  Line& setLineJoin(utils::SmallString join) {
    return setStyle("stroke-linejoin", join);
  }
  Line& setOpacity(float op) { return setStyle("opacity", to_ss(op)); }

  void render_impl(std::ostream& os) const {
    write_open(os);
    attr(os, "x1", x1_);
    attr(os, "y1", y1_);
    attr(os, "x2", x2_);
    attr(os, "y2", y2_);
    write_styles(os);
    write_close(os);
  }
};

} // namespace svg

#endif // SVG_LINE_H