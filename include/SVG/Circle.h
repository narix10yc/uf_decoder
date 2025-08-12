#ifndef SVG_CIRCLE_H
#define SVG_CIRCLE_H

#include "SVG/SVGObject.h"

namespace svg {

class Circle : public SVGObject<Circle> {
  float cx_, cy_, r_;

public:
  Circle(float cx, float cy, float r)
      : SVGObject<Circle>(CircleKind), cx_(cx), cy_(cy), r_(r) {}

  static constexpr const char* element_name = "circle";

  static bool classof(const SVGObjectBase* base) {
    return base->getKind() == CircleKind;
  }

  Circle& setFill(utils::SmallString color) {
    return setStyle("fill", std::move(color));
  }

  Circle& setStroke(utils::SmallString color) {
    return setStyle("stroke", std::move(color));
  }

  Circle& setStrokeWidth(float w) {
    return setStyle("stroke-width", w);
  }

  void render_impl(std::ostream& os) const {
    write_open(os);
    attr(os, "cx", cx_);
    attr(os, "cy", cy_);
    attr(os, "r", r_);
    write_styles(os);
    write_close(os);
  }
};

} // namespace svg

#endif // SVG_CIRCLE_H