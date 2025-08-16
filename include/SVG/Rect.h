#ifndef SVG_RECT_H
#define SVG_RECT_H

#include "SVG/SVGObject.h"

namespace svg {

class Rect : public SVGObject<Rect> {
  float x_, y_, width_, height_;
  float rx_{0.f}, ry_{0.f}; // optional rounded corners
public:
  /// @brief Construct a new Rect object
  /// @param rx The x-axis radius for rounded corners
  /// @param ry The y-axis radius for rounded corners
  Rect(float x,
       float y,
       float width,
       float height,
       float rx = 0.f,
       float ry = 0.f)
      : SVGObject<Rect>(RectKind), x_(x), y_(y), width_(width), height_(height),
        rx_(rx), ry_(ry) {}

  static constexpr const char* element_name = "rect";

  static bool classof(const SVGObjectBase* base) {
    return base->getKind() == RectKind;
  }

  Rect& setFill(utils::SmallString color) {
    return setStyle("fill", std::move(color));
  }
  Rect& setStroke(utils::SmallString color) {
    return setStyle("stroke", std::move(color));
  }
  Rect& setStrokeWidth(float w) {
    return setStyle("stroke-width", to_ss(w));
  }
  Rect& setDashArray(utils::SmallString pattern) {
    return setStyle("stroke-dasharray", std::move(pattern));
  }
  // miter/round/bevel
  Rect& setLineJoin(utils::SmallString v) {
    return setStyle("stroke-linejoin", std::move(v));
  }
  // butt/round/square
  Rect& setLineCap(utils::SmallString v) {
    return setStyle("stroke-linecap", std::move(v));
  }
  Rect& setOpacity(float o) { return setStyle("opacity", to_ss(o)); }

  void render_impl(std::ostream& os) const {
    write_open(os);
    attr(os, "x", x_);
    attr(os, "y", y_);
    attr(os, "width", width_);
    attr(os, "height", height_);
    if (rx_ > 0.0f || ry_ > 0.0f) {
      attr(os, "rx", rx_);
      attr(os, "ry", ry_);
    }
    write_styles(os);
    write_close(os);
  }

  Bounds computeBounds() const override {
    Bounds b;
    b.include(x_, y_);
    b.include(x_ + width_, y_ + height_);
    return b;
  }
};

} // namespace svg

#endif // SVG_RECT_H