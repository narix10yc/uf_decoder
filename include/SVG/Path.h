#ifndef SVG_PATH_H
#define SVG_PATH_H

#include "SVG/SVGObject.h"

namespace svg {

class Path : public SVGObject<Path> {
  std::string d_; // path data

  // append helpers
  void sp() { d_.push_back(' '); }
  void appendNumber(float v) {
    auto ss = to_ss(v);
    d_.append(ss.view().data(), ss.view().size());
  }
  void appendCmd(char c) { sp(); d_.push_back(c); }

public:
  Path() : SVGObject<Path>(PathKind), d_() {
    setStyle("fill", "none");
    setStyle("stroke", "#888");
  }

  static constexpr const char* element_name = "path";

  static bool classof(const SVGObjectBase* base) {
    return base->getKind() == PathKind;
  }

  // Geometry commands (absolute)
  Path& moveTo(float x, float y) {
    appendCmd('M');
    sp();
    appendNumber(x);
    sp();
    appendNumber(y);
    return *this;
  }

  Path& lineTo(float x, float y) {
    appendCmd('L');
    sp();
    appendNumber(x);
    sp();
    appendNumber(y);
    return *this;
  }

  Path& hTo(float x) {
    appendCmd('H');
    sp();
    appendNumber(x);
    return *this;
  }

  Path& vTo(float y) {
    appendCmd('V');
    sp();
    appendNumber(y);
    return *this;
  }

  Path& closePath() {
    appendCmd('Z');
    return *this;
  }

  // Common styling sugar
  Path& setStroke(utils::SmallString color) {
    return setStyle("stroke", std::move(color));
  }
  Path& setStrokeWidth(float w) { return setStyle("stroke-width", w); }
  Path& setStrokeDashArray(utils::SmallString dasharray) {
    return setStyle("stroke-dasharray", std::move(dasharray));
  }
  Path& setLineCap(utils::SmallString cap) {
    return setStyle("stroke-linecap", std::move(cap));
  }
  Path& setLineJoin(utils::SmallString join) {
    return setStyle("stroke-linejoin", std::move(join));
  }
  Path& setOpacity(float o) { return setStyle("opacity", o); }

  void render_impl(std::ostream& os) const {
    // materialize 'd', remove leading space if any
    std::string_view d_view;
    if (!d_.empty() && d_[0] == ' ')
      d_view = std::string_view(d_.data() + 1, d_.size() - 1);
    else
      d_view = std::string_view(d_.data(), d_.size());
    const_cast<Path*>(this)->styles_["d"] = utils::SmallString(d_view);

    this->write_open(os);
    this->write_styles(os);
    this->write_close(os);
  }
}; // class Path

} // namespace svg

#endif // SVG_PATH_H