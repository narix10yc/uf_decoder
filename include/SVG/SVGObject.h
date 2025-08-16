#ifndef SVG_SVGOBJECT_H
#define SVG_SVGOBJECT_H

#include "Utils/SmallString.h"
#include <map>

namespace svg {

utils::SmallString to_ss(float v);

struct Bounds {
  float xmin = std::numeric_limits<float>::infinity();
  float ymin = std::numeric_limits<float>::infinity();
  float xmax = -std::numeric_limits<float>::infinity();
  float ymax = -std::numeric_limits<float>::infinity();

  void include_x(float x) {
    xmin = std::min(xmin, x);
    xmax = std::max(xmax, x);
  }

  void include_y(float y) {
    ymin = std::min(ymin, y);
    ymax = std::max(ymax, y);
  }
  
  void include(float x, float y) {
    include_x(x);
    include_y(y);
  }

  void include(const Bounds& b) {
    include(b.xmin, b.ymin);
    include(b.xmax, b.ymax);
  }

  void pad(float p) {
    xmin -= p;
    ymin -= p;
    xmax += p;
    ymax += p;
  }
  
  bool valid() const { return xmin <= xmax && ymin <= ymax; }
  float width() const { return valid() ? (xmax - xmin) : 0.f; }
  float height() const { return valid() ? (ymax - ymin) : 0.f; }
};

class SVGObjectBase {
protected:
  enum Kind { LineKind, PathKind, RectKind, CircleKind, TextKind, EndKind_ };
  Kind kind_;

public:
  explicit SVGObjectBase(Kind kind) : kind_(kind) {}

  Kind getKind() const { return kind_; }

  virtual ~SVGObjectBase() = default;
  virtual void render(std::ostream& os) const = 0;
  virtual Bounds computeBounds() const = 0;

  // ----- Helpers -----
  static void attr(std::ostream& os, const char* key, float v) {
    os << ' ' << key << "=\"" << v << '"';
  }

  static void
  attr(std::ostream& os, const char* key, const utils::SmallString& v) {
    os << ' ' << key << "=\"" << v << '"';
  }
};

template <typename Derived> class SVGObject : public SVGObjectBase {
protected:
  std::map<utils::SmallString, utils::SmallString> styles_;

public:
  SVGObject(Kind kind) : SVGObjectBase(kind), styles_() {}

  Derived& setStyle(utils::SmallString key, utils::SmallString value) {
    styles_[key] = value;
    return static_cast<Derived&>(*this);
  }

  Derived& setStyle(utils::SmallString key, float value) {
    styles_[key] = to_ss(value);
    return static_cast<Derived&>(*this);
  }

  // Writes the opening bracket and the element name
  void write_open(std::ostream& os) const {
    os << '<' << Derived::element_name;
  }

  // Writes the style
  void write_styles(std::ostream& os) const {
    for (const auto& [key, value] : styles_)
      os << ' ' << key << "=\"" << value << '"';
  }

  // Writes the no-content closing tag (/>)
  void write_close(std::ostream& os) const { os << "/>\n"; }

  void render(std::ostream& os) const override {
    static_cast<const Derived&>(*this).render_impl(os);
  }
};

} // namespace svg

#endif // SVG_SVGOBJECT_H