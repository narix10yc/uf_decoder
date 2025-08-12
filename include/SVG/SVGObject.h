#ifndef SVG_SVGOBJECT_H
#define SVG_SVGOBJECT_H

#include "Utils/SmallString.h"
#include <map>

namespace svg {

utils::SmallString to_ss(float v);

class SVGObjectBase {
protected:
  enum Kind { LineKind, PathKind, RectKind, CircleKind, TextKind, EndKind_ };
  Kind kind_;

public:
  explicit SVGObjectBase(Kind kind) : kind_(kind) {}

  Kind getKind() const { return kind_; }

  virtual ~SVGObjectBase() = default;
  virtual void render(std::ostream& os) const = 0;

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