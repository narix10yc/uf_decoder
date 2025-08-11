#ifndef SVG_DOCUMENT_H
#define SVG_DOCUMENT_H

#include "SVG/Style.h"

namespace svg {

// ----- Helpers -----
inline void attr(std::ostream& os, const char* key, float v) {
  os << ' ' << key << "=\"" << v << '"';
}

inline void
attr(std::ostream& os, const char* key, const utils::SmallString& v) {
  os << ' ' << key << "=\"" << v << '"';
}

class SVGElementBase {
public:
  virtual ~SVGElementBase() = default;
  virtual void render(std::ostream& os) const = 0;
};

template <typename Derived> class SVGElement : public SVGElementBase {
protected:
  utils::SmallString tag_;
  svg::Style style_;

public:
  SVGElement(utils::SmallString tag, const svg::Style& style = {})
      : tag_(tag), style_(style) {}

  // Writes: <tag ... style="..."/>   (Derived adds its own attributes)
  void write_open(std::ostream& os) const { os << '<' << tag_; }
  void write_close(std::ostream& os) const { os << ' ' << style_ << " />\n"; }

  Derived& setStyle(const svg::Style& newStyle) {
    style_ = newStyle;
    return static_cast<Derived&>(*this);
  }

  Derived& setTag(utils::SmallString newTag) {
    tag_ = std::move(newTag);
    return static_cast<Derived&>(*this);
  }

  void render(std::ostream& os) const override {
    static_cast<const Derived&>(*this).render_impl(os);
  }
};

class Line : public SVGElement<Line> {
  float x1_, y1_, x2_, y2_;

public:
  Line(float x1, float y1, float x2, float y2)
      : SVGElement<Line>("line"), x1_(x1), y1_(y1), x2_(x2), y2_(y2) {}

  void render_impl(std::ostream& os) const {
    this->write_open(os);
    attr(os, "x1", x1_);
    attr(os, "y1", y1_);
    attr(os, "x2", x2_);
    attr(os, "y2", y2_);
    this->write_close(os);
  }
};

class Rect : public SVGElement<Rect> {
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
      : SVGElement<Rect>("rect"), x_(x), y_(y), width_(width), height_(height),
        rx_(rx), ry_(ry) {}

  Rect& setCornerRadius(float rx, float ry = -1.f) {
    rx_ = rx;
    ry_ = (ry < 0.f ? rx : ry);
    return *this;
  }

  void render_impl(std::ostream& os) const {
    this->write_open(os);
    attr(os, "x", x_);
    attr(os, "y", y_);
    attr(os, "width", width_);
    attr(os, "height", height_);
    if (rx_ > 0.0f || ry_ > 0.0f) {
      attr(os, "rx", rx_);
      attr(os, "ry", ry_);
    }
    this->write_close(os);
  }
};

class Document {
  int width_ = 0;
  int height_ = 0;
  utils::SmallString bg_; // optional background (empty = none)

  std::vector<std::unique_ptr<SVGElementBase>> elems_;

public:
  Document(int w, int h) : width_(w), height_(h) {}

  Document& background(utils::SmallString color) {
    bg_ = std::move(color);
    return *this;
  }

  Line& addLine(float x1, float y1, float x2, float y2) {
    elems_.emplace_back(std::make_unique<Line>(x1, y1, x2, y2));
    return static_cast<Line&>(*elems_.back());
  }

  Rect& addRect(float x,
                float y,
                float width,
                float height,
                float rx = 0.f,
                float ry = 0.f) {
    elems_.emplace_back(std::make_unique<Rect>(x, y, width, height, rx, ry));
    return static_cast<Rect&>(*elems_.back());
  }

  void render(std::ostream& os) const {
    os << R"(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)"
       << "\n<svg xmlns=\"http://www.w3.org/2000/svg\" "
       << "width=\"" << width_ << "\" height=\"" << height_ << "\" "
       << "viewBox=\"0 0 " << width_ << ' ' << height_ << "\">\n";
    if (!bg_.empty()) {
      Style s;
      s.fill = bg_;
      s.stroke = "none";
      s.stroke_width = 0;
      Rect(0, 0, float(width_), float(height_)).setStyle(s).render(os);
    }
    for (const auto& n : elems_)
      n->render(os);
    os << "</svg>\n";
  }
};
}; // namespace svg

#endif // SVG_DOCUMENT_H