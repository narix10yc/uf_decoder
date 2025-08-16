#ifndef SVG_DOCUMENT_H
#define SVG_DOCUMENT_H

#include "SVG/Circle.h"
#include "SVG/Line.h"
#include "SVG/Path.h"
#include "SVG/Rect.h"
#include "SVG/Text.h"

#include <iostream>
#include <list>

namespace svg {

class Document {
  int width_ = 0;
  int height_ = 0;
  bool autosize_ = false;
  float margin_ = 5.0f;        // padding around the content box
  utils::SmallString bgColor_; // optional background (empty = none)

  std::vector<std::unique_ptr<SVGObjectBase>> objects_;

  Bounds compute_content_bounds_() const {
    Bounds b{0.0f, 0.0f, 0.0f, 0.0f};
    for (const auto& e : objects_) {
      auto objBounds = e->computeBounds();
      // Some objects may not have valid bounds (e.g., empty paths)
      if (objBounds.valid())
        b.include(objBounds);
    }
    b.pad(margin_);
    // We zero-initialized the bounds so it should always be valid
    assert(b.valid());
    return b;
  }

  void render_header_(std::ostream& os) const {
    Bounds box{0.0f, 0.0f, float(width_), float(height_)};
    if (width_ <= 0 && height_ <= 0 && !autosize_) {
      std::cerr
          << "Warning: SVG Document has no size set. Will be auto-sized.\n";
      box = compute_content_bounds_();
    }
    if (autosize_)
      box = compute_content_bounds_();

    os << R"(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)"
       << "\n<svg xmlns=\"http://www.w3.org/2000/svg\" "
       << "width=\"" << static_cast<int>(std::ceil(box.width())) << "\" "
       << "height=\"" << static_cast<int>(std::ceil(box.height())) << "\" "
       << "viewBox=\"" << box.xmin << ' ' << box.ymin << ' ' << box.width()
       << ' ' << box.height() << "\">\n";

    // Render background
    if (!bgColor_.empty()) {
      Rect(box.xmin, box.ymin, box.width(), box.height())
          .setFill(bgColor_)
          .setStroke("none")
          .setStrokeWidth(0)
          .render(os);
    }
  }

  void render_background_(std::ostream& os) const {
    if (!bgColor_.empty()) {
      Rect(0, 0, float(width_), float(height_))
          .setFill(bgColor_)
          .setStroke("none")
          .setStrokeWidth(0)
          .render(os);
    }
  }

public:
  Document() = default;

  Document(int w, int h) : width_(w), height_(h), bgColor_() {}

  Document& autosize(bool enable = true) {
    autosize_ = enable;
    return *this;
  }

  Document& setMargin(float m) {
    margin_ = m;
    return *this;
  }

  Document& setBackgroundColor(utils::SmallString color) {
    bgColor_ = std::move(color);
    return *this;
  }

  // Clear all elements in the document. This will not remove the background.
  void clear() { objects_.clear(); }

  void add(std::unique_ptr<SVGObjectBase> obj) {
    objects_.emplace_back(std::move(obj));
  }

  Line& addLine(float x1, float y1, float x2, float y2) {
    objects_.emplace_back(std::make_unique<Line>(x1, y1, x2, y2));
    return static_cast<Line&>(*objects_.back());
  }

  Rect& addRect(float x,
                float y,
                float width,
                float height,
                float rx = 0.f,
                float ry = 0.f) {
    objects_.emplace_back(std::make_unique<Rect>(x, y, width, height, rx, ry));
    return static_cast<Rect&>(*objects_.back());
  }

  /// @param cx, cy Center coordinates.
  /// @param r Radius.
  Circle& addCircle(float cx, float cy, float r) {
    objects_.emplace_back(std::make_unique<Circle>(cx, cy, r));
    return static_cast<Circle&>(*objects_.back());
  }

  Text& addText(float x, float y, utils::SmallString content) {
    objects_.emplace_back(std::make_unique<Text>(x, y, content));
    return static_cast<Text&>(*objects_.back());
  }

  Path& addPath() {
    objects_.emplace_back(std::make_unique<Path>());
    return static_cast<Path&>(*objects_.back());
  }

  void render(std::ostream& os) const {
    render_header_(os);

    for (const auto& n : objects_)
      n->render(os);
    os << "</svg>\n";
  }
};
}; // namespace svg

#endif // SVG_DOCUMENT_H