#ifndef SVG_DOCUMENT_H
#define SVG_DOCUMENT_H

#include "SVG/Circle.h"
#include "SVG/Line.h"
#include "SVG/Rect.h"
#include "SVG/Text.h"
#include "SVG/Path.h"

namespace svg {

class Document {
  int width_ = 0;
  int height_ = 0;
  utils::SmallString bgColor_; // optional background (empty = none)

  std::vector<std::unique_ptr<SVGObjectBase>> objects_;

public:
  Document(int w, int h) : width_(w), height_(h), bgColor_() {}

  Document& setBackgroundColor(utils::SmallString color) {
    bgColor_ = std::move(color);
    return *this;
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
    os << R"(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)"
       << "\n<svg xmlns=\"http://www.w3.org/2000/svg\" "
       << "width=\"" << width_ << "\" height=\"" << height_ << "\" "
       << "viewBox=\"0 0 " << width_ << ' ' << height_ << "\">\n";
    if (!bgColor_.empty()) {
      Rect(0, 0, float(width_), float(height_))
          .setFill(bgColor_)
          .setStroke("none")
          .setStrokeWidth(0)
          .render(os);
    }
    for (const auto& n : objects_)
      n->render(os);
    os << "</svg>\n";
  }
};
}; // namespace svg

#endif // SVG_DOCUMENT_H