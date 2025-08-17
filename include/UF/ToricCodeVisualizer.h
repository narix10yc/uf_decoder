#ifndef UF_TORICCODE_VISUALIZER_H
#define UF_TORICCODE_VISUALIZER_H

#include "UF/ToricCode.h"

#include "SVG/Document.h"

namespace uf {

template <int L> class ToricCodeVisualizer {
  using Error = typename ToricCode<L>::Error;
  using Syndrome = typename ToricCode<L>::Syndrome;

private:
  svg::Document doc_;
  float canvasPadding_ = 5.0f; // Padding around the canvas
  utils::SmallString canvasBackgroundColor_ = "white";

public:
  float latticePadding = 15.0f;  // Padding around the lattice
  float cellSize = 40.0f;        // Size of each cell in the lattice
  float edgeStrokeWidth = 1.5f;  // Stroke width for edges
  float errorStrokeWidth = 3.0f; // Stroke width for error edges

  utils::SmallString latticeEdgeColor = "#888";
  utils::SmallString xErrColor = "red";
  utils::SmallString yErrColor = "green";
  utils::SmallString zErrColor = "blue";
  utils::SmallString erasureErrColor = "orange";

  utils::SmallString xSyndromeColor = "red";
  utils::SmallString zSyndromeColor = "blue";

  // Radius for X syndrome circles, in unit of cell size
  float xSyndromeRadius = 0.1f;
  // Size for Z syndrome squares, in unit of cell size
  float zSyndromeSize = 0.8f;

public:
  ToricCodeVisualizer() : doc_() {
    doc_.setBackgroundColor(canvasBackgroundColor_)
        .autosize()
        .setMargin(canvasPadding_);
  }

  svg::Document& getDocument() { return doc_; }
  const svg::Document& getDocument() const { return doc_; }

  void clearDocument() { doc_.clear(); }

  void setCanvasPadding(float padding) {
    canvasPadding_ = padding;
    doc_.setMargin(padding);
  }

  void setCanvasBackgroundColor(utils::SmallString color) {
    canvasBackgroundColor_ = std::move(color);
    doc_.setBackgroundColor(canvasBackgroundColor_);
  }

  // Get the position of a star (vertex in the lattice) in the canvas
  std::pair<float, float> getStarPosition(int r, int c) const {
    assert(0 <= r && r < L);
    assert(0 <= c && c < L);
    return {latticePadding + c * cellSize, latticePadding + r * cellSize};
  }

  // Get the position of the center of a plaquette (face of the lattice) in the
  // canvas
  std::pair<float, float> getPlaquettePosition(int r, int c) const {
    assert(0 <= r && r < L);
    assert(0 <= c && c < L);
    return {latticePadding + (c + 0.5f) * cellSize,
            latticePadding + (r + 0.5f) * cellSize};
  }

  // Get the start position of an edge (qubit) in the canvas. The start position
  // of (H,r,c) and (V,r,c) edges is the same.
  std::pair<float, float> getEdgeStartPosition(int r, int c) const {
    assert(0 <= r && r < L);
    assert(0 <= c && c < L);
    return {latticePadding + c * cellSize, latticePadding + r * cellSize};
  }

  void drawLattice() {
    auto& path = doc_.addPath()
                     .setStrokeWidth(edgeStrokeWidth)
                     .setStroke(latticeEdgeColor);
    // horizontal lines
    for (int row = 0; row < L; ++row) {
      path.moveTo(latticePadding, latticePadding + row * cellSize)
          .hTo(latticePadding + L * cellSize);
    }
    // vertical lines
    for (int col = 0; col < L; ++col) {
      path.moveTo(latticePadding + col * cellSize, latticePadding)
          .vTo(latticePadding + L * cellSize);
    }
  }

  void drawError(const Error& error, utils::SmallString color) {
    auto& path =
        doc_.addPath().setStrokeWidth(errorStrokeWidth).setStroke(color);

    for (int r = 0; r < L; ++r) {
      for (int c = 0; c < L; ++c) {
        // horizontal edges
        int idx_H = ToricCode<L>::qubitIdx_HG(r, c);
        if (error.test(idx_H)) {
          path.moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .hTo(latticePadding + (c + 1) * cellSize);
        }
        // vertical edges
        int idx_V = ToricCode<L>::qubitIdx_VG(r, c);
        if (error.test(idx_V)) {
          path.moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .vTo(latticePadding + (r + 1) * cellSize);
        }
      }
    }
  }

  // A more succinct way to draw Pauli errors (X, Y, and Z)
  void drawPauliError(const Error& xErr, const Error& zErr) {
    auto& pathX =
        doc_.addPath().setStrokeWidth(errorStrokeWidth).setStroke(xErrColor);
    auto& pathY =
        doc_.addPath().setStrokeWidth(errorStrokeWidth).setStroke(yErrColor);
    auto& pathZ =
        doc_.addPath().setStrokeWidth(errorStrokeWidth).setStroke(zErrColor);

    for (int r = 0; r < L; ++r) {
      for (int c = 0; c < L; ++c) {
        int x, z;

        // horizontal edges
        int idx_H = ToricCode<L>::qubitIdx_HG(r, c);
        x = xErr.test(idx_H);
        z = zErr.test(idx_H);
        if (x && z) {
          // Y error
          pathY
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .hTo(latticePadding + (c + 1) * cellSize);
        } else if (x) {
          // X error
          pathX
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .hTo(latticePadding + (c + 1) * cellSize);
        } else if (z) {
          // Z error
          pathZ
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .hTo(latticePadding + (c + 1) * cellSize);
        }

        // vertical edges
        int idx_V = ToricCode<L>::qubitIdx_VG(r, c);
        x = xErr.test(idx_V);
        z = zErr.test(idx_V);
        if (x && z) {
          // Y error
          pathY
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .vTo(latticePadding + (r + 1) * cellSize);
        } else if (x) {
          // X error
          pathX
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .vTo(latticePadding + (r + 1) * cellSize);
        } else if (z) {
          // Z error
          pathZ
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .vTo(latticePadding + (r + 1) * cellSize);
        }
      } // for c
    } // for r
  }

  // A more succinct way to draw Pauli and erasure errors at the same time
  void drawPauliErasureError(const Error& xErr,
                             const Error& zErr,
                             const Error& erasureErr) {
    auto& pathX =
        doc_.addPath().setStrokeWidth(errorStrokeWidth).setStroke(xErrColor);
    auto& pathY =
        doc_.addPath().setStrokeWidth(errorStrokeWidth).setStroke(yErrColor);
    auto& pathZ =
        doc_.addPath().setStrokeWidth(errorStrokeWidth).setStroke(zErrColor);
    auto& pathE = doc_.addPath()
                      .setStrokeWidth(errorStrokeWidth)
                      .setStroke(erasureErrColor);

    for (int r = 0; r < L; ++r) {
      for (int c = 0; c < L; ++c) {
        int x, z;

        // horizontal edges
        int idx_H = ToricCode<L>::qubitIdx_HG(r, c);
        x = xErr.test(idx_H);
        z = zErr.test(idx_H);
        if (x && z) {
          // Y error
          pathY
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .hTo(latticePadding + (c + 1) * cellSize);
        } else if (x) {
          // X error
          pathX
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .hTo(latticePadding + (c + 1) * cellSize);
        } else if (z) {
          // Z error
          pathZ
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .hTo(latticePadding + (c + 1) * cellSize);
        } else if (erasureErr.test(idx_H)) {
          // Erasure error
          pathE
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .hTo(latticePadding + (c + 1) * cellSize);
        }

        // vertical edges
        int idx_V = ToricCode<L>::qubitIdx_VG(r, c);
        x = xErr.test(idx_V);
        z = zErr.test(idx_V);
        if (x && z) {
          // Y error
          pathY
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .vTo(latticePadding + (r + 1) * cellSize);
        } else if (x) {
          // X error
          pathX
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .vTo(latticePadding + (r + 1) * cellSize);
        } else if (z) {
          // Z error
          pathZ
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .vTo(latticePadding + (r + 1) * cellSize);
        } else if (erasureErr.test(idx_V)) {
          // Erasure error
          pathE
              .moveTo(latticePadding + c * cellSize,
                      latticePadding + r * cellSize)
              .vTo(latticePadding + (r + 1) * cellSize);
        }
      } // for c
    } // for r
  }

  void drawSyndromeX(const Syndrome& syndrome) {
    for (int r = 0; r < L; ++r) {
      for (int c = 0; c < L; ++c) {
        int idx = ToricCode<L>::qubitIdx_HG(r, c);
        if (syndrome.test(idx)) {
          doc_.addCircle(latticePadding + c * cellSize,
                         latticePadding + r * cellSize,
                         cellSize * xSyndromeRadius)
              .setFill(xSyndromeColor);
        }
      }
    }
  }

  void drawSyndromeZ(const Syndrome& syndrome) {
    float shift = 0.5f * (1.0f - zSyndromeSize);
    for (int r = 0; r < L; ++r) {
      for (int c = 0; c < L; ++c) {
        int idx = ToricCode<L>::qubitIdx_HG(r, c);
        if (syndrome.test(idx)) {
          doc_.addRect(latticePadding + (c + shift) * cellSize,
                       latticePadding + (r + shift) * cellSize,
                       cellSize * zSyndromeSize,
                       cellSize * zSyndromeSize)
              .setFill(zSyndromeColor)
              .setOpacity(0.5f);
        }
      }
    }
  }

  void drawSyndromes(const Syndrome& xSyndrome, const Syndrome& zSyndrome) {
    drawSyndromeX(xSyndrome);
    drawSyndromeZ(zSyndrome);
  }

  void render(std::ostream& os) const { doc_.render(os); }

}; // class ToricCodeVisualizer

} // namespace uf

#endif // UF_TORICCODE_VISUALIZER_H