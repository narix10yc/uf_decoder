#include "SVG/Document.h"
#include <fstream>

using namespace svg;

int main(int argc, char* argv[]) {
  Document doc(800, 600);
  doc.background("#f0f0f0");
  doc.addLine(100, 100, 200, 200);
  doc.addRect(300, 300, 150, 100, 30, 30);

  // write SVG to file
  std::ofstream ofs(argv[1]);
  if (!ofs) {
    std::cerr << "Error opening output file!\n";
    return 1;
  }
  doc.render(ofs);
  ofs.close();
  return 0;
}