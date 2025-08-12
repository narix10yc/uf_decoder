#include "UF/SurfaceCode.h"
#include <fstream>
#include <iostream>

#include <chrono>

using namespace uf;

int main(int argc, char* argv[]) {
  svg::Document doc(600, 600);
  doc.setBackgroundColor("white");
  auto t0 = std::chrono::high_resolution_clock::now();

  SurfaceCode<11> toric(true); // 11x11 grid, periodic
  auto t1 = std::chrono::high_resolution_clock::now();

  std::mt19937_64 rng;
  toric.injectPauliNoise(0.1, 0.1, rng);
  auto t2 = std::chrono::high_resolution_clock::now();

  auto zSyn = toric.computeZSyndrome(); // plaquette syndromes (X errors)
  auto xSyn = toric.computeXSyndrome(); // star syndromes (Z errors)
  auto t3 = std::chrono::high_resolution_clock::now();

  std::cout
      << "SurfaceCode<11> toric construction time: "
      << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()
      << " ns\n";
  std::cout
      << "Noise injection time: "
      << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()
      << " ns\n";
  std::cout
      << "Syndrome computation time: "
      << std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count()
      << " ns\n";

  std::cout << "#qubits=" << toric.numQubits()
            << "  #plaquettes=" << toric.numPlaquettes()
            << "  #stars=" << toric.numStars() << "\n";

  toric.xErr_.display(std::cerr << "X errors:\n");
  toric.zErr_.display(std::cerr << "Z errors:\n");

  zSyn.display(std::cerr << "Z-syndrome:\n");
  xSyn.display(std::cerr << "X-syndrome:\n");

  toric.generateSVG(doc);
  doc.addText(10, 20, "Toric Code Example").setFontSize(16);
  std::ofstream ofs(argv[1]);
  if (!ofs) {
    std::cerr << "Error opening output file!\n";
    return 1;
  }
  doc.render(ofs);
  ofs.close();
  return 0;
}