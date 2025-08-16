#include "UF/ToricCode.h"
#include "UF/ToricCodeVisualizer.h"
#include "UF/UFDecoder.h"
#include <fstream>
#include <iostream>

#include <chrono>

using namespace uf;

int main(int argc, char* argv[]) {
  assert(argc > 1 && "Usage: uf_decoder <output_file_prefix>");
  
  constexpr int L = 5; // 11x11 toric code
  auto t0 = std::chrono::high_resolution_clock::now();

  ToricCode<L> toric;
  ToricCodeVisualizer<L> visualizer;

  std::random_device rd;
  std::mt19937_64 rng(rd());
  // toric.injectPauliNoise(0.00, 0.03, rng);
  toric.injectErasureNoise(0.1, rng);
  toric.xErr.clear(); // ignore X errors
  // toric.zErr.flip(0);
  // toric.zErr.flip(1);
  // toric.zErr.flip(2);
  // toric.zErr.flip(3);
  // toric.zErr.flip(7);
  // toric.zErr.flip(9);
  auto t2 = std::chrono::high_resolution_clock::now();

  auto zSyn = toric.computeZSyndrome(); // plaquette syndromes (X errors)
  auto xSyn = toric.computeXSyndrome(); // star syndromes (Z errors)
  auto t3 = std::chrono::high_resolution_clock::now();

  std::cout
      << "Noise injection time: "
      << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t0).count()
      << " ns\n";
  std::cout
      << "Syndrome computation time: "
      << std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count()
      << " ns\n";

  toric.zErr.display(std::cerr << "Z errors:\n");

  xSyn.display(std::cerr << "X-syndrome:\n");

  std::string fileName(argv[1]);
  std::ofstream ofs;

  visualizer.clearDocument();
  visualizer.drawLattice();
  visualizer.drawPauliErasureError(toric.xErr, toric.zErr, toric.erasureErr);
  visualizer.drawSyndromes(xSyn, zSyn);
  ofs.open(fileName + ".original.svg");
  visualizer.render(ofs);
  ofs.close();

  visualizer.clearDocument();
  visualizer.drawLattice();
  visualizer.drawError(toric.erasureErr, visualizer.erasureErrColor);
  visualizer.drawSyndromes(xSyn, zSyn);
  ofs.open(fileName + ".original.syndrome_only.svg");
  visualizer.render(ofs);
  ofs.close();
  
  std::vector<decltype(xSyn)> xSyndromes;
  xSyndromes.push_back(xSyn);

  std::vector<decltype(zSyn)> zSyndromes;
  zSyndromes.push_back(zSyn);

  auto t4 = std::chrono::high_resolution_clock::now();
  UF3D_ToricDecoder<L> decoder;
  auto zCorrection = decoder.decode(xSyndromes);
  auto xCorrection = decoder.decode(zSyndromes);
  auto t5 = std::chrono::high_resolution_clock::now();

  std::cout
      << "Decode time: "
      << std::chrono::duration_cast<std::chrono::nanoseconds>(t5 - t4).count()
      << " ns\n";

  zCorrection.display(std::cerr << "Z correction:\n");

  for (int i = 0; i < L * L; ++i) {
    if (zCorrection.test(i))
      toric.zErr.flip(i);
    if (xCorrection.test(i))
      toric.xErr.flip(i);
  }

  return 0;
}