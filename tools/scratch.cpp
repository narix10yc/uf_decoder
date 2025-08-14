#include "UF/ToricCode.h"
#include "UF/UFDecoder.h"
#include <fstream>
#include <iostream>

#include <chrono>

using namespace uf;

int main(int argc, char* argv[]) {
  constexpr int GRID_SIZE = 17; // 11x11 toric code
  svg::Document doc(750, 750);
  doc.setBackgroundColor("white");
  auto t0 = std::chrono::high_resolution_clock::now();

  ToricCode<GRID_SIZE> toric;

  std::mt19937_64 rng;
  toric.injectPauliNoise(0.00, 0.03, rng);
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

  doc.clear();
  toric.generateSVG(doc);
  ofs.open(fileName + ".original.svg");
  doc.render(ofs);
  ofs.close();
  
  doc.clear();
  toric.generateSVG(doc, false, true);
  ofs.open(fileName + ".original.syndrome_only.svg");
  doc.render(ofs);
  ofs.close();

  std::vector<decltype(xSyn)> xSyndromes;
  xSyndromes.push_back(xSyn);

  std::vector<decltype(zSyn)> zSyndromes;
  zSyndromes.push_back(zSyn);

  auto t4 = std::chrono::high_resolution_clock::now();
  UF3D_ToricDecoder<GRID_SIZE> decoder;
  auto zCorrection = decoder.decode(xSyndromes);
  auto xCorrection = decoder.decode(zSyndromes);
  auto t5 = std::chrono::high_resolution_clock::now();

  std::cout
      << "Decode time: "
      << std::chrono::duration_cast<std::chrono::nanoseconds>(t5 - t4).count()
      << " ns\n";

  zCorrection.display(std::cerr << "Z correction:\n");

  for (int i = 0; i < GRID_SIZE * GRID_SIZE; ++i) {
    if (zCorrection.test(i))
      toric.zErr.flip(i);
    if (xCorrection.test(i))
      toric.xErr.flip(i);
  }

  doc.clear();
  toric.generateSVG(doc);
  ofs.open(fileName + ".decoded.svg");
  doc.render(ofs);
  ofs.close();

  doc.clear();
  toric.generateSVG(doc, false, true);
  ofs.open(fileName + ".decoded.syndrome_only.svg");
  doc.render(ofs);
  ofs.close();
  return 0;
}