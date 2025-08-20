#include "UF/ToricCode.h"
#include "UF/ToricCodeVisualizer.h"
#include "UF/UFDecoder.h"
#include "UF/UFDecoder3D.h"
#include <fstream>
#include <iostream>

#include <chrono>

using namespace uf;

int main(int argc, char* argv[]) {
  assert(argc > 1 && "Usage: uf_decoder <output_file_prefix>");

  constexpr int L = 5;
  auto t0 = std::chrono::high_resolution_clock::now();
  auto t1 = std::chrono::high_resolution_clock::now();

  ToricCode<L> toric;
  ToricCodeVisualizer<L> visualizer;

  std::random_device rd;
  // std::mt19937_64 rng(rd());
  std::mt19937_64 rng(41); // fixed seed for reproducibility

  t0 = std::chrono::high_resolution_clock::now();
  // toric.injectPauliNoise(0.00, 0.3, rng);
  toric.injectErasureNoise(0.4f, rng);
  toric.xErr.clear(); // ignore X errors
                      //   toric.zErr.set(1);
                      //   toric.zErr.set(2);
                      //   toric.erasureErr.flip(1);
                      //   toric.erasureErr.flip(2);
  // toric.erasureErr.flip(3);
  // toric.erasureErr.flip(7);
  // toric.erasureErr.flip(8);
  // toric.erasureErr.flip(28);
  // toric.erasureErr.flip(29);
  t1 = std::chrono::high_resolution_clock::now();
  std::cout
      << "Noise injection time: "
      << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()
      << " ns\n";

  t0 = std::chrono::high_resolution_clock::now();
  auto xSyn = toric.computeXSyndrome(); // star syndromes (Z errors)
  t1 = std::chrono::high_resolution_clock::now();

  std::cout
      << "Syndrome computation time: "
      << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()
      << " ns\n";

  toric.zErr.display(std::cerr << "Z errors:\n");

  xSyn.display(std::cerr << "X-syndrome:\n");

  std::string fileName(argv[1]);
  std::ofstream ofs;

  visualizer.clearDocument();
  visualizer.drawLattice();
  visualizer.drawPauliErasureError(toric.xErr, toric.zErr, toric.erasureErr);
  visualizer.drawSyndromeX(xSyn);
  ofs.open(fileName + ".original.svg");
  visualizer.render(ofs);
  ofs.close();

  visualizer.clearDocument();
  visualizer.drawLattice();
  visualizer.drawError(toric.erasureErr, visualizer.erasureErrColor);
  visualizer.drawSyndromeX(xSyn);
  ofs.open(fileName + ".original.syndrome_only.svg");
  visualizer.render(ofs);
  ofs.close();

  t0 = std::chrono::high_resolution_clock::now();
  auto adjacencies = build_star_adjacencies<L>();
  t1 = std::chrono::high_resolution_clock::now();
  std::cout
      << "Adjacency matrix build time: "
      << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()
      << " ns\n";

  t0 = std::chrono::high_resolution_clock::now();
  auto correction =
      peeling_decode_Z<L>(toric.erasureErr, xSyn, adjacencies.get());
  t1 = std::chrono::high_resolution_clock::now();
  std::cout
      << "Peeling decode time: "
      << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()
      << " ns\n";

  std::cerr << "Peeling decode correction:\n";
  correction.display(std::cerr);

  //   toric.zErr ^= correction;
  //   xSyn = toric.computeXSyndrome(); // recompute X syndrome after correction

  //   visualizer.clearDocument();
  //   visualizer.drawLattice();
  //   visualizer.drawPauliError(toric.xErr, toric.zErr);
  //   visualizer.drawSyndromeX(xSyn);
  //   ofs.open(fileName + ".decoded.syndrome_only.svg");
  //   visualizer.render(ofs);
  //   ofs.close();

  int T = 2;
  int nEdges = 3 * L * L * T - L * L;
  int nStars = L * L * T;

  BitArray correction3d(nEdges);
  BitArray erasure3d(nEdges);
  BitArray syndromes3d(nStars);

  for (int i = 0; i < 2 * L * L; ++i) {
    if (toric.erasureErr.test(i))
      erasure3d.set(i);
    if (correction.test(i))
      correction3d.set(i);
  }

  for (int i = 0; i < nStars; ++i) {
    if (xSyn.test(i))
      syndromes3d.set(i);
  }

  std::cerr << "Decode3D erasure:\n" << erasure3d;
  std::cerr << "Decode3D syndromes:\n" << syndromes3d;
  t0 = std::chrono::high_resolution_clock::now();
  uf::peelingDecode(correction3d, erasure3d, syndromes3d, L, T);
  t1 = std::chrono::high_resolution_clock::now();
  std::cout
      << "Peeling decode 3D time: "
      << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()
      << " ns\n";

  std::cerr << "Peeling decode correction 3D:\n" << correction3d << '\n';

  return 0;
}