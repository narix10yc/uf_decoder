#include "UF/UFDecoder3D.h"
#include <iostream>
using namespace uf;

int main(int argc, char** argv) {

  int L = 5;
  int T = 3;
  BitArray erasure(3 * L * L * T);
  BitArray syndrome(L * L * T);

  syndrome.set(0);
  syndrome.set(1);

  uf::clustering(erasure, syndrome, L, T);

  std::cerr << "Erasure is now\n" << erasure << "\n";

  return 0;
}