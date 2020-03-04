#include <fstream>
#include <iostream>

#include "calculation.h"

using namespace std;

int main() {
  fstream file;
  file.precision(5);
  file.open("Result.txt",
            std::fstream::in | std::fstream::out | std::fstream::trunc);
  initParticles(216);
  startSimulation(2000000, file);

  file.close();
}
