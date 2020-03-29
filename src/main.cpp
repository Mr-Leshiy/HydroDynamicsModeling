#include <fstream>
#include <iostream>
#include <random>

#include "calculation.h"

using namespace std;

int main() {
  fstream file;
  file.precision(5);
  file.open("Result2.txt",
            std::fstream::in | std::fstream::out | std::fstream::trunc);
  initParticles(216);
  startSimulation(1000, file);

  file.close();
}
