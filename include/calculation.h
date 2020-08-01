#ifndef CALCULATION_H
#define CALCULATION_H

#include <fstream>
#include <memory>

#include "triangulation.h"

struct HParticle {
  constexpr static double INIT_VELOCITY = -5e-1;
  // constexpr static double INIT_VELOCITY = 0;
  // constexpr static double INIT_DENSITY = 996.3 / 1.66e3;
  constexpr static double INIT_DENSITY = 1400 / 1.66e3;
  constexpr static double INIT_MASS = 0;
  constexpr static double INIT_VOLUME = 0;
  constexpr static double INIT_TEMPERATURE = 0;

  Gt::Point_3 coordinates;
  Gt::Point_3 velocity;
  double mass;
  double density;
  double volume;
  double temperature;
  std::vector<std::shared_ptr<HParticle>> neighbours_points;
  std::vector<Tetrahedron> tets;

  HParticle() {}
  HParticle(double x, double y, double z);

  // return a copy with new coordinates
  HParticle get_shifted_copy(const double& x, const double& y, const double& z);

  //~HParticle() = default;
  // HParticle& operator=(const HParticle&) = default;
};

void initParticles(int);

void startSimulation(const int&, std::fstream&);
#endif  // !
