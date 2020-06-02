#ifndef UTILS_H
#define UTILS_H

#include <time.h>

#include <cmath>
#include <iostream>

template <typename FuncType>
void workingTime(FuncType func) {
  float fTimeStart = clock() / (float)CLOCKS_PER_SEC;
  func();
  float fTimeStop = clock() / (float)CLOCKS_PER_SEC;
  printf("function time %f sec\n", fTimeStop - fTimeStart);
}

template <typename FuncType, typename FuncParamType>
void workingTime(FuncType func, FuncParamType arg) {
  float fTimeStart = clock() / (float)CLOCKS_PER_SEC;
  func(arg);
  float fTimeStop = clock() / (float)CLOCKS_PER_SEC;
  printf("function time %f sec\n", fTimeStop - fTimeStart);
}

inline double fRand(double fMin, double fMax) {
  double f = (double)rand() / RAND_MAX;
  f = fMin + f * (fMax - fMin);
  // if (f == fMax) return f - 0.0000001;
  return f;
}

inline double calculate_absolute_value(const Gt::Point_3& point) {
  return sqrt(point.x() * point.x() + point.y() * point.y() +
              point.z() * point.z());
}

inline int positive_mod(int x, int y) { return (x % y + y) % y; }

inline void operator/=(Gt::Point_3& p, const double divider) {
  p = Gt::Point_3(p.x() / divider, p.y() / divider, p.z() / divider);
}

inline void operator*=(Gt::Point_3& p, const double multiplier) {
  p = Gt::Point_3(p.x() * multiplier, p.y() * multiplier, p.z() * multiplier);
}

inline void operator+=(Gt::Point_3& p1, const Gt::Point_3& p2) {
  p1 = Gt::Point_3(p1.x() + p2.x(), p1.y() + p2.y(), p1.z() + p2.z());
}

inline void operator-=(Gt::Point_3& p1, const Gt::Point_3& p2) {
  p1 = Gt::Point_3(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
}

inline Gt::Point_3 operator+(const Gt::Point_3& p1, const Gt::Point_3& p2) {
  return Gt::Point_3(p1.x() + p2.x(), p1.y() + p2.y(), p1.z() + p2.z());
}

inline Gt::Point_3 operator-(const Gt::Point_3& p1, const Gt::Point_3& p2) {
  return Gt::Point_3(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
}

inline Gt::Point_3 operator/(const Gt::Point_3& p1, const double& divider) {
  return Gt::Point_3(p1.x() / divider, p1.y() / divider, p1.z() / divider);
}

inline Gt::Point_3 operator*(const Gt::Point_3& p1, const double& multiplier) {
  return Gt::Point_3(
      p1.x() * multiplier, p1.y() * multiplier, p1.z() * multiplier);
}

inline Gt::Point_3 operator*(const Gt::Point_3& p1, const int& multiplier) {
  return Gt::Point_3(
      p1.x() * multiplier, p1.y() * multiplier, p1.z() * multiplier);
}

inline double operator*(const Gt::Point_3& p1, const Gt::Point_3& p2) {
  return p1.x() * p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
}

inline bool has_corner(const P3DT3::Tetrahedron& tet,
                       const HParticle& particle,
                       int& out_param_corner) {
  for (int corner = 0; corner < 4; corner++) {
    if (tet[corner].x() == particle.coordinates.x() &&
        tet[corner].y() == particle.coordinates.y() &&
        tet[corner].z() == particle.coordinates.z()) {
      out_param_corner = corner;
      return true;
    }
  }
  return false;
}

inline bool isInside(const std::vector<HParticle>& particles,
                     const HParticle& p) {
  for (auto particle : particles) {
    if (particle.coordinates.x() == p.coordinates.x() &&
        particle.coordinates.y() == p.coordinates.y() &&
        particle.coordinates.z() == p.coordinates.z())
      return true;
  }
  return false;
}

inline bool isInside(const std::vector<std::shared_ptr<HParticle>>& particles,
                     const Gt::Point_3& p) {
  for (auto particle : particles) {
    if (particle->coordinates.x() == p.x() &&
        particle->coordinates.y() == p.y() &&
        particle->coordinates.z() == p.z())
      return true;
  }
  return false;
}

inline bool compare_double(const double& a, const double& b) {
  const static double eps = 0.0000001;
  return abs(a - b) < eps;
}
#endif
