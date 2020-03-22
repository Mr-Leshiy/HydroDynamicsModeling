#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <memory>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K> Gt;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<Gt> P3DT3;

typedef Gt::Tetrahedron_3 Tetrahedra;

struct HParticle;

struct Tetrahedron {
  P3DT3::Tetrahedron tetrahedron;
  double density;
  double temperature;
  Gt::Point_3 velocity;
  Gt::Point_3 vectorB;

  static const int matrix_size = 3;
  double dW[matrix_size][matrix_size];

  std::vector<std::shared_ptr<HParticle>> particles;

  Tetrahedron(const P3DT3::Tetrahedron& tetrahedron)
      : tetrahedron(tetrahedron),
        density(0),
        temperature(0),
        velocity(0, 0, 0),
        vectorB(0, 0, 0) {}

  void initVectorB(const int& corner);

 private:
  double calculate_area(const Gt::Point_3&,
                        const Gt::Point_3&,
                        const Gt::Point_3&);
};

class Triangulation {
 private:
  P3DT3 delaunay_triangulation;

 public:
  Triangulation() = default;
  Triangulation(const double& x1,
                const double& x2,
                const double& y1,
                const double& y2,
                const double& z1,
                const double& z2);

  void initPoints(std::vector<HParticle>&);

  void calculateParametrs(const double& box_size, const double& time_step);
};

#endif  // !TRIANGULATION_H
