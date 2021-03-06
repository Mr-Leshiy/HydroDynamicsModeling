#include <math.h>

#include <iostream>
#include <random>

#include "calculation.h"
#include "triangulation.h"
#include "utils.h"

double Tetrahedron::calculate_area(const Gt::Point_3& p1,
                                   const Gt::Point_3& p2,
                                   const Gt::Point_3& p3) {
  double a = calculate_absolute_value(p1 - p2);
  double b = calculate_absolute_value(p2 - p3);
  double c = calculate_absolute_value(p3 - p1);

  double p = (a + b + c) / 2.0;

  return sqrt(p * (p - a) * (p - b) * (p - c));
}

void Tetrahedron::initVectorB(const int& corner) {
  Gt::Point_3 p_0 = tetrahedron[corner];
  Gt::Point_3 p_1 = tetrahedron[positive_mod(corner + 1, 4)];
  Gt::Point_3 p_2 = tetrahedron[positive_mod(corner + 2, 4)];
  Gt::Point_3 p_3 = tetrahedron[positive_mod(corner + 3, 4)];

  vectorB = p_0 - (p_1 + p_2 + p_3) / 3;
  double s = calculate_area(p_0, p_1, p_2) + calculate_area(p_0, p_2, p_3) +
             calculate_area(p_0, p_1, p_3) + calculate_area(p_1, p_2, p_3);
  vectorB /= (2 * s);
}

Triangulation::Triangulation(const double& x1,
                             const double& x2,
                             const double& y1,
                             const double& y2,
                             const double& z1,
                             const double& z2) {
  P3DT3::Iso_cuboid domain_size(Gt::Point_3(x1, y1, z1),
                                Gt::Point_3(x2, y2, z2));
  this->delaunay_triangulation.set_domain(domain_size);
}

void Triangulation::initPoints(std::vector<HParticle>& init_points) {
  this->delaunay_triangulation.clear();
  for (size_t i = 0; i < init_points.size(); ++i) {
    delaunay_triangulation.insert(Gt::Point_3(init_points[i].coordinates.x(),
                                              init_points[i].coordinates.y(),
                                              init_points[i].coordinates.z(),
                                              &init_points[i]));
  }
}

void Triangulation::calculateParametrs(const double& box_size,
                                       const double& time_step) {

  P3DT3::Point_3 devider(delaunay_triangulation.domain().xmax(),
                         delaunay_triangulation.domain().ymax(),
                         delaunay_triangulation.domain().zmax());

  auto& domain = delaunay_triangulation.domain();
  P3DT3::Tetrahedron tetrahedron;
  for (P3DT3::Periodic_tetrahedron_iterator ptit =
           delaunay_triangulation.periodic_tetrahedra_begin(
               P3DT3::UNIQUE_COVER_DOMAIN);
       ptit != delaunay_triangulation.periodic_tetrahedra_end(
                   P3DT3::UNIQUE_COVER_DOMAIN);
       ++ptit) {
    tetrahedron = delaunay_triangulation.construct_tetrahedron(*ptit);

    double sum_density = 0;
    double sum_tempreture = 0;
    Gt::Point_3 sum_velocity(0, 0, 0);
    for (int corner = 0; corner < 4; ++corner) {
      sum_density += tetrahedron[corner].store_el->density;
      sum_tempreture += tetrahedron[corner].store_el->temperature;
      sum_velocity += tetrahedron[corner].store_el->velocity;
    }

    static std::default_random_engine generator;
    static std::normal_distribution<double> distribution(0.0, 1.0);
    double dW[Tetrahedron::matrix_size][Tetrahedron::matrix_size];

    for (int i = 0; i < Tetrahedron::matrix_size; ++i) {
      for (int j = 0; j < Tetrahedron::matrix_size; ++j) {
        dW[i][j] = std::sqrt(time_step) * distribution(generator);
      }
    }

    for (int corner1 = 0; corner1 < 4; ++corner1) {
      if (tetrahedron[corner1].x() < domain.xmax() &&
          tetrahedron[corner1].x() >= domain.xmin() &&
          tetrahedron[corner1].y() < domain.ymax() &&
          tetrahedron[corner1].y() >= domain.ymin() &&
          tetrahedron[corner1].z() < domain.zmax() &&
          tetrahedron[corner1].z() >= domain.zmin()) {
        Tetrahedron new_tet(tetrahedron);
        new_tet.initVectorB(corner1);
        new_tet.density = sum_density / 4;
        new_tet.temperature = sum_tempreture / 4;
        new_tet.velocity = sum_velocity / 4;

        for (int i = 0; i < Tetrahedron::matrix_size; ++i) {
          for (int j = 0; j < Tetrahedron::matrix_size; ++j) {
            new_tet.dW[i][j] = dW[i][j];
          }
        }

        for (int corner2 = 0; corner2 < 4; ++corner2) {
          if (corner1 != corner2 &&
              !isInside(tetrahedron[corner1].store_el->neighbours_points,
                        tetrahedron[corner2])) {
            std::shared_ptr<HParticle> temp = std::make_shared<HParticle>();
            temp->coordinates = tetrahedron[corner2];
            temp->density = tetrahedron[corner2].store_el->density;
            temp->mass = tetrahedron[corner2].store_el->mass;
            temp->temperature = tetrahedron[corner2].store_el->temperature;
            temp->density = tetrahedron[corner2].store_el->density;
            temp->velocity = tetrahedron[corner2].store_el->velocity;

            Tetrahedron n_tet(tetrahedron);
            n_tet.initVectorB(corner2);
            n_tet.temperature = new_tet.temperature;
            n_tet.velocity = new_tet.velocity;
            temp->tets = {n_tet};

            for (int i = 0; i < Tetrahedron::matrix_size; ++i) {
              for (int j = 0; i < Tetrahedron::matrix_size; ++i) {
                n_tet.dW[i][j] = dW[i][j];
              }
            }

            tetrahedron[corner1].store_el->neighbours_points.push_back(temp);
            new_tet.particles.push_back(temp);
          }
        }
        tetrahedron[corner1].store_el->tets.push_back(new_tet);
        tetrahedron[corner1].store_el->volume += tetrahedron.volume() / 4.0;
      }
    }
  }
}
