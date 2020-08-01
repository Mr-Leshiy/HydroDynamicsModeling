#include "calculation.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "triangulation.h"
#include "utils.h"

#define NUM_THREADS 1

using namespace std;

// Static variable
static vector<HParticle> particles;
static vector<HParticle> new_particles;
static vector<HParticle>
    image_particles;  // images of particles for periodic conditions
static Gt::Point_3 system_velocity(0, 0, 0);
static double system_mass = 0;

static Triangulation triangulation;

static const double volume = 60 * 60 * 60;
static const double time_step = 1;
static double elapsed_time = 0;
static const double box_size = std::cbrt(volume);
static const double bulk_viscosity = 3.0272 / 166;
static const double shear_viscosity = 9.0898 / 166;
static const double boltzmann_constant = 1.380649 / 1.66;
static const double system_tempreture = 86.5;

class Timer {
 public:
  Timer() : beg_(clock_::now()) {}
  void reset() { beg_ = clock_::now(); }
  double elapsed() const {
    return std::chrono::duration_cast<second_>(clock_::now() - beg_).count();
  }

 private:
  typedef std::chrono::high_resolution_clock clock_;
  typedef std::chrono::duration<double, std::ratio<1> > second_;
  std::chrono::time_point<clock_> beg_;
};

HParticle::HParticle(double x, double y, double z)
    : coordinates(x, y, z),
      mass(INIT_MASS),
      density(INIT_DENSITY),
      volume(INIT_VOLUME),
      temperature(INIT_TEMPERATURE),
      velocity(Gt::Point_3(fRand(0, INIT_VELOCITY),
                           fRand(0, INIT_VELOCITY),
                           fRand(0, INIT_VELOCITY))) {}

HParticle HParticle::get_shifted_copy(const double& shift_x,
                                      const double& shift_y,
                                      const double& shift_z) {
  HParticle result = *this;
  result.coordinates = Gt::Point_3(this->coordinates.x() + shift_x,
                                   this->coordinates.y() + shift_y,
                                   this->coordinates.z() + shift_z);
  return result;
}

static double calcPressure(double density, Gt::Point_3 velocity) {
  static double a, b, c, d, e, f;
  a = -1.86789e-4;
  b = 1.9878e-1;
  c = -9.6160e1;
  d = 8.2666e4;
  e = -1.5356e6;
  f = 1.1584e-7;

  density *= 1.66e3;
  velocity *= 1000;

  /*return (f * pow(density, 5) + a * pow(density, 4) + b * pow(density, 3) +
          c * pow(density, 2) + d * density + e) / 1.66e9;*/

  return density * pow(calculate_absolute_value(velocity), 2) / 3 / 1.66e9;
}

static void pereodicConditions(HParticle& particle) {
  if (particle.coordinates.x() >= box_size)
    particle.coordinates -= Gt::Point_3(box_size, 0, 0) *
                            (int32_t)(particle.coordinates.x() / box_size);
  if (particle.coordinates.x() < 0)
    particle.coordinates += Gt::Point_3(box_size, 0, 0) *
                            (int32_t)(1 - particle.coordinates.x() / box_size);

  if (particle.coordinates.y() >= box_size)
    particle.coordinates -= Gt::Point_3(0, box_size, 0) *
                            (int32_t)(particle.coordinates.y() / box_size);
  if (particle.coordinates.y() < 0)
    particle.coordinates += Gt::Point_3(0, box_size, 0) *
                            (int32_t)(-particle.coordinates.y() / box_size + 1);

  if (particle.coordinates.z() >= box_size)
    particle.coordinates -= Gt::Point_3(0, 0, box_size) *
                            (int32_t)(particle.coordinates.z() / box_size);
  if (particle.coordinates.z() < 0)
    particle.coordinates += Gt::Point_3(0, 0, box_size) *
                            (int32_t)(1 - particle.coordinates.z() / box_size);
}

static Gt::Point_3 calcForce(const int& index, fstream& file, bool display) {
  Gt::Point_3 term1(0, 0, 0);
  Gt::Point_3 term2(0, 0, 0);

  double random_term_x = 0;
  double random_term_y = 0;
  double random_term_z = 0;

  double matrix[Tetrahedron::matrix_size][Tetrahedron::matrix_size];

  if (display) {
    file << endl;
    file << "index:( " << index
         << " ), || tempreture: " << particles[index].temperature
         << ", || pressure: "
         << calcPressure(particles[index].density, particles[index].velocity)
         << " || tets count: " << particles[index].tets.size() << " || ";
  }

  for (const auto& tet : particles[index].tets) {
    if (display) {
      file << "||||||||||||||||||||||||" << endl;
      file << "vec. B: (" << tet.vectorB.x() << ", " << tet.vectorB.y() << ", "
           << tet.vectorB.z()
           << " ) || pressure: " << calcPressure(tet.density, tet.velocity)
           << " || volume: " << tet.tetrahedron.volume()
           << " || density: " << tet.density << " || tet. velocity: ("
           << tet.velocity.x() << ", " << tet.velocity.y() << ", "
           << tet.velocity.z() << " ) || ";
    }

    Gt::Point_3 s1 = tet.vectorB * calcPressure(tet.density, tet.velocity);
    Gt::Point_3 s2 = tet.velocity * tet.density * (tet.vectorB * tet.velocity);
    term1 += (s1 + s2) * tet.tetrahedron.volume();

    double coef =
        boltzmann_constant * tet.temperature * tet.tetrahedron.volume();
    double coef1 = std::sqrt(4 * coef * bulk_viscosity);
    double coef2 = std::sqrt(3 * coef * shear_viscosity);

    double dW_avarage1 = (tet.dW[0][0] + tet.dW[1][1] + tet.dW[2][2]) / 3;
    double dW_avarage2 = (tet.dW[2][0] + tet.dW[1][1] + tet.dW[0][2]) / 3;

    for (int i = 0; i < Tetrahedron::matrix_size; ++i) {
      for (int j = 0; j < Tetrahedron::matrix_size; ++j) {
        matrix[i][j] = (tet.dW[i][j] + tet.dW[j][i]) / 2;
        if (i == j) {
          matrix[i][j] -= dW_avarage1;
        }

        matrix[i][j] *= coef1;

        if (j == Tetrahedron::matrix_size - i - 1) {
          matrix[i][j] += coef2 * dW_avarage2;
        }
      }
    }

    random_term_x += tet.vectorB.x() * matrix[0][0] +
                     tet.vectorB.y() * matrix[1][0] +
                     tet.vectorB.z() * matrix[2][0];

    random_term_y += tet.vectorB.x() * matrix[0][1] +
                     tet.vectorB.y() * matrix[1][1] +
                     tet.vectorB.z() * matrix[2][1];

    random_term_z += tet.vectorB.x() * matrix[0][2] +
                     tet.vectorB.y() * matrix[1][2] +
                     tet.vectorB.z() * matrix[2][2];

    for (const auto& neighbour : tet.particles) {
      for (const auto& n_tet : neighbour->tets) {
        if (n_tet.tetrahedron == tet.tetrahedron) {
          double scalar_product1 =
              (tet.velocity - neighbour->velocity) * n_tet.vectorB;
          double scalar_product2 = tet.vectorB * n_tet.vectorB;
          double scalar_product3 =
              (tet.velocity - neighbour->velocity) * tet.vectorB;

          if (tet.temperature == 0 && neighbour->temperature == 0) {
            term2 += (tet.vectorB * bulk_viscosity * scalar_product1 +
                      ((tet.velocity - neighbour->velocity) * scalar_product2 +
                       n_tet.vectorB * scalar_product3 -
                       tet.vectorB * 2.0 / 3.0 * scalar_product1) *
                          shear_viscosity) *
                     tet.tetrahedron.volume();
          } else {
            term2 += (tet.vectorB * bulk_viscosity * scalar_product1 +
                      ((tet.velocity - neighbour->velocity) * scalar_product2 +
                       n_tet.vectorB * scalar_product3 -
                       tet.vectorB * 2.0 / 3.0 * scalar_product1) *
                          shear_viscosity) *
                     tet.temperature * tet.tetrahedron.volume() /
                     neighbour->temperature;
          }
          break;
        }
      }
    }
  }

  auto result = (term1 + term2 /*+
                 Gt::Point_3(random_term_x, random_term_y, random_term_z)*/) /
                particles[index].volume;

  if (display) {
    file << "||||||||||||||||||||||||||||" << endl;
    file << " term1: (" << term1.x() << ", " << term1.y() << ", " << term1.z()
         << " ) || term2: ( " << term2.x() << ", " << term2.y() << ", "
         << term2.z() << " || " << endl;
  }

  return result;
}

static void calcNewVelocity(const int& index, fstream& file, bool display) {
  Gt::Point_3 force = calcForce(index, file, display);

  new_particles[index].velocity =
      particles[index].velocity + force * time_step / particles[index].density;
  system_velocity += new_particles[index].velocity;
}

static void calcNewDencity(const int& index) {
  double sum = 0;
  for (const auto& tet : particles[index].tets) {
    sum +=
        tet.tetrahedron.volume() * (tet.vectorB * tet.velocity) * tet.density;
  }
  new_particles[index].density =
      particles[index].density + time_step * sum / particles[index].volume;
}

static void display_particle(fstream& file, const HParticle& particle) {
  auto momentum = particle.velocity * (particle.volume * particle.density);

  file << fixed;
  file << "coords: (" << particle.coordinates.x() << ", "
       << particle.coordinates.y() << ", " << particle.coordinates.z()
       << ")  \t";
  file << scientific;
  file << " mass:" << particle.volume * particle.density << " \t"
       << " temp:" << particle.temperature << " \t"
       << "density:" << particle.density << " \t"
       << "volume:" << particle.volume << " \t"
       << "vel: (" << particle.velocity.x() << ", " << particle.velocity.y()
       << ", " << particle.velocity.z() << ") \t"
       << "momentum: (" << momentum.x() << ", " << momentum.y() << ", "
       << momentum.z() << ") \t"
       << "abs. mom.:"
       << std::sqrt(momentum.x() * momentum.x() + momentum.y() * momentum.y() +
                    momentum.z() * momentum.z())
       << endl
       << "---------" << endl;
}

static void display_results(fstream& file, const int& iteration) {
  file << "__________________ Series " << iteration << "__________________"
       << endl;
  double den_avg = 0, sys_mass = 0, temp_avg = 0, vel_avg = 0, mom_avg = 0,
         sys_volume = 0;
  size_t tets_count = 0;
  Gt::Point_3 sys_momentum(0, 0, 0);
  for (int i = 0; i < particles.size(); i++) {
    sys_mass += (particles[i].density * particles[i].volume);
    den_avg += particles[i].density / particles.size();
    temp_avg += particles[i].temperature / particles.size();
    sys_volume += particles[i].volume;
    sys_momentum +=
        particles[i].velocity * (particles[i].volume * particles[i].density);
    file << "-----[" << i << "]-- sys. mom.: (" << sys_momentum.x() << ", "
         << sys_momentum.y() << ", " << sys_momentum.z() << ")" << endl;

    tets_count += particles[i].tets.size();
    display_particle(file, particles[i]);
  }
  file << endl
       << "avg density: " << den_avg << "\t sys volume: " << sys_volume
       << "\t tets count: " << tets_count << "\t sys mass: " << sys_mass
       << "\t avg temp: " << temp_avg << "\t vel: " << vel_avg
       << "\t mom: " << mom_avg << "\t sys_momentum: (" << sys_momentum.x()
       << ", " << sys_momentum.y() << ", " << sys_momentum.z() << ") \t"
       << "abs. mom.:"
       << std::sqrt(sys_momentum.x() * sys_momentum.x() +
                    sys_momentum.y() * sys_momentum.y() +
                    sys_momentum.z() * sys_momentum.z())
       << endl;
}

static void calculateStep(fstream& file, bool display) {
  system_velocity = Gt::Point_3(0, 0, 0);

  Timer timer;

  for (int i = 0; i < particles.size(); ++i) {
    particles[i].temperature =
        calcPressure(particles[i].density, particles[i].velocity) * 1e6 /
        particles[i].density / 208.13;
  }

  printf("Temperature count, time : %f \n", timer.elapsed());
  timer.reset();

  triangulation.calculateParametrs(box_size, time_step);

  printf("CalculateParametrs count, time : %f \n", timer.elapsed());
  timer.reset();

  for (int i = 0; i < particles.size(); ++i) {
    calcNewVelocity(i, file, display);
    calcNewDencity(i);
  }

  printf("Velocity count, time : %f \n", timer.elapsed());
  timer.reset();

  system_velocity /= particles.size();

  for (int index = 0; index < particles.size(); ++index) {
    new_particles[index].velocity -= system_velocity;
    new_particles[index].coordinates =
        particles[index].coordinates +
        new_particles[index].velocity * time_step * 2;
    pereodicConditions(new_particles[index]);
  }

  printf("Coordinates count, time : %f \n", timer.elapsed());
  timer.reset();
}

static void swap_and_clear() {
  elapsed_time += time_step;
  for (int i = 0; i < particles.size(); ++i) {
    particles[i].neighbours_points.clear();
    particles[i].tets.clear();
    particles[i].density = new_particles[i].density;
    particles[i].coordinates = new_particles[i].coordinates;
    particles[i].velocity = new_particles[i].velocity;
    particles[i].volume = 0;
  }
}

void initParticles(int n) {
  particles.resize(n);
  new_particles.resize(n);

  std::default_random_engine generator;
  std::normal_distribution<double> distribution;

  for (int i = 0; i < n; ++i) {
    double x = fRand(0, box_size);
    double y = fRand(0, box_size);
    double z = fRand(0, box_size);

    particles[i] = HParticle(x, y, z);
    particles[i].velocity = Gt::Point_3(distribution(generator),
                                        distribution(generator),
                                        distribution(generator));

    system_velocity += particles[i].velocity;
  }

  system_velocity /= particles.size();

  for (int i = 0; i < n; ++i) {
    particles[i].velocity -= system_velocity;
  }
  triangulation = Triangulation(0, box_size, 0, box_size, 0, box_size);

  triangulation.initPoints(particles);
}

void startSimulation(const int& num_iterations, fstream& file) {
  for (int i = 0; i < num_iterations; ++i) {
    bool display = false;
    if (i > 390) {
      display = true;
    }

    printf("-- %d --\n", i);
    Timer timer;
    file << "__________________ Series " << i << "__________________" << endl;
    printf("---- values --- \n");
    calculateStep(file, display);
    printf("--------");
    if (display) {
      display_results(file, i);
    }
    swap_and_clear();
    printf("Triangulation... \n");
    triangulation.initPoints(particles);
    printf("Iteration , time : %f \n", timer.elapsed());
    timer.reset();
  }
}
