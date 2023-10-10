#include "header.h"

int main (int argc, char **argv) {
  Eigen::initParallel();
  Config config;

  // Finite Differences Laplace Solver Matrices, we use RowMajor for Eigen parallel algorithms
  Eigen::SparseMatrix<double, Eigen::RowMajor> I(config.N, config.N);
  I.setIdentity();
  Eigen::SparseMatrix<double, Eigen::RowMajor> Laplacian(config.N, config.N); // Finite differences Laplace operator
  set_up_finite_laplacian(config, Laplacian);
  Eigen::SparseMatrix<double, Eigen::RowMajor> Partial_dx(config.N, config.N); // Finite differences Derivative operator
  set_up_finite_partial_dx(config, Partial_dx);
  Eigen::SparseMatrix<double, Eigen::RowMajor> Partial_dy(config.N, config.N); // Finite differences Derivative operator
  set_up_finite_partial_dy(config, Partial_dy);
  Eigen::SparseMatrix<double, Eigen::RowMajor> Domain(config.N, config.N); // = Diag[dmain(x,y)] --> domain = 1 in boundary, 0 else where

  // Vectors for electric field calculation
  Eigen::VectorXd phi(config.N);
  Eigen::VectorXd electric_field_x(config.N);
  Eigen::VectorXd electric_field_y(config.N);
  Eigen::VectorXd density(config.N); density.setZero();
  std::vector<Eigen::VectorXd> electric_field = {electric_field_x, electric_field_y};

  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;
  auto t1 = high_resolution_clock::now();

  std::vector<double> particles;
  Crandom random(config.seed);
  int iteration = 0;
  double t = 0.;
  double dt = config.dt_init;
  double dt_old = dt;
  bool last = false;
  int vis_iteration = 0;
  int vis_steps = config.vis_steps_max;
  int vis_print = 0;
  int relax_iter = 0;
  double relax_res = 0.;
  int n_particles = particles.size();

  if (config.verbose) {
    std::cout << std::left << std::setw(12)
              << "-------------------------------------------------------------------------------\n"
              << std::left << std::setw(9)
              << "Progress" << std::setw(8)
              << "Step" << std::setw(9)
              << "Time" << std::setw(10)
              << "Dt" << std::setw(9)
              << "Printed" << std::setw(12)
              << "Iterations" << std::setw(12)
              << "Converged" << std::setw(9)
              << "Particles"
              << std::left << std::setw(16)
              << "\n-------------------------------------------------------------------------------\n";
  }

  n_particles = initialization(config, Domain, phi, electric_field, particles, density);
  compute_electric_potential(phi, Laplacian, Domain, I);
  compute_electric_field(electric_field, phi, Partial_dx, Partial_dy); // Still need to calculate integral and vary Voltage
  {
      //Update corresponding time_step
      double EMax = config.V/(config.Rint*log(config.Rext/config.Rint));
      double a = config.diffusivity*std::pow(EMax/config.V_ref, 2);
      double b = EMax*config.l/(2*config.V_ref);
      dt = std::min(config.dt_init, (4+b+std::sqrt(16+8*b))/a);
      dt_old = dt;
  }

  
  print(config, t, Domain, phi, electric_field, density);

  if (config.verbose){
      std::cout << std::left << std::setw(9)
                << "0%" << std::setw(8)
                << iteration << std::setw(9)
                << t  << std::setw(10)
                << dt_old << std::setw(9)
                << vis_print << std::setw(12)
                << relax_iter << std::setw(12)
                << relax_res << std::setw(9)
                << n_particles
                << "\r";
      std::cout.flush();
  }

  for (iteration = 1, vis_iteration = 1; !last; iteration++, vis_iteration++){

      //Update iteration parameters
      last = (t >= config.t_final - 1e-8*config.dt_init);
      dt = std::min(dt, config.t_final - t);

      //Perform a time step
      n_particles = system_evolve(config, dt, random, Domain, phi, electric_field, particles, density);
      t += dt;

      //Update visualization steps
      vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

      //Calculate new electric field
      dt_old = dt;
      compute_electric_potential(phi, Laplacian, Domain, I);
      compute_electric_field(electric_field, phi, Partial_dx, Partial_dy);

      if (last || vis_steps <= vis_iteration){
          //Update parameters
          vis_iteration = 0;
          vis_print += 1;

          print(config, t, Domain, phi, electric_field, density, "results/data/data_"+std::to_string(vis_print));
      }

      if (config.verbose){
          std::cout << std::left << std::setw(9)
                    << std::to_string((int)(100*t/config.t_final))+"%" << std::setw(8)
                    << iteration << std::setw(9)
                    << t  << std::setw(10)
                    << dt_old  << std::setw(9)
                    << vis_print << std::setw(12)
                    << relax_iter << std::setw(12)
                    << relax_res << std::setw(9)
                    << n_particles
                    << "\r";
          std::cout.flush();
      }

      if (n_particles == 0){
          if (config.verbose)
              std::cout << "No more particles.\n";
          break;
      }
  }
  auto t2 = high_resolution_clock::now();
  duration<double> ms_double = t2 - t1;
  std::cout << "\n\nExecution time: " << ms_double.count() << " s\n";

  return 0;
}
