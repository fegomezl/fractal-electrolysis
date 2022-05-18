#include "header.h"

int main (int argc, char **argv){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    auto t1 = high_resolution_clock::now();

    Config config;
    std::vector<bool> domain(config.N);
    std::vector<double> particles;
    std::vector<double> phi(config.N);
    std::vector<double> electric_field_x(config.N,0.0);
    std::vector<double> electric_field_y(config.N,0.0);
    std::vector<std::vector<double>> electric_field = {electric_field_x, electric_field_y};

    Crandom random(config.seed);
    int ii = 0, iter = 0;
    int printed = 0;
    int n_particles = particles.size();
    double t = 0.;
    double total_res = 0.;

    if (config.verbose){
        std::cout << std::left << std::setw(12)
                  << "-------------------------------------------------------------------------------------\n"
                  << std::left << std::setw(12)
                  << "Progress" << std::setw(12)
                  << "Step" << std::setw(12)
                  << "Time" << std::setw(12)
                  << "Printed" << std::setw(12)
                  << "Iterations" << std::setw(15)
                  << "Converged" << std::setw(15)
                  << "Particles"
                  << std::left << std::setw(12)
                  << "\n-------------------------------------------------------------------------------------\n";
    }

    n_particles = initialization(config, domain, phi, electric_field, particles);
    print(config, domain, phi, electric_field, particles);

    if (config.verbose){
        std::cout << std::left << std::setw(12)
                  << "0%" << std::setw(12)
                  << ii << std::setw(12)
                  << t  << std::setw(12)
                  << printed << std::setw(12)
                  << iter << std::setw(15)
                  << total_res << std::setw(15)
                  << n_particles 
                  << "\r";
        std::cout.flush();
    }

    for (ii = 1; ii <= config.iterations; ii++){

        n_particles = system_evolve(config, random, domain, phi, electric_field, particles);
        auto [iter,total_res] = relaxation(config, domain, phi, electric_field);

        if (ii%config.vis_iterations == 0){
            printed += 1;
            print(config, domain, phi, electric_field, particles, "results/data/data_"+std::to_string(printed));
        }

        if (config.verbose){
            std::cout << std::left << std::setw(12)
                      << std::to_string((int)100*ii/config.iterations)+"%" << std::setw(12)
                      << ii << std::setw(12)
                      << ii*config.dt  << std::setw(12)
                      << printed << std::setw(12)
                      << iter << std::setw(15)
                      << total_res << std::setw(15)
                      << n_particles 
                      << "\r";
            std::cout.flush();
        }

        if (n_particles == 0)
            break;
    }

    if ((ii-1)%config.vis_iterations != 0){
        printed += 1;
        print(config, domain, phi, electric_field, particles, "results/data/data_"+std::to_string(printed));

        if (config.verbose){
            std::cout << std::left << std::setw(12)
                      << std::to_string((int)100*ii/config.iterations)+"%" << std::setw(12)
                      << ii << std::setw(12)
                      << ii*config.dt  << std::setw(12)
                      << printed << std::setw(12)
                      << iter << std::setw(15)
                      << total_res << std::setw(15)
                      << n_particles 
                      << "\r";
            std::cout.flush();
        }
    }

    if (n_particles == 0 && config.verbose)
        std::cout << "No more particles.\n";

    auto t2 = high_resolution_clock::now();
    duration<double> ms_double = t2 - t1;
    std::cout << "\n\nExecution time: " << ms_double.count() << " s\n";

    return 0;
}
