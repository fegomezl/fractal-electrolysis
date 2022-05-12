#include "header.h"

int main (int argc, char **argv){

    Config config;
    std::vector<bool> domain(config.N);
    std::vector<double> particles;
    std::vector<double> phi(config.N);
    std::vector<double> electric_field_x(config.N,0.0);
    std::vector<double> electric_field_y(config.N,0.0);
    std::vector<std::vector<double>> electric_field = {electric_field_x, electric_field_y};

    initialization(config, domain, particles, phi, electric_field);

    print_fields(config, domain, phi, electric_field, "results/data/fields_"+std::to_string(0)+".dat");
    print_particles(config, particles, "results/data/particles_"+std::to_string(0)+".dat");

    Crandom random(config.seed);

    std::cout << std::left << std::setw(12)
              << "------------------------------------------------\n"
              << std::left << std::setw(12)
              << "Step" << std::setw(12)
              << "Time" << std::setw(12)
              << "Printed" << std::setw(12)
              << "Progress"
              << std::left << std::setw(12)
              << "\n------------------------------------------------\n";

    for (int ii = 1; ii <= config.iterations; ii++){

        double percentage = 100*ii/config.iterations;
        std::cout << std::left << std::setw(12)
             << ii << std::setw(12)
             << ii*config.dt  << std::setw(12)
             << ii/config.vis_iterations << std::setw(12)
             << std::to_string((int)percentage)+"%" << "\r";
        std::cout.flush();

        system_evolve(config, random, domain, particles, phi, electric_field);
        relaxation(config, domain, phi,0);
        get_electric_field(config, phi, electric_field);

        if (ii%config.vis_iterations == 0){
            print_fields(config, domain, phi, electric_field, "results/data/fields_"+std::to_string(ii/config.vis_iterations)+".dat");
            print_particles(config, particles, "results/data/particles_"+std::to_string(ii/config.vis_iterations)+".dat");
        }

        if (particles.size() == 0){
            std::cout << "No more particles.";
            break;
        }
    }

    /*using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    auto t1 = high_resolution_clock::now();
    relaxation(config, domain, phi, 1);
    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << "Execution time: " << ms_double.count() << "ms \n";*/

    return 0;
}
