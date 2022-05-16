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

    //benchmark(config, domain, particles, phi, electric_field);

    std::filesystem::create_directories("results/data/data_0");
    print_fields(config, domain, phi, electric_field, "results/data/data_0/fields_");
    print_particles(config, particles, "results/data/data_0/particles_");

    Crandom random(config.seed);
    int ii = 0;
    int printed = 0;
    int n_particles = particles.size();
    double t = 0.;
    double total_res = 0.;

    if (config.verbose){
        std::cout << std::left << std::setw(12)
                  << "----------------------------------------------------------------------\n"
                  << std::left << std::setw(12)
                  << "Progress" << std::setw(12)
                  << "Step" << std::setw(12)
                  << "Time" << std::setw(12)
                  << "Printed" << std::setw(12)
                  << "Converged" << std::setw(12)
                  << "Particles"
                  << std::left << std::setw(12)
                  << "\n----------------------------------------------------------------------\n";

        std::cout << std::left << std::setw(12)
                  << "0%" << std::setw(12)
                  << ii << std::setw(12)
                  << t  << std::setw(12)
                  << printed << std::setw(12)
                  << total_res << std::setw(12)
                  << n_particles 
                  << "\r";
        std::cout.flush();
    }

    for (ii = 1; ii <= config.iterations; ii++){

        n_particles = system_evolve(config, random, domain, particles, phi, electric_field);
        total_res = relaxation(config, domain, phi);
        get_electric_field(config, phi, electric_field);

        if (ii%config.vis_iterations == 0){
            printed += 1;
            std::filesystem::create_directories("results/data/data_"+std::to_string(printed));
            print_fields(config, domain, phi, electric_field, "results/data/data_"+std::to_string(printed)+"/fields_");
            print_particles(config, particles, "results/data/data_"+std::to_string(printed)+"/particles_");
        }

        if (config.verbose){
            std::cout << std::left << std::setw(12)
                      << std::to_string((int)100*ii/config.iterations)+"%" << std::setw(12)
                      << ii << std::setw(12)
                      << ii*config.dt  << std::setw(12)
                      << printed << std::setw(12)
                      << total_res << std::setw(12)
                      << n_particles 
                      << "\r";
            std::cout.flush();
        }

        if (n_particles == 0)
            break;
    }

    if ((ii-1)%config.vis_iterations != 0){
        printed += 1;
        print_fields(config, domain, phi, electric_field, "results/data/fields_"+std::to_string(printed)+".dat");
        print_particles(config, particles, "results/data/particles_"+std::to_string(printed)+".dat");

        if (config.verbose){
            std::cout << std::left << std::setw(12)
                      << std::to_string((int)100*ii/config.iterations)+"%" << std::setw(12)
                      << ii << std::setw(12)
                      << ii*config.dt  << std::setw(12)
                      << printed << std::setw(12)
                      << total_res << std::setw(12)
                      << n_particles 
                      << "\r";
            std::cout.flush();
        }
    }

    if (n_particles == 0 && config.verbose)
        std::cout << "No more particles.\n";

    return 0;
}
