#include "header.h"

int main (int argc, char **argv){

    Config config;
    std::vector<bool> domain(config.N);
    std::vector<double> particles;
    std::vector<double> phi(config.N);
    std::vector<double> electric_field_x(config.N);
    std::vector<double> electric_field_y(config.N);
    std::vector<std::vector<double>> electric_field = {electric_field_x, electric_field_y};

    initialization(config, domain, particles, phi, electric_field);

    print_fields(config, domain, phi, electric_field, "results/data/fields_"+std::to_string(0)+".dat");
    print_particles(config, particles, "results/data/particles_"+std::to_string(0)+".dat");

    for (int ii = 1; ii <= config.iterations; ii++){
        if (ii%config.vis_iterations == 0){
            print_fields(config, domain, phi, electric_field, "results/data/fields_"+std::to_string(ii/config.vis_iterations)+".dat");
            print_particles(config, particles, "results/data/particles_"+std::to_string(ii/config.vis_iterations)+".dat");
        }
    } 

    return 0;
}
