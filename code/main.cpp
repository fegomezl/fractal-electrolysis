#include "header.h"

int main (int argc, char **argv){

    Config config;
    std::vector<bool> domain(config.N);
    std::vector<double> particles;
    std::vector<double> phi(config.N);
    std::vector<double> electric_field_x(config.N);
    std::vector<double> electric_field_y(config.N);
    std::vector<std::vector<double>> electric_field = {electric_field_x, electric_field_y};

    initialization(config, domain, particles, phi);
    relaxation(config, domain, phi);
    get_electric_field(config, phi, electric_field);

    print_field(config, domain, "results/domain.dat");
    print_particles(config, particles, "results/particles.dat");

    print_field(config, phi, "results/phi.dat");
    print_field(config, electric_field[0], "results/electric_field_x.dat");
    print_field(config, electric_field[1], "results/electric_field_y.dat");

    std::cout << "nx: " << config.nx << " ny: " << config.ny << "\n";
    std::cout << "N: " << config.N << "\n";

    return 0;
}
