#include "header.h"

template<typename T>
void print_field(const Config &config, const std::vector<T> &field, const std::string name = "results/data.dat"); 
void print_particles(const Config &config, const std::vector<double> &particles, const std::string name = "results/data.dat"); 

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

template<typename T>
void print_field(const Config &config, const std::vector<T> &field, const std::string name){   
    
    /****
     * Print values of field in each position.
     ****/ 

    std::ofstream fout;
    fout.open(name);

    //Internal values
    for(int ii = 0; ii < config.N; ii++)
        fout << field[ii] << "\n";

    fout.close();
}

void print_particles(const Config &config, const std::vector<double> &particles, const std::string name){ 
    
    /****
     * Print position of particles.
     ****/ 

    std::ofstream fout;
    fout.open(name);

    for (long unsigned int ii = 0; ii < particles.size()/2; ii++)
        fout << particles[2*ii] << "\t" << particles[2*ii+1] << "\n";

    fout.close();
}
