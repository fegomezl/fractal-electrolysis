#include "header.h"

template<typename T>
void print_field(const Config &config, const std::vector<T> &field, const std::string name = "results/data.dat"); 
void print_particles(const Config &config, const std::vector<double> &particles, const std::string name = "results/data.dat"); 

int main (int argc, char **argv){

    Config config;
    std::vector<bool> boundary(config.N);
    std::vector<double> phi(config.N);
    std::vector<double> particles;

    initialization(config, boundary, particles, phi);
    relaxation(config, boundary, phi);
    auto gradient = get_gradient(config, phi);

    print_field(config, boundary, "results/boundary.dat");
    print_particles(config, particles, "results/particles.dat");

    print_field(config, phi, "results/phi.dat");
    print_field(config, gradient[0], "results/gradientx.dat");
    print_field(config, gradient[1], "results/gradienty.dat");

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
        fout << particles[ii] << "\t" << particles[ii+1] << "\n";

    fout.close();
}
