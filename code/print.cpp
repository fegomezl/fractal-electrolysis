#include "header.h"

void print_fields(const Config &config, const std::vector<bool> &domain, const std::vector<double> &phi, std::vector<std::vector<double>> &electric_field, const std::string name){
    
    /****
     * Print values of field in each position.
     ****/ 

    std::ofstream fout;
    fout.open(name);

    //internal values
    for(int ii = 0; ii < config.N; ii++)
        fout << domain[ii] << "\t" << phi[ii] << "\t" << electric_field[0][ii] << "\t" << electric_field[1][ii] << "\n";

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

