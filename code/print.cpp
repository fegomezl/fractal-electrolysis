#include "header.h"

void print_field(const Config &config, const std::vector<bool> &field, const std::string name){   
    
    /****
     * print values of field in each position.
     ****/ 

    std::ofstream fout;
    fout.open(name);

    //internal values
    for(int ii = 0; ii < config.N; ii++)
        fout << field[ii] << "\n";

    fout.close();
}

void print_field(const Config &config, const std::vector<double> &field, const std::string name){   
    
    /****
     * print values of field in each position.
     ****/ 

    std::ofstream fout;
    fout.open(name);

    //internal values
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
