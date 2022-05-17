#include "header.h"

void print_fields(const Config &config, const std::vector<bool> &domain, const std::vector<double> &phi, std::vector<std::vector<double>> &electric_field, const std::string name){
    
    /****
     * Print values of field in each position.
     ****/ 

    std::vector<int> print_range(config.nproc+1, config.N/config.nproc);
    print_range[0] = 0;

    for (long unsigned int ii = 0; ii < config.N%config.nproc; ii++)
        print_range[ii+1] += 1;

    for (long unsigned int ii = 1; ii < print_range.size(); ii++)
        print_range[ii] += print_range[ii-1];

    #pragma omp parallel
    {
        int pid = omp_get_thread_num();

        std::ofstream fout;
        fout.open(name+std::to_string(pid)+".dat");

        //internal values
        for(long unsigned int ii = print_range[pid]; ii < print_range[pid+1]; ii++)
            fout << domain[ii] << "\t" << phi[ii] << "\t" << electric_field[0][ii] << "\t" << electric_field[1][ii] << "\n";

        fout.close();
    }
}

void print_particles(const Config &config, const std::vector<double> &particles, const std::string name){ 
    
    /****
     * Print position of particles.
     ****/ 

    std::vector<int> print_range(config.nproc+1, (particles.size()/2)/config.nproc);
    print_range[0] = 0;

    for (long unsigned int ii = 0; ii < particles.size()%config.nproc; ii++)
        print_range[ii+1] += 1;

    for (long unsigned int ii = 1; ii < print_range.size(); ii++)
        print_range[ii] += print_range[ii-1];

    #pragma omp parallel
    {
        int pid = omp_get_thread_num();

        std::ofstream fout;
        fout.open(name+std::to_string(pid)+".dat");

        for (long unsigned int ii = print_range[pid]; ii < print_range[pid+1]; ii++)
            fout << particles[2*ii] << "\t" << particles[2*ii+1] << "\n";

        fout.close();
    }
}
