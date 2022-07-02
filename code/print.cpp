#include "header.h"

void print(const Config &config, double t, const std::vector<bool> &domain, const std::vector<double> &phi, const std::vector<std::vector<double>> &electric_field, const std::vector<double> &particles, const std::vector<int> &density, const std::string folder){

    std::filesystem::create_directory(folder);

    //Print time
    std::ofstream print_time;
    print_time.open(folder+"/time.txt", std::ios::trunc);
    print_time << t;
    print_time.close();
    
    /****
     * Print values of field in each position.
     ****/ 

    std::vector<int> print_fields_range(config.nproc+1, config.N/config.nproc);
    print_fields_range[0] = 0;

    for (int ii = 0; ii < config.N%config.nproc; ii++)
        print_fields_range[ii+1] += 1;

    for (long unsigned int ii = 1; ii < print_fields_range.size(); ii++)
        print_fields_range[ii] += print_fields_range[ii-1];

    std::vector<int> print_particles_range(config.nproc+1, (particles.size()/2)/config.nproc);
    print_particles_range[0] = 0;

    for (long unsigned int ii = 0; ii < particles.size()%config.nproc; ii++)
        print_particles_range[ii+1] += 1;

    for (long unsigned int ii = 1; ii < print_particles_range.size(); ii++)
        print_particles_range[ii] += print_particles_range[ii-1];

    #pragma omp parallel
    {
        int pid = omp_get_thread_num();

        std::ofstream print_fields;
        print_fields.open(folder+"/fields_"+std::to_string(pid)+".dat");

        for(int ii = print_fields_range[pid]; ii < print_fields_range[pid+1]; ii++)
            print_fields << domain[ii] << "\t" << phi[ii] << "\t" << electric_field[0][ii] << "\t" << electric_field[1][ii] << "\t" << density[ii] << "\n";

        print_fields.close();

        //Print only the fractal bitmap
        print_fields.open(folder+"/bit_map_"+std::to_string(pid)+".dat");
        for(int ii = print_fields_range[pid]; ii < print_fields_range[pid+1]; ii++){
        		if (phi[ii]==config.V)
            		print_fields << 1 << "\n";
            	else
            		print_fields << domain[ii] << "\n";
            }
        print_fields.close();

        std::ofstream print_particles;
        print_particles.open(folder+"/particles_"+std::to_string(pid)+".dat");

        for (int ii = print_particles_range[pid]; ii < print_particles_range[pid+1]; ii++)
            print_particles << particles[2*ii] << "\t" << particles[2*ii+1] << "\n";

        print_particles.close();
    }
}
