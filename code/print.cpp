#include "header.h"

void print(const Config &config, double t, const std::vector<bool> &domain, const std::vector<double> &phi, const std::vector<std::vector<double>> &electric_field, const std::vector<int> &density, const std::string folder){

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

    #pragma omp parallel
    {
        int pid = omp_get_thread_num();

        std::ofstream print_fields;
        print_fields.open(folder+"/fields_"+std::to_string(pid)+".dat");

        for(int ii = print_fields_range[pid]; ii < print_fields_range[pid+1]; ii++)
            print_fields << domain[ii] << "," << phi[ii] << "," << electric_field[0][ii] << "," << electric_field[1][ii] << "," << density[ii] << "\n";

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
    }
}
