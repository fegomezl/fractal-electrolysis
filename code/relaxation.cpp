#include "header.h"

void relaxation(const Config &config, const std::vector<bool> &domain, std::vector<double> &phi, int cores){
    int i,j,iter;
    auto phi_new = phi;
    double R,TotalRes;

    //array[j][i] -> array[i+nx*j] j->y, i->x
    for (iter = 0; iter < config.max_iter_relax; ++iter)
    {
        TotalRes=-1.0;

        #pragma omp parallel for private(j,i,R), reduction(max:TotalRes), schedule(static), num_threads(cores)
        for(j = 1; j < config.ny-1; j++)
            for(i = 1; i < config.ny-1; i++) {
                R = domain[i+config.nx*j]*(4*phi[i+config.nx*j] - phi[i+config.nx*(j+1)] - phi[i+config.nx*(j-1)] - phi[i+1+config.nx*j] - phi[i-1+config.nx*j]);
                phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                TotalRes = std::max(TotalRes,std::abs(R));
            }

        //Check if method Converged
        if (TotalRes<config.res_relax)
            break;
    
        phi = phi_new;
    }
}

void relaxation(const Config &config, const std::vector<bool> &domain, std::vector<double> &phi, bool verbose, int cores){
    int i,j,iter,size,start;
    auto phi_new = phi;
    double R,TotalRes;

    //array[j][i] -> array[i+nx*j] j->y, i->x
    for (iter = 0; iter < config.max_iter_relax; ++iter)
    {
        TotalRes=-1.0;
        #pragma omp parallel private(j,i,R,size,start), reduction(max:TotalRes), num_threads(cores)
        {
            size = (config.ny-2)/omp_get_num_threads();
            start = 1 + omp_get_thread_num()*size;

            for(j = 1; j < config.ny-1; j++)
                for(i = start; i < start+size; i++) {
                    R = domain[i+config.nx*j]*(4*phi[i+config.nx*j] - phi[i+config.nx*(j+1)] - phi[i+config.nx*(j-1)] - phi[i+1+config.nx*j] - phi[i-1+config.nx*j]);
                    phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                    TotalRes = std::max(TotalRes,std::abs(R));
                }
        }
        //Check if method Converged
        if (TotalRes<config.res_relax)
        {
            std::cout << "Relaxation converged after " << iter << " steps. Residue: "<< TotalRes <<"\n";
            break;
        }
        if (verbose)
            std::cout << "Iteration: " << iter << " Residue: " << TotalRes <<"\n";
    
        phi = phi_new;
    }
    if(iter==config.max_iter_relax)
        std::cout << "Relaxation dint converge after " << iter << " steps. Residue: "<< TotalRes <<"\n";
}

void get_electric_field(const Config &config, const std::vector<double> &phi, std::vector<std::vector<double>> &electric_field){
    //size of the grid
    int i,j,size,start;
    double partial_x=0, partial_y=0;

    #pragma omp parallel private(j,i,size,start,partial_x,partial_y)
    {
        size = (config.ny-2)/omp_get_num_threads();
        start = 1 + omp_get_thread_num()*size;

        for(j = 1; j < config.ny-1; j++)
            for(i = start; i < start+size; i++) {
                partial_x = (phi[i+1+config.nx*j]-phi[i-1+config.nx*j])/(2*config.lx); //center deriv
                partial_y = (phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)])/(2*config.ly); //center deriv
                electric_field[0][i+config.nx*j] = -partial_x;
                electric_field[1][i+config.nx*j] = -partial_y;
            }
    }
}
