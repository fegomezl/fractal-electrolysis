#include "header.h"

void relaxation(const Config &config, const std::vector<bool> &domain, std::vector<double> &phi){
    int i,j,iter;
    auto phi_new = phi;
    double R,TotalRes;

    //array[j][i] -> array[i+nx*j] j->y, i->x
    for (iter = 0; iter < config.max_iter_relax; ++iter)
    {
        TotalRes=-1.0;
        //No he encontrado diferencia significativa entre static y guided.
        //Guided parece ser mejor que static cuando el numero de procesos no es multiplo y en casos raros.
        //Solo hay 1 caso en el que static es mejor que guided
        #pragma omp parallel for private(j,i,R), reduction(max:TotalRes), schedule(guided)
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

void relaxation(const Config &config, const std::vector<bool> &domain, std::vector<double> &phi, bool verbose){
    int i,j,iter;
    auto phi_new = phi;
    double R,TotalRes;

    //array[j][i] -> array[i+nx*j] j->y, i->x
    for (iter = 0; iter < config.max_iter_relax; ++iter)
    {
        TotalRes=-1.0;
        #pragma omp parallel for private(j,i,R), reduction(max:TotalRes), schedule(guided)
        for(j = 1; j < config.ny-1; j++)
            for(i = 1; i < config.ny-1; i++) {
                R = domain[i+config.nx*j]*(4*phi[i+config.nx*j] - phi[i+config.nx*(j+1)] - phi[i+config.nx*(j-1)] - phi[i+1+config.nx*j] - phi[i-1+config.nx*j]);
                phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                TotalRes = std::max(TotalRes,std::abs(R));
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
    int i,j;
    double partial_x=0, partial_y=0;

    //#pragma omp parallel private(j,i,size,start,partial_x,partial_y), num_threads(cores)
    #pragma omp parallel for private(j,i,partial_x,partial_y), schedule(guided)
        for(j = 1; j < config.ny-1; j++)
            for(i = 1; i < config.nx-1; i++) {
                partial_x = (phi[i+1+config.nx*j]-phi[i-1+config.nx*j])/(2*config.lx); //center deriv
                partial_y = (phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)])/(2*config.ly); //center deriv
                electric_field[0][i+config.nx*j] = -partial_x;
                electric_field[1][i+config.nx*j] = -partial_y;
            }
}
