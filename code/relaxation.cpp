#include "header.h"

double relaxation(const Config &config, const std::vector<bool> &domain, std::vector<double> &phi){
    int i=0,j=0,iter=0;
    double R=0,TotalRes=0;
    auto phi_new = phi;

    //array[j][i] -> array[i+n*j] j->y, i->x
    for (iter = 0; iter < config.relax_max_iter; ++iter)
    {
        TotalRes=-1.0;
        //No he encontrado diferencia significativa entre static y guided.
        //Guided parece ser mejor que static cuando el numero de procesos no es multiplo y en casos raros.
        //Solo hay 1 caso en el que static es mejor que guided
        #pragma omp parallel for private(j,i,R) reduction(max:TotalRes) schedule(guided)
        for(j = 1; j < config.n-1; j++)
            for(i = 1; i < config.n-1; i++) {
                R = domain[i+config.n*j]*(4*phi[i+config.n*j] - phi[i+config.n*(j+1)] - phi[i+config.n*(j-1)] - phi[i+1+config.n*j] - phi[i-1+config.n*j]);
                phi_new[i+config.n*j] = phi[i+config.n*j]-config.relax_alpha*R*0.25;
                TotalRes = std::max(TotalRes,std::abs(R));
            }

        //Check if method Converged
        if (TotalRes<config.relax_res)
            break;
    
        phi = phi_new;
    }

    return TotalRes;
}

double relaxation(const Config &config, const std::vector<bool> &domain, std::vector<double> &phi, bool verbose){
    int i=0,j=0,iter=0;
    double R=0,TotalRes=0;
    auto phi_new = phi;

    //array[j][i] -> array[i+n*j] j->y, i->x
    for (iter = 0; iter < config.relax_max_iter; ++iter)
    {
        TotalRes=-1.0;
        #pragma omp parallel for private(j,i,R) reduction(max:TotalRes) schedule(guided)
        for(j = 1; j < config.n-1; j++)
            for(i = 1; i < config.n-1; i++) {
                R = domain[i+config.n*j]*(4*phi[i+config.n*j] - phi[i+config.n*(j+1)] - phi[i+config.n*(j-1)] - phi[i+1+config.n*j] - phi[i-1+config.n*j]);
                phi_new[i+config.n*j] = phi[i+config.n*j]-config.relax_alpha*R*0.25;
                TotalRes = std::max(TotalRes,std::abs(R));
            }

        //Check if method Converged
        if (TotalRes<config.relax_res)
        {
            std::cout << "Relaxation converged after " << iter << " steps. Residue: "<< TotalRes <<"\n";
            break;
        }
        if (verbose)
            std::cout << "Iteration: " << iter << " Residue: " << TotalRes <<"\n";
    
        phi = phi_new;
    }
    if(iter==config.relax_max_iter)
        std::cout << "Relaxation dint converge after " << iter << " steps. Residue: "<< TotalRes <<"\n";
    return TotalRes;
}

void get_electric_field(const Config &config, const std::vector<double> &phi, std::vector<std::vector<double>> &electric_field){
    //size of the grid
    int i=0,j=0;

    //#pragma omp parallel private(j,i,size,start,partial_x,partial_y), num_threads(cores)
    #pragma omp parallel for private(j,i) schedule(guided)
        for(j = 1; j < config.n-1; j++)
            for(i = 1; i < config.n-1; i++) {
                electric_field[0][i+config.n*j] = -(phi[i+1+config.n*j]-phi[i-1+config.n*j])/(2*config.l); //center deriv
                electric_field[1][i+config.n*j] = -(phi[i+config.n*(j+1)]-phi[i+config.n*(j-1)])/(2*config.l); //center deriv
            }
}
