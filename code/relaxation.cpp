#include "header.h"

std::tuple<int,double> relaxation(const Config &config, double &dt, const std::vector<bool> &domain, std::vector<double> &phi, std::vector<std::vector<double>> &electric_field){
    int i=0,j=0,iter=0;
    double R=0,TotalRes=0,EMax=0;
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
            break;
    
        phi = phi_new;
    }

    #pragma omp parallel for private(j,i) schedule(guided) reduction(max: EMax)
    for(j = 1; j < config.n-1; j++){
        for(i = 1; i < config.n-1; i++) {
            electric_field[0][i+config.n*j] = -(phi[i+1+config.n*j]-phi[i-1+config.n*j])/(2*config.l); //center deriv
            electric_field[1][i+config.n*j] = -(phi[i+config.n*(j+1)]-phi[i+config.n*(j-1)])/(2*config.l); //center deriv
            double ELocal = std::hypot(electric_field[0][i+config.n*j], electric_field[1][i+config.n*j]);
            EMax = EMax > ELocal ? EMax : ELocal;
        }
    }

    //Update corresponding time_step
    double a = config.diffusivity*std::pow(EMax/config.V_ref, 2);
    double b = EMax*config.l/config.V_ref;
    dt = std::min(config.dt_init, (1+b+std::sqrt(1+2*b))/a);

    return {iter,TotalRes};
}

std::tuple<int,double> relaxation(const Config &config, double &dt, const std::vector<bool> &domain, std::vector<double> &phi, std::vector<std::vector<double>> &electric_field, bool verbose){
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
        std::cout << "Relaxation doesnt converge after " << iter << " steps. Residue: "<< TotalRes <<"\n";

    #pragma omp parallel for private(j,i) schedule(guided)
    for(j = 1; j < config.n-1; j++)
        for(i = 1; i < config.n-1; i++) {
            electric_field[0][i+config.n*j] = -(phi[i+1+config.n*j]-phi[i-1+config.n*j])/(2*config.l); //center deriv
            electric_field[1][i+config.n*j] = -(phi[i+config.n*(j+1)]-phi[i+config.n*(j-1)])/(2*config.l); //center deriv
        }

    return {iter,TotalRes};
}
