#include "header.h"

std::tuple<int,double> relaxation(const Config &config, double &dt, const std::vector<bool> &domain, std::vector<double> &phi, std::vector<std::vector<double>> &electric_field, const std::vector<int> &density){
    int i=0,j=0,iter=0;
    double R=0,TotalRes=0,EMax=0,sumx=0,sumy=0;
    auto phi_new = phi;
    std::vector<std::vector<double>> Grad_N = {electric_field[0], electric_field[1]};

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

    //Calculate Electric Field and Density Gradient
    #pragma omp parallel for private(j,i,sumx,sumy) schedule(guided) reduction(max: EMax)
    for(j = 1; j < config.n-1; j++){
        for(i = 1; i < config.n-1; i++) {
            electric_field[0][i+config.n*j] = -(phi[i+1+config.n*j]-phi[i-1+config.n*j])/(2*config.l); //center deriv
            electric_field[1][i+config.n*j] = -(phi[i+config.n*(j+1)]-phi[i+config.n*(j-1)])/(2*config.l); //center deriv
            Grad_N[0][i+config.n*j] = (density[i+1+config.n*j]-density[i-1+config.n*j])/2.0;
            Grad_N[1][i+config.n*j] = (density[i+config.n*(j+1)]-density[i+config.n*(j-1)])/2.0;
            sumx = 0;
            sumy = 0;
            for(int jj = 1; jj < config.n-1; jj++)
                for(int ii = 1; ii < config.n-1; ii++) {
                double r2 = (i-ii)*(i-ii)+(j-jj)*(j-jj);
                r2 = std::sqrt(r2)/(r2+1e-8);
                sumx+=Grad_N[0][ii+config.n*jj]*r2;
                sumy+=Grad_N[1][ii+config.n*jj]*r2;
            }

            electric_field[0][i+config.n*j]+=config.E_cte*sumx;
            electric_field[1][i+config.n*j]+=config.E_cte*sumy;

            EMax = std::max(EMax,std::hypot(electric_field[0][i+config.n*j], electric_field[1][i+config.n*j]));
        }
    }

    //Update corresponding time_step
    double a = config.diffusivity*std::pow(EMax/config.V_ref, 2);
    double b = EMax*config.l/(2*config.V_ref);
    dt = std::min(config.dt_init, (4+b+std::sqrt(16+8*b))/a);

    return {iter,TotalRes};
}

std::tuple<int,double> relaxation(const Config &config, double &dt, const std::vector<bool> &domain, std::vector<double> &phi, std::vector<std::vector<double>> &electric_field, const std::vector<int> &density, bool verbose){
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
