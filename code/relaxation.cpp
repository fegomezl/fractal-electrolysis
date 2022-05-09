#include "header.h"

void relaxation(Config config, vector<double> &phi, vector<int> boundary, bool verbose){
    int i,j,iter;
    auto phi_new = phi;
    double R=0.0;
    double TotalRes=0.0;


    for (iter = 0; iter < config.max_iter_relax; ++iter)
    {
        TotalRes=0.0;
        

        for (i = 0; i < config.nx; ++i)
        {
            for (j = 0; j < config.ny; ++j){   
                R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)]-phi[(i+1)+config.nx*j]-phi[(i-1)+config.nx*j]);
                phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                TotalRes+=abs(R);

                //Left colum
                if(i ==0){
                    R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)]-phi[(i+1)+config.nx*j]-config.V);
                    phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                    TotalRes+=abs(R);

                    //Corners
                    if(j==0){
                        R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j+1)]-phi[(i+1)+config.nx*j]-2*config.V);
                        phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                        TotalRes+=abs(R);
                    }
                    if(j==config.ny-1){
                        R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j-1)]-phi[(i+1)+config.nx*j]-2*config.V);
                        phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                        TotalRes+=abs(R);
                    }
                }
                //Right colum
                if(i==config.ny-1){
                    R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)]-phi[(i-1)+config.nx*j]-config.V);
                    phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                    TotalRes+=abs(R);

                    //Corners
                    if(j==0)
                    {
                        R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j+1)]-phi[(i-1)+config.nx*j]-2*config.V);
                        phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                        TotalRes+=abs(R);
                    }
                    if(j==config.ny-1)
                    {
                        R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j-1)]-phi[(i-1)+config.nx*j]-2*config.V);
                        phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                        TotalRes+=abs(R);
                    }
                }
            }

            //Bottom row
            j = 0;
            R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j+1)]-phi[(i+1)+config.nx*j]-phi[(i-1)+config.nx*j]-config.V);
            phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
            TotalRes+=abs(R);
            //Top row
            j = config.ny-1;
            R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j-1)]-phi[(i+1)+config.nx*j]-phi[(i-1)+config.nx*j]-config.V);
            phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
            TotalRes+=abs(R);
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
