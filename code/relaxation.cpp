#include "header.h"

void relaxation(const Config &config, const std::vector<bool> &boundary, std::vector<double> &phi, bool verbose){
    int i,j,iter;
    std::vector<double> phi_new(config.N);
    phi_new = phi;
    double R=0.0;
    double TotalRes=0.0;

    //array[j][i] -> array[i+nx*j] j->y, i->x
    for (iter = 0; iter < config.max_iter_relax; ++iter)
    {
        TotalRes=0.0;
        for (j = 0; j < config.ny; ++j)
        {
            for (i = 1; i < config.nx-1; ++i){  
                //Bottom row
                if(j==0){
                    R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j+1)]-phi[(i+1)+config.nx*j]-phi[(i-1)+config.nx*j]-config.V);
                    phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                    TotalRes+=abs(R);
                }
                //Top row
                else if(j==config.ny-1){
                    R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j-1)]-phi[(i+1)+config.nx*j]-phi[(i-1)+config.nx*j]-config.V);
                    phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                    TotalRes+=abs(R);
                }
                //Center
                else{
                    R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)]-phi[(i+1)+config.nx*j]-phi[(i-1)+config.nx*j]);
                    phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                    TotalRes+=abs(R);
                }
            }
            //Bottom Corners
            if(j==0){
                //left
                i=0;
                R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j+1)]-phi[(i+1)+config.nx*j]-2*config.V);
                phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                TotalRes+=abs(R);
                //right
                i=config.nx-1;
                R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j+1)]-phi[(i-1)+config.nx*j]-2*config.V);
                phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                TotalRes+=abs(R);
            }
            //Top Corners
            else if(j==config.ny-1){
                //left
                i=0;
                R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j-1)]-phi[(i+1)+config.nx*j]-2*config.V);
                phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                TotalRes+=abs(R);
                //right
                i=config.nx-1;
                R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j-1)]-phi[(i-1)+config.nx*j]-2*config.V);
                phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                TotalRes+=abs(R);
            }
            //Left and Right colums
            else{
                //Left
                i=0;
                R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)]-phi[(i+1)+config.nx*j]-config.V);
                phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                TotalRes+=abs(R);
                //Right
                i=config.nx-1;
                R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)]-phi[(i-1)+config.nx*j]-config.V);
                phi_new[i+config.nx*j] = phi[i+config.nx*j]-config.alpha_relax*R*0.25;
                TotalRes+=abs(R);
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

std::vector<std::vector<double>> get_gradient(const Config &config, const std::vector<double> &phi){   
    //size of the grid
    int i,j;
    std::vector<double> grad_field_x(config.N);
    std::vector<double> grad_field_y(config.N);
    double partial_x=0, partial_y=0;

    for (j = 0; j < config.ny; ++j)
    {
        for (i = 1; i < config.nx-1; ++i){  
            //Bottom row
            if(j==0){
                partial_x = (phi[i+1+config.nx*j]-phi[i-1+config.nx*j])/(2*config.lx); //center deriv
                partial_y = (phi[i+config.nx*(j+1)]-phi[i+config.nx*j])/(config.ly); //forward deriv
                grad_field_x[i+config.nx*j] = -partial_x;
                grad_field_y[i+config.nx*j] = -partial_y;
            }
            //Top row
            else if(j==config.ny-1){
                partial_x = (phi[i+1+config.nx*j]-phi[i-1+config.nx*j])/(2*config.lx); //center deriv
                partial_y = (phi[i+config.nx*j]-phi[i+config.nx*(j-1)])/(config.ly); //backward deriv
                grad_field_x[i+config.nx*j] = -partial_x;
                grad_field_y[i+config.nx*j] = -partial_y;
            }
            //Center
            else{
                partial_x = (phi[i+1+config.nx*j]-phi[i-1+config.nx*j])/(2*config.lx); //center deriv
                partial_y = (phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)])/(2*config.ly); //center deriv
                grad_field_x[i+config.nx*j] = -partial_x;
                grad_field_y[i+config.nx*j] = -partial_y;
            }
        }
        //Bottom Corners
        if(j==0){
            //left
            i=0;
            partial_x = (phi[i+1+config.nx*j]-phi[i+config.nx*j])/(config.lx); //forward deriv
            partial_y = (phi[i+config.nx*(j+1)]-phi[i+config.nx*j])/(config.ly); //forward deriv
            grad_field_x[i+config.nx*j] = -partial_x;
            grad_field_y[i+config.nx*j] = -partial_y;
            //right
            i=config.nx-1;
            partial_x = (phi[i+config.nx*j]-phi[i-1+config.nx*j])/(config.lx); //backward deriv
            partial_y = (phi[i+config.nx*(j+1)]-phi[i+config.nx*j])/(config.ly); //forward deriv
            grad_field_x[i+config.nx*j] = -partial_x;
            grad_field_y[i+config.nx*j] = -partial_y;
        }
        //Top Corners
        else if(j==config.ny-1){
            //left
            i=0;
            partial_x = (phi[i+1+config.nx*j]-phi[i+config.nx*j])/(config.lx); //forward deriv
            partial_y = (phi[i+config.nx*j]-phi[i+config.nx*(j-1)])/(config.ly); //backward deriv
            grad_field_x[i+config.nx*j] = -partial_x;
            grad_field_y[i+config.nx*j] = -partial_y;
            //right
            i=config.nx-1;
            partial_x = (phi[i+config.nx*j]-phi[i-1+config.nx*j])/(config.lx); //backward deriv
            partial_y = (phi[i+config.nx*j]-phi[i+config.nx*(j-1)])/(config.ly); //backward deriv
            grad_field_x[i+config.nx*j] = -partial_x;
            grad_field_y[i+config.nx*j] = -partial_y;
        }
        //Left and Right colums
        else{
            //Left
            i=0;
            partial_x = (phi[i+1+config.nx*j]-phi[i+config.nx*j])/(config.lx); //forward deriv
            partial_y = (phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)])/(2*config.ly); //center deriv
            grad_field_x[i+config.nx*j] = -partial_x;
            grad_field_y[i+config.nx*j] = -partial_y;
            //Right
            i=config.nx-1;
            partial_x = (phi[i+config.nx*j]-phi[i-1+config.nx*j])/(config.lx); //backward deriv
            partial_y = (phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)])/(2*config.ly); //center deriv
            grad_field_x[i+config.nx*j] = -partial_x;
            grad_field_y[i+config.nx*j] = -partial_y;
        }
    }
    std::vector<std::vector<double>> grad_field = {grad_field_x, grad_field_y};
    return grad_field;
}
