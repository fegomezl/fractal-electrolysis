#include "header.h"

void relaxation(Config config, vector<double> &phi, vector<int> boundary, bool verbose){
    int i,j,iter,t;
    auto phi_new = phi;
    double R=0.0;
    double TotalRes=0.0;


    for (iter = 0; iter < config.max_iter_relax; ++iter)
    {
        TotalRes=0.0;
        for (i = 0; i < config.nx; ++i)
        {
            for (j = 1; j < config.ny-1; ++j){   
                //Left colum
                if(i ==0){
                    R = boundary[config.nx*j]*(4*phi[config.nx*j]-phi[config.nx*(j+1)]-phi[config.nx*(j-1)]-phi[1+config.nx*j]-config.V);
                    phi_new[config.nx*j] = phi[config.nx*j]-config.alpha_relax*R*0.25;
                    TotalRes+=abs(R);
                }
                //Right colum
                else if(i==config.ny-1){
                    R = boundary[i+config.nx*j]*(4*phi[i+config.nx*j]-phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)]-phi[(i-1)+config.nx*j]-config.V);
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

            //left corners
            if(i == 0){
                //bottom j=0
                R = boundary[0]*(4*phi[0]-phi[config.nx]-phi[1]-2*config.V);
                phi_new[0] = phi[0]-config.alpha_relax*R*0.25;
                TotalRes+=abs(R);
                //top
                t = config.ny-1; //j
                R = boundary[config.nx*t]*(4*phi[config.nx*t]-phi[config.nx*(t-1)]-phi[1+config.nx*t]-2*config.V);
                phi_new[config.nx*t] = phi[config.nx*t]-config.alpha_relax*R*0.25;
                TotalRes+=abs(R);
            }
            //right corners
            else if (i == config.nx-1){
                //bottom j=0
                R = boundary[i]*(4*phi[i]-phi[i+config.nx]-phi[i-1]-2*config.V);
                phi_new[i] = phi[i]-config.alpha_relax*R*0.25;
                TotalRes+=abs(R);
                //top
                t = config.ny-1; //j
                R = boundary[i+config.nx*t]*(4*phi[i+config.nx*t]-phi[i+config.nx*(t-1)]-phi[(i-1)+config.nx*t]-2*config.V);
                phi_new[i+config.nx*t] = phi[i+config.nx*t]-config.alpha_relax*R*0.25;
                TotalRes+=abs(R);
            }
            else{
                //Bottom row j=0
                R = boundary[i]*(4*phi[i]-phi[i+config.nx]-phi[i+1]-phi[i-1]-config.V);
                phi_new[i] = phi[i]-config.alpha_relax*R*0.25;
                TotalRes+=abs(R);
                //Top row
                t = config.ny-1; //j
                R = boundary[i+config.nx*t]*(4*phi[i+config.nx*t]-phi[i+config.nx*(t-1)]-phi[(i+1)+config.nx*t]-phi[(i-1)+config.nx*t]-config.V);
                phi_new[i+config.nx*t] = phi[i+config.nx*t]-config.alpha_relax*R*0.25;
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


vector<vector<double>> get_gradient(Config config, vector<double> &phi)
{   
    //size of the grid
    int i,j,t;
    vector<double> grad_field_x(config.N);
    vector<double> grad_field_y(config.N);
    double partial_x,partial_y;

    for (i = 0; i < config.nx; ++i)
    {
        for (j = 1; j < config.ny-1; ++j){   
            //Left colum
            if(i ==0){
                partial_x = (phi[1+config.nx*j]-phi[config.nx*j])/(config.lx); //forward deriv
                partial_y = (phi[config.nx*(j+1)]-phi[config.nx*(j-1)])/(2*config.ly); //center deriv
                grad_field_x[config.nx*j] = partial_x;
                grad_field_y[config.nx*j] = partial_y;
            }
            //Right colum
            else if(i==config.ny-1){
                partial_x = (phi[i+config.nx*j]-phi[i-1+config.nx*j])/(config.lx); //backward deriv
                partial_y = (phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)])/(2*config.ly); //center deriv
                grad_field_x[i+config.nx*j] = partial_x;
                grad_field_y[i+config.nx*j] = partial_y;
            }
            //Center
            else{
                partial_x = (phi[i+1+config.nx*j]-phi[i-1+config.nx*j])/(2*config.lx); //center deriv
                partial_y = (phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)])/(2*config.ly); //center deriv
                grad_field_x[i+config.nx*j] = partial_x;
                grad_field_y[i+config.nx*j] = partial_y;
            }
        }

        //left corners
        if(i == 0){
            //bottom j=0
            partial_x = (phi[1]-phi[0])/(config.lx); //forward deriv
            partial_y = (phi[config.nx]-phi[0])/(config.ly); //forward deriv
            grad_field_x[0] = partial_x;
            grad_field_y[0] = partial_y;
            //top
            t = config.ny-1; //j
            partial_x = (phi[1+config.nx*t]-phi[config.nx*t])/(config.lx); //forward deriv
            partial_y = (phi[config.nx*t]-phi[config.nx*(t-1)])/(config.ly); //backward deriv
            grad_field_x[config.nx*t] = partial_x;
            grad_field_y[config.nx*t] = partial_y;
        }
        //right corners
        else if (i == config.nx-1){
            //bottom j=0
            partial_x = (phi[i]-phi[i-1])/(config.lx); //backward deriv
            partial_y = (phi[i+config.nx]-phi[i])/(config.ly); //forward deriv
            grad_field_x[i] = partial_x;
            grad_field_y[i] = partial_y;
            //top
            t = config.ny-1; //j
            partial_x = (phi[i+config.nx*t]-phi[i-1+config.nx*t])/(config.lx); //backward deriv
            partial_y = (phi[i+config.nx*t]-phi[i+config.nx*(t-1)])/(config.ly); //backward deriv
            grad_field_x[i+config.nx*t] = partial_x;
            grad_field_y[i+config.nx*t] = partial_y;
        }
        else{
            //Bottom row j=0
            partial_x = (phi[i+1]-phi[i-1])/(2*config.lx); //center deriv
            partial_y = (phi[i+config.nx]-phi[i])/(config.ly); //forward deriv
            grad_field_x[i] = partial_x;
            grad_field_y[i] = partial_y;
            //Top row
            t = config.ny-1; //j
            partial_x = (phi[i+1+config.nx*t]-phi[i-1+config.nx*t])/(2*config.lx); //center deriv
            partial_y = (phi[i+config.nx*t]-phi[i+config.nx*(t-1)])/(config.ly); //backward deriv
            grad_field_x[i+config.nx*t] = partial_x;
            grad_field_y[i+config.nx*t] = partial_y;
        }
    }
    vector<vector<double>> grad_field = {grad_field_x,grad_field_x};
    return grad_field;
}