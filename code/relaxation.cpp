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


vector<vector3D> get_gradient(Config config, vector<double> &phi)
{   
    //size of the grid
    int i,j,t;
    vector <vector3D> grad_field(config.N);
    vector3D grad; grad.load(0,0,0);
    double partial_x,partial_y;

    for (i = 0; i < config.nx; ++i)
    {
        for (j = 1; j < config.ny-1; ++j){   
            //Left colum
            if(i ==0){
                partial_x = (phi[1+config.nx*j]-phi[config.nx*j])/(config.lx); //forward deriv
                partial_y = (phi[config.nx*(j+1)]-phi[config.nx*(j-1)])/(2*config.ly); //center deriv
                grad.load(partial_x,partial_y,0);
                grad_field[config.nx*j] = grad;
            }
            //Right colum
            else if(i==config.ny-1){
                partial_x = (phi[i+config.nx*j]-phi[i-1+config.nx*j])/(config.lx); //backward deriv
                partial_y = (phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)])/(2*config.ly); //center deriv
                grad.load(partial_x,partial_y,0);
                grad_field[i+config.nx*j] = grad;
            }
            //Center
            else{
                partial_x = (phi[i+1+config.nx*j]-phi[i-1+config.nx*j])/(2*config.lx); //center deriv
                partial_y = (phi[i+config.nx*(j+1)]-phi[i+config.nx*(j-1)])/(2*config.ly); //center deriv
                grad.load(partial_x,partial_y,0);
                grad_field[i+config.nx*j] = grad;
            }
        }

        //left corners
        if(i == 0){
            //bottom j=0
            partial_x = (phi[1]-phi[0])/(config.lx); //forward deriv
            partial_y = (phi[config.nx]-phi[0])/(config.ly); //forward deriv
            grad.load(partial_x,partial_y,0);
            grad_field[0] = grad;
            //top
            t = config.ny-1; //j
            partial_x = (phi[1+config.nx*t]-phi[config.nx*t])/(config.lx); //forward deriv
            partial_y = (phi[config.nx*t]-phi[config.nx*(t-1)])/(config.ly); //backward deriv
            grad.load(partial_x,partial_y,0);
            grad_field[config.nx*t] = grad;
        }
        //right corners
        else if (i == config.nx-1){
            //bottom j=0
            partial_x = (phi[i]-phi[i-1])/(config.lx); //backward deriv
            partial_y = (phi[i+config.nx]-phi[i])/(config.ly); //forward deriv
            grad.load(partial_x,partial_y,0);
            grad_field[i] = grad;
            //top
            t = config.ny-1; //j
            partial_x = (phi[i+config.nx*t]-phi[i-1+config.nx*t])/(config.lx); //backward deriv
            partial_y = (phi[i+config.nx*t]-phi[i+config.nx*(t-1)])/(config.ly); //backward deriv
            grad.load(partial_x,partial_y,0);
            grad_field[i+config.nx*t] = grad;
        }
        else{
            //Bottom row j=0
            partial_x = (phi[i+1]-phi[i-1])/(2*config.lx); //center deriv
            partial_y = (phi[i+config.nx]-phi[i])/(config.ly); //forward deriv
            grad.load(partial_x,partial_y,0);
            grad_field[i] = grad;
            //Top row
            t = config.ny-1; //j
            partial_x = (phi[i+1+config.nx*t]-phi[i-1+config.nx*t])/(2*config.lx); //center deriv
            partial_y = (phi[i+config.nx*t]-phi[i+config.nx*(t-1)])/(config.ly); //backward deriv
            grad.load(partial_x,partial_y,0);
            grad_field[i+config.nx*t] = grad;
        }
    }
    return grad_field;
}