#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <algorithm>

template<std::size_t SIZE>
void boundary_conditions(std::array<double, SIZE> & phi,std::array<int, SIZE> & boundary,int nx,int ny);
template<std::size_t SIZE>
void relaxation(std::array<double, SIZE> & phi,std::array<int, SIZE> boundary,int nx,int ny, int max_iter=1000, bool verbose=false, double alpha=1.0, double res=1e-6);
template<typename T, std::size_t SIZE>
void print_array(std::array<T, SIZE> array, int nx, int ny, std::string name="data.dat");

int main(int argc, char const *argv[])
{
    const int nx = 100;
    const int ny = 100;
    const int N = nx*ny;
    std::array<double, N> phi; for(auto &i : phi) i =  (rand() % 20);
    std::array<int, N> boundary; boundary.fill(1);//0 if a cell is a boundary

    //Impose boundary conditions
    boundary_conditions(phi,boundary,nx,ny);
    
    //Solve Laplace Equation
    relaxation(phi,boundary,nx,ny);

    //Save solution to file
    print_array(phi,nx,ny);

    return 0;
}

template<std::size_t SIZE>
void relaxation(std::array<double, SIZE> & phi,std::array<int, SIZE> boundary,int nx,int ny, int max_iter, bool verbose, double alpha, double res)
{
    int i,j, iter;
    auto phi_old = phi;
    auto phi_new = phi;
    double R=0.0;
    double TotalRes=0.0;

    for (iter = 0; iter < max_iter; ++iter)
    {
        TotalRes=0.0;
        //Center of Domain
        for (i = 1; i < nx-1; ++i)
            for (j = 1; j < ny-1; ++j){   
                R = boundary[i+nx*j]*(4*phi_old[i+nx*j]-phi_old[i+nx*(j+1)]-phi_old[i+nx*(j-1)]-phi_old[(i+1)+nx*j]-phi_old[(i-1)+nx*j]);
                phi_new[i+nx*j] = phi_old[i+nx*j]-alpha*R*0.25;
                TotalRes+=abs(R);
            }

        for (i = 1; i < nx-1; ++i)
        {
            //Bottom row
            j =0;
            R = boundary[i+nx*j]*(4*phi_old[i+nx*j]-phi_old[i+nx*(j+1)]-phi_old[(i+1)+nx*j]-phi_old[(i-1)+nx*j]);
            phi_new[i+nx*j] = phi_old[i+nx*j]-alpha*R*0.25;
            TotalRes+=abs(R);
            //Top row
            j = ny-1;
            R = boundary[i+nx*j]*(4*phi_old[i+nx*j]-phi_old[i+nx*(j-1)]-phi_old[(i+1)+nx*j]-phi_old[(i-1)+nx*j]);
            phi_new[i+nx*j] = phi_old[i+nx*j]-alpha*R*0.25;
            TotalRes+=abs(R);
        }

        for (j = 1; j < ny-1; ++j)
        {
            //Left colum
            i =0;
            R = boundary[i+nx*j]*(4*phi_old[i+nx*j]-phi_old[i+nx*(j+1)]-phi_old[i+nx*(j-1)]-phi_old[(i+1)+nx*j]);
            phi_new[i+nx*j] = phi_old[i+nx*j]-alpha*R*0.25;
            TotalRes+=abs(R);
            //Right colum
            i = nx-1;
            R = boundary[i+nx*j]*(4*phi_old[i+nx*j]-phi_old[i+nx*(j+1)]-phi_old[i+nx*(j-1)]-phi_old[(i-1)+nx*j]);
            phi_new[i+nx*j] = phi_old[i+nx*j]-alpha*R*0.25;
            TotalRes+=abs(R);
        }
        //Corners
        i=0; j=0;
        R = boundary[i+nx*j]*(4*phi_old[i+nx*j]-phi_old[i+nx*(j+1)]-phi_old[(i+1)+nx*j]);
        phi_new[i+nx*j] = phi_old[i+nx*j]-alpha*R*0.25;
        TotalRes+=abs(R);
        i=nx-1; j=0;
        R = boundary[i+nx*j]*(4*phi_old[i+nx*j]-phi_old[i+nx*(j+1)]-phi_old[(i-1)+nx*j]);
        phi_new[i+nx*j] = phi_old[i+nx*j]-alpha*R*0.25;
        TotalRes+=abs(R);
        i=0; j=ny-1;
        R = boundary[i+nx*j]*(4*phi_old[i+nx*j]-phi_old[i+nx*(j-1)]-phi_old[(i+1)+nx*j]);
        phi_new[i+nx*j] = phi_old[i+nx*j]-alpha*R*0.25;
        TotalRes+=abs(R);
        i=nx-1; j=ny-1;
        R = boundary[i+nx*j]*(4*phi_old[i+nx*j]-phi_old[i+nx*(j-1)]-phi_old[(i-1)+nx*j]);
        phi_new[i+nx*j] = phi_old[i+nx*j]-alpha*R*0.25;
        TotalRes+=abs(R);

        //Check if method Converged
        //TotalRes/=1.0*(nx*ny);
        if (TotalRes<res)
        {
            std::cout << "Relaxation converged after " << iter << " steps. Residue: "<< TotalRes <<"\n";
            break;
        }

        if (verbose)
            std::cout << "Iteration: " << iter << " Residue: " << TotalRes <<"\n";
    
        phi_old = phi_new;
    }
    if(iter==max_iter)
        std::cout << "Relaxation dint converge after " << iter << " steps. Residue: "<< TotalRes <<"\n";

    phi = phi_new;
}

template<std::size_t SIZE>
void boundary_conditions(std::array<double, SIZE> & phi,std::array<int, SIZE> & boundary,int nx,int ny){
    int i, j;
    for (i = 0; i < nx; ++i)
    {
        //Bottom row
        j =0;
        phi[i+j*nx]=20.0;
        boundary[i+j*nx] = 0;
        //Top row
        j = ny-1;
        phi[i+j*nx]=20.0;
        boundary[i+j*nx] =0;
    }
    for (j = 1; j < ny-1; ++j) //Dont re-writte boundaries in corners
    {
        //Left colum
        i =0;
        phi[i+j*nx]=20.0;
        boundary[i+j*nx] = 0;
        //Right colum
        i = nx-1;
        phi[i+j*nx]=20.0;
        boundary[i+j*nx] = 0;
    }
    for (i = int(2*nx/5); i < int(3*nx/5); ++i)
        for (j = int(2*ny/5); j < int(3*ny/5); ++j)
        {
            phi[i+j*nx]=-1.0;
            boundary[i+j*nx] = 0;
        }
}
template<typename T, std::size_t SIZE>
void print_array(std::array<T, SIZE> array, int nx, int ny, std::string name)
{   
    std::ofstream fout;
    fout.open(name);
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            fout << i << "\t" << j << "\t" << array[i+nx*j] << "\n";
        }
    }
    fout.close();
}