#include "header.h"

template<typename T>
void print_array(Config config, vector<T> array, string name = "results/data.dat"); 

int main (int argc, char **argv){

    Config config;

    vector<double> phi(config.N);
    vector<int> boundary(config.N);
    vector<int> dissociation(config.N);

    initialization(config, phi, boundary, dissociation);

    relaxation(config, phi, boundary);
    auto gradient = get_gradient(config,phi);

    print_array(config, phi, "results/phi.txt");
    print_array(config, boundary, "results/boundary.txt");
    print_array(config, dissociation, "results/dissociation.txt");

    std::cout << "nx: " << config.nx << " ny: " << config.ny << "\n";

    return 0;
}

template<typename T>
void print_array(Config config, vector<T> array, string name){   

    /****
     * Print values of array in each position.
     ****/ 
    ofstream fout;
    fout.open(name);
    for(int ii = 0; ii < config.N; ii++){
        double x = (ii%config.nx-(config.nx-1)/2)*config.lx;
        double y = (ii/config.nx-(config.ny-1)/2)*config.ly;
        double r = sqrt(x*x+y*y);

        fout << x << "\t" << y << "\t" << array[ii] << "\n";
    }
    fout.close();
}
