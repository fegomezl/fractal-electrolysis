#include "header.h"

template<typename T>
void print_array(const Config &config, const std::vector<T> &array, const std::string name = "results/data.dat"); 

int main (int argc, char **argv){

    Config config;

    std::vector<bool> boundary(config.N);
    std::vector<double> phi(config.N);
    std::vector<double> particles;

    initialization(config, boundary, particles, phi);

    relaxation(config, boundary, phi);
    auto gradient = get_gradient(config, phi);

    print_array(config, phi, "results/phi.txt");
    print_array(config, boundary, "results/boundary.txt");
    print_array(config, gradient[0], "results/gradientx.txt");
    print_array(config, gradient[1], "results/gradienty.txt");

    std::cout << "nx: " << config.nx << " ny: " << config.ny << "\n";
    std::cout << "N: " << config.N << "\n";

    return 0;
}

template<typename T>
void print_array(const Config &config, const std::vector<T> &array, const std::string name){   
    /****
     * Print values of array in each position.
     ****/ 
    std::ofstream fout;
    fout.open(name);
    for (int j = 0; j < config.ny; ++j)
        for (int i = 0; i < config.nx; ++i){   
            double x = i*config.lx;
            double y = j*config.ly;
            fout << x << "\t" << y << "\t" << array[i+config.nx*j] << "\n";
        }
    fout.close();
}
