#include "header.h"

void initialization(Config config, vector<double> &phi, vector<int> &boundary, vector<int> &dissociation);

template<typename T>
void print_array(Config config, vector<T> array, string name="results/data.dat");

int main (int argc, char **argv){

    Config config;

    vector<double> phi(config.N);
    vector<int> boundary(config.N);
    vector<int> dissociation(config.N);

    initialization(config, phi, boundary, dissociation);

    print_array(config, phi, "results/phi.txt");
    print_array(config, boundary, "results/boundary.txt");
    print_array(config, dissociation, "results/dissociation.txt");

    cout << "Number of elements: " << config.N << ".\n";

    return 0;
}

void initialization(Config config, vector<double> &phi, vector<int> &boundary, vector<int> &dissociation){

    /****
     * Initialization of electric potentiall, boundary 
     * conditions and dissociation probability.
     ****/ 
    for(int ii = 0; ii < config.N; ii++){
        double x = (ii%config.nx-(config.nx-1)/2)*config.lx;
        double y = (ii/config.nx-(config.ny-1)/2)*config.ly;
        double r = sqrt(x*x+y*y);

        if (r <= config.Rint) {
            phi[ii] = 0.;
            boundary[ii] = 0;
            dissociation[ii] = 0;
        } else if (r >= config.Rext) {
            phi[ii] = 1.;
            boundary[ii] = 0;
            dissociation[ii] = 0;
        } else {
            phi[ii] = config.V*log(r/config.Rint)/log(config.Rext/config.Rint); 
            boundary[ii] = 1;
            dissociation[ii] = (phi[ii] > config.V_dis) ? 1 : 0;
        }
    }
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
