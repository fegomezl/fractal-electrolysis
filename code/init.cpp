#include "header.h"

int int_pow(int base, unsigned int exponent){
  if (exponent == 0) return 1;
  if (exponent == 1) return base;
  
  int results = int_pow(base, exponent/2);
  if (exponent%2 == 0) return results*results;
  else return base*results*results;
}

Config::Config(){
    base = 3;
    refinements = 4;

    Lx = 11.;
    Ly = 11.;

    Rint = 1.;
    Rext = 5.;

    V = 1.;
    V_dis = 0.5;

    max_iter_relax = 1000;
    alpha_relax = 1.0;
    res_relax = 1e-6;

    nx = int_pow(base, refinements);
    ny = int_pow(base, refinements);
    N = nx*ny;
    lx = Lx/nx;
    ly = Ly/ny;
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
