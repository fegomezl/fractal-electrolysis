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
    refinements = 6;

    Lx = 16.;
    Ly = 16.;

    Rint = 1.;
    Rext = 6.;

    V = 1.;
    V_dis = 0.5;

    max_iter_relax = 10000;
    alpha_relax = 1.0;
    res_relax = 1e-8;

    molar_volume = 0.04;
    molarity = 1.;

    nx = int_pow(base, refinements); 
    ny = int_pow(base, refinements); 
    N = nx*ny;
    lx = Lx/nx; 
    ly = Ly/ny;
    particle_proportion = molarity*molar_volume;
}

void initialization(const Config &config, std::vector<bool> &boundary, std::vector<double> &particles, std::vector<double> &phi){

    /****
     * Initialization of electric potentiall, boundary 
     * conditions and dissociation probability.
     ****/
    std::vector<int> dissociation;    
    for(int ii = 0; ii < config.N; ii++){
        double x = (ii%config.nx-(config.nx-1)/2)*config.lx;
        double y = (ii/config.nx-(config.ny-1)/2)*config.ly;
        double r = hypot(x, y);

        if (r <= config.Rint) {
            phi[ii] = 0.;
            boundary[ii] = 0;
        } else if (r >= config.Rext) {
            phi[ii] = config.V;
            boundary[ii] = 0;
        } else {
            phi[ii] = config.V*log(r/config.Rint)/log(config.Rext/config.Rint);
            boundary[ii] = 1;
            dissociation.push_back(ii);
        }
    }

    Crandom random(config.seed);
    int N_sites = dissociation.size();
    int N_particles = config.particle_proportion*dissociation.size();
    for(int ii = 0; ii < N_particles; ii++){
        int jj = (N_sites-1)*random.r();
        int kk = dissociation[jj];
        particles.push_back((kk%config.nx-(config.nx-1)/2)*config.lx);
        particles.push_back((kk/config.nx-(config.ny-1)/2)*config.ly);
    }
}
