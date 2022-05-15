#include "header.h"

Config::Config(){

    Lx = Ly = 9.;
    nx = ny = 6081; 

    Rint = 0.15;
    Rext = 4.5;

    V = 1.;

    seed = 50;

    max_iter_relax = 5;
    alpha_relax = 1.0;
    res_relax = 1e-6;

    molar_volume = 0.04;
    molarity = 1.;

    N = nx*ny;
    lx = Lx/nx; 
    ly = Ly/ny;
    particle_proportion = molarity*molar_volume;

    double T_celsius = 20;

    iterations = 3600*2;
    vis_iterations = 120;
    dt = 1.;
    T = T_celsius + 273.15;
    diffusivity = 0.733E-5;
    oxidation = 2;
    electro_boltzmann = 1.160451E+4;

    sigma = sqrt(2*diffusivity*dt);
    mu = electro_boltzmann*oxidation*diffusivity*dt/T; 
}

void initialization(const Config &config, std::vector<bool> &domain, std::vector<double> &particles, std::vector<double> &phi, std::vector<std::vector<double>> &electric_field){

    /****
     * Initialization of electric potential, domain 
     * conditions and dissociation probability.
     ****/
    std::vector<int> dissociation;    
    for(int ii = 0; ii < config.N; ii++){
        double x = (ii%config.nx-(config.nx-1)/2)*config.lx;
        double y = (ii/config.nx-(config.ny-1)/2)*config.ly;
        double r = hypot(x, y);

        if (r <= config.Rint) {
            phi[ii] = 0.;
            domain[ii] = 0;
        } else if (r >= config.Rext) {
            phi[ii] = config.V;
            domain[ii] = 0;
        } else {
            phi[ii] = config.V*log(r/config.Rint)/log(config.Rext/config.Rint);
            domain[ii] = 1;
            dissociation.push_back(ii);
        }
    }

    /****
     * Random setting of particles in domain. The distribution is uniform
     * along the free cells.
     ****/
    Crandom random(config.seed);
    int N_sites = dissociation.size();
    int N_particles = config.particle_proportion*dissociation.size();
    for(int ii = 0; ii < N_particles; ii++){
        int jj = (N_sites-1)*random.r();
        int kk = dissociation[jj];
        particles.push_back((kk%config.nx-(config.nx-1)/2)*config.lx);
        particles.push_back((kk/config.nx-(config.ny-1)/2)*config.ly);
    }

    /****
     * Calculate initial electric field
     ****/
    get_electric_field(config, phi, electric_field);
}
