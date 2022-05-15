#include "header.h"

Config::Config(){

    /****
     * Read parameters of the system and 
     * stablish constants for the program.
     ****/

    std::string line;
    std::ifstream parameters("settings/parameters.txt");

    for (int ii = 0; ii < 14; ii++)
        std::getline(parameters,line);

    std::getline(parameters,line);
    nproc = std::stoi(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    iterations = std::stoi(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    vis_iterations = std::stoi(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    nx = ny = std::stoi(line.erase(line.find('#'), line.size()));

    for (int ii = 0; ii < 2; ii++)
        std::getline(parameters,line);

    std::getline(parameters,line);
    seed = std::stoi(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    alpha_relax = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    max_iter_relax = (int)std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    res_relax = std::stod(line.erase(line.find('#'), line.size()));

    for (int ii = 0; ii < 2; ii++)
        std::getline(parameters,line);

    std::getline(parameters,line);
    dt = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    Lx = Ly = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    Rint = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    Rext = std::stod(line.erase(line.find('#'), line.size()));

    for (int ii = 0; ii < 2; ii++)
        std::getline(parameters,line);

    std::getline(parameters,line);
    molarity = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    molar_volume = std::stod(line.erase(line.find('#'), line.size()));

    for (int ii = 0; ii < 2; ii++)
        std::getline(parameters,line);

    std::getline(parameters,line);
    V = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    T = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    diffusivity = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    oxidation = std::stoi(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    electro_boltzmann = std::stod(line.erase(line.find('#'), line.size()));

    parameters.close();

    N = nx*ny;
    lx = Lx/nx; 
    ly = Ly/ny;
    particle_proportion = molarity*molar_volume;

    T += 273.15;
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
