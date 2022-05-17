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
    verbose = std::stoi(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    iterations = std::stoi(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    vis_iterations = std::stoi(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    n = std::stoi(line.erase(line.find('#'), line.size()));

    for (int ii = 0; ii < 2; ii++)
        std::getline(parameters,line);

    std::getline(parameters,line);
    seed = std::stoi(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    relax_alpha = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    relax_max_iter = (int)std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    relax_res = std::stod(line.erase(line.find('#'), line.size()));

    for (int ii = 0; ii < 2; ii++)
        std::getline(parameters,line);

    std::getline(parameters,line);
    dt = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    L = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    Rint = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    Rext = std::stod(line.erase(line.find('#'), line.size()));

    for (int ii = 0; ii < 2; ii++)
        std::getline(parameters,line);

    std::getline(parameters,line);
    double molarity = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    double molar_volume = std::stod(line.erase(line.find('#'), line.size()));

    for (int ii = 0; ii < 2; ii++)
        std::getline(parameters,line);

    std::getline(parameters,line);
    V = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    double T = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    double diffusivity = std::stod(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    double oxidation = std::stoi(line.erase(line.find('#'), line.size()));

    std::getline(parameters,line);
    double electro_boltzmann = std::stod(line.erase(line.find('#'), line.size()));

    parameters.close();

    /****
     * Set derived constants
     ****/

    N = n*n;
    l = L/n; 
    particle_proportion = molarity*molar_volume;
    sigma = sqrt(2*diffusivity*dt);
    mu = electro_boltzmann*oxidation*diffusivity*dt/(T+273.15); 
}

double initialization(const Config &config, std::vector<bool> &domain, std::vector<double> &phi, std::vector<std::vector<double>> &electric_field, std::vector<double> &particles){

    /****
     * Initialization of electric potential, domain 
     * conditions and dissociation probability.
     ****/
    std::vector<int> dissociation;    
    for(int ii = 0; ii < config.N; ii++){
        double x = (ii%config.n-(config.n-1)/2)*config.l;
        double y = (ii/config.n-(config.n-1)/2)*config.l;
        double r = hypot(x, y);

        if (r <= config.Rint) {
            domain[ii] = 0;
            phi[ii] = 0.;
            electric_field[0][ii] = 0.;
            electric_field[1][ii] = 0.;
        } else if (r >= config.Rext) {
            domain[ii] = 0;
            phi[ii] = config.V;
            electric_field[0][ii] = 0.;
            electric_field[1][ii] = 0.;
        } else {
            domain[ii] = 1;
            phi[ii] = config.V*log(r/config.Rint)/log(config.Rext/config.Rint);
            electric_field[0][ii] = -config.V*x/(r*r*log(config.Rext/config.Rint));
            electric_field[1][ii] = -config.V*y/(r*r*log(config.Rext/config.Rint));
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
        particles.push_back((kk%config.n-(config.n-1)/2)*config.l);
        particles.push_back((kk/config.n-(config.n-1)/2)*config.l);
    }

    return particles.size();
}
