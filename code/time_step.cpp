#include "header.h"

void system_evolve(const Config &config, Crandom &random, std::vector<bool> &domain, std::vector<double> &particles, std::vector<double> &phi, const std::vector<std::vector<double>> &electric_field){
    /****
     * Move particles according to the Smoluchowski Diffusion Equation.
     *
     * Diagram for estimation of electric field
     *      
     *  I01=(x0,y1)   I11=(x1,y1)
     *      o---------------o
     *      |               |
     *      |               |
     *      |       o       |
     *      |     (x,y)     |
     *      |               |
     *      |               |
     *      o---------------o
     *  I00=(x0,y0)   I10=(x1,y0)
     *
     ****/ 

    auto particles_old = particles;

    for(long unsigned int ii = 0; ii < particles.size()/2; ii++){
        double x = particles[2*ii]/config.lx;
        double y = particles[2*ii+1]/config.ly;
        int x0 = floor(x);
        int x1 = x0 + 1;
        int y0 = floor(y);
        int y1 = y0 + 1;

        double Ex00 = 0., Ex01 = 0., Ex10 = 0., Ex11 = 0.;
        double Ey00 = 0., Ey01 = 0., Ey10 = 0., Ey11 = 0.;

        if (abs(x0) <= (config.nx - 1)/2 && abs(y0) <= (config.ny - 1)/2){
            Ex00 = electric_field[0][x0+(config.nx-1)/2 + (y0+(config.ny-1)/2)*config.ny];
            Ey00 = electric_field[1][x0+(config.nx-1)/2 + (y0+(config.ny-1)/2)*config.ny];
        }
        if (abs(x0) <= (config.nx - 1)/2 && abs(y1) <= (config.ny - 1)/2){
            Ex01 = electric_field[0][x0+(config.nx-1)/2 + (y1+(config.ny-1)/2)*config.ny];
            Ey01 = electric_field[1][x0+(config.nx-1)/2 + (y1+(config.ny-1)/2)*config.ny];
        }
        if (abs(x1) <= (config.nx - 1)/2 && abs(y0) <= (config.ny - 1)/2){
            Ex10 = electric_field[0][x1+(config.nx-1)/2 + (y0+(config.ny-1)/2)*config.ny];
            Ey10 = electric_field[1][x1+(config.nx-1)/2 + (y0+(config.ny-1)/2)*config.ny];
        }
        if (abs(x1) <= (config.nx - 1)/2 && abs(y1) <= (config.ny - 1)/2){
            Ex11 = electric_field[0][x1+(config.nx-1)/2 + (y1+(config.ny-1)/2)*config.ny];
            Ey11 = electric_field[1][x1+(config.nx-1)/2 + (y1+(config.ny-1)/2)*config.ny];
        }

        double Ex = Ex00*(x1-x)*(y1-y) + Ex01*(x1-x)*(y-y0) + Ex10*(x-x0)*(y1-y) + Ex11*(x-x0)*(y-y0);
        double Ey = Ey00*(x1-x)*(y1-y) + Ey01*(x1-x)*(y-y0) + Ey10*(x-x0)*(y1-y) + Ey11*(x-x0)*(y-y0);

        particles[2*ii]   += config.mu*Ex + config.sigma*random.gauss(0., 1.);
        particles[2*ii+1] += config.mu*Ey + config.sigma*random.gauss(0., 1.);
    }

    int k = 0;
    for(long unsigned int ii = 0; ii < particles.size()/2; ii++){
        double x = particles[2*(ii-k)]/config.lx;
        double y = particles[2*(ii-k)+1]/config.ly;
        int x0 = floor(x);
        int x1 = x0 + 1;
        int y0 = floor(y);
        int y1 = y0 + 1;

        std::vector<int> neighbors;
        if (abs(x0) <= (config.nx - 1)/2 && abs(y0) <= (config.ny - 1)/2)
            neighbors.push_back(x0+(config.nx-1)/2 + (y0+(config.ny-1)/2)*config.ny);
        if (abs(x0) <= (config.nx - 1)/2 && abs(y1) <= (config.ny - 1)/2)
            neighbors.push_back(x0+(config.nx-1)/2 + (y1+(config.ny-1)/2)*config.ny);
        if (abs(x1) <= (config.nx - 1)/2 && abs(y0) <= (config.ny - 1)/2)
            neighbors.push_back(x1+(config.nx-1)/2 + (y0+(config.ny-1)/2)*config.ny);
        if (abs(x1) <= (config.nx - 1)/2 && abs(y1) <= (config.ny - 1)/2)
            neighbors.push_back(x1+(config.nx-1)/2 + (y1+(config.ny-1)/2)*config.ny);        

        bool liquid = 1;
        double V_new = 1.;
        for (auto &ii : neighbors){
            liquid &= domain[ii];
            V_new *= phi[ii];
        }

        if (!liquid){
            if (V_new == 0) {
                for (auto &ii : neighbors){
                    domain[ii] = 0;
                    phi[ii] = 0.;
                }
                particles.erase(particles.begin()+2*(ii-k), particles.begin()+2*(ii-k)+1);
                k += 1;
            } else {
                particles[2*(ii-k)] = particles_old[2*ii];
                particles[2*(ii-k)+1] = particles_old[2*ii+1];
            }
        }
    }
}
