#include "header.h"

double system_evolve(const Config &config, const double dt, Crandom &random, Eigen::SparseMatrix<double, Eigen::RowMajor> &Domain, Eigen::VectorXd &phi, const std::vector<Eigen::VectorXd> &electric_field, std::vector<double> &particles, Eigen::VectorXd &density){
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

    double x=0., y=0.;
    int x0=0, x1=0, y0=0, y1=0;
    double Ex=0., Ey=0.;

    //Reset Density counts
    density.setZero();
    Eigen::VectorXd domain = Domain.diagonal();

    for (long unsigned int ii = 0; ii < particles.size()/2; ii++){

        x = particles[2*ii]/config.l;
        y = particles[2*ii+1]/config.l;
        x0 = floor(x);
        x1 = x0 + 1;
        y0 = floor(y);
        y1 = y0 + 1;

        Ex = Ey = 0.;            
        if (((1-2*std::signbit(x0))*(x0+x1) <= config.n - 2) && ((1-2*std::signbit(y0))*(y0+y1) <= config.n - 2)){
            Ex = (x1-x)*(y1-y)*electric_field[0][x0+(config.n-1)/2 + (y0+(config.n-1)/2)*config.n]
               + (x1-x)*(y-y0)*electric_field[0][x0+(config.n-1)/2 + (y1+(config.n-1)/2)*config.n]
               + (x-x0)*(y1-y)*electric_field[0][x1+(config.n-1)/2 + (y0+(config.n-1)/2)*config.n]
               + (x-x0)*(y-y0)*electric_field[0][x1+(config.n-1)/2 + (y1+(config.n-1)/2)*config.n];
            Ey = (x1-x)*(y1-y)*electric_field[1][x0+(config.n-1)/2 + (y0+(config.n-1)/2)*config.n]
               + (x1-x)*(y-y0)*electric_field[1][x0+(config.n-1)/2 + (y1+(config.n-1)/2)*config.n]
               + (x-x0)*(y1-y)*electric_field[1][x1+(config.n-1)/2 + (y0+(config.n-1)/2)*config.n]
               + (x-x0)*(y-y0)*electric_field[1][x1+(config.n-1)/2 + (y1+(config.n-1)/2)*config.n];
        }

        particles[2*ii]   += config.mu*dt*Ex + config.sigma*std::sqrt(dt)*random.gauss(0., 1.);
        particles[2*ii+1] += config.mu*dt*Ey + config.sigma*std::sqrt(dt)*random.gauss(0., 1.);

        //Density count
        density[x0+(config.n-1)/2 + (y0+(config.n-1)/2)*config.n]+=1;
    }

    bool liquid = 0;
    double V_new = 0.;
    int k = 0;
    std::vector<int> neighbors = {0, 0, 0, 0};

    for(long unsigned int ii = 0; ii < particles.size()/2; ii++){

        x = particles[2*(ii-k)]/config.l;
        y = particles[2*(ii-k)+1]/config.l;
        x0 = floor(x);
        x1 = x0 + 1;
        y0 = floor(y);
        y1 = y0 + 1;

        liquid = 1;
        V_new = 1.;
        if (((1-2*std::signbit(x0))*(x0+x1) <= config.n - 2) && ((1-2*std::signbit(y0))*(y0+y1) <= config.n - 2)){

            neighbors[0] = x0+(config.n-1)/2 + (y0+(config.n-1)/2)*config.n;
            neighbors[1] = x0+(config.n-1)/2 + (y1+(config.n-1)/2)*config.n;
            neighbors[2] = x1+(config.n-1)/2 + (y0+(config.n-1)/2)*config.n;
            neighbors[3] = x1+(config.n-1)/2 + (y1+(config.n-1)/2)*config.n;

            for (auto &jj : neighbors){
                liquid &= 1 - int(domain[jj]);
                V_new *= phi[jj];
            }

            if (!liquid){
                if (V_new == 0) {
                    for (auto &jj : neighbors){
                        domain[jj] = 1;
                        if (Domain.coeff(jj, jj) == 0)
                            Domain.insert(jj, jj) = 1;
                        phi[jj] = 0.0;
                    }
                    particles.erase(particles.begin()+2*(ii-k), particles.begin()+2*(ii-k)+1);
                    k += 1;
                } else {
                    particles[2*(ii-k)] = particles_old[2*ii];
                    particles[2*(ii-k)+1] = particles_old[2*ii+1];
                }
            }
        } else {
            particles[2*(ii-k)] = particles_old[2*ii];
            particles[2*(ii-k)+1] = particles_old[2*ii+1];
        }
    }

    return particles.size();
}