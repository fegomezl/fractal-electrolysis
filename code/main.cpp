#include "header.h"

int main (int argc, char **argv){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    auto t1 = high_resolution_clock::now();

    Config config;
    std::vector<bool> domain(config.N);
    std::vector<double> particles;
    std::vector<double> phi(config.N);
    std::vector<double> electric_field_x(config.N, 0.0);
    std::vector<double> electric_field_y(config.N, 0.0);
    std::vector<std::vector<double>> electric_field = {electric_field_x, electric_field_y};

    Crandom random(config.seed);
    int iteration = 0;
    double t = 0.;
    double dt = config.dt_init;
    double dt_old = dt;
    bool last = false;
    int vis_iteration = 0;
    int vis_steps = config.vis_steps_max;
    int vis_print = 0;
    int relax_iter = 0;
    double relax_res = 0.;
    int n_particles = particles.size();

    if (config.verbose){
        std::cout << std::left << std::setw(12)
                  << "-------------------------------------------------------------------------------------\n"
                  << std::left << std::setw(12)
                  << "Progress" << std::setw(12)
                  << "Step" << std::setw(12)
                  << "Time" << std::setw(12)
                  << "Dt" << std::setw(12)
                  << "Printed" << std::setw(12)
                  << "Iterations" << std::setw(15)
                  << "Converged" << std::setw(15)
                  << "Particles"
                  << std::left << std::setw(12)
                  << "\n-------------------------------------------------------------------------------------\n";
    }

    n_particles = initialization(config, domain, phi, electric_field, particles);

    {
        //Update corresponding time_step
        double EMax = config.V/(config.Rint*log(config.Rext/config.Rint));
        double a = config.diffusivity*std::pow(EMax/config.V_ref, 2);
        double b = EMax*config.l/config.V_ref;
        dt = std::min(config.dt_init, (1+b+std::sqrt(1+2*b))/a);
        dt_old = dt;
    }

    print(config, t, domain, phi, electric_field, particles);

    if (config.verbose){
        std::cout << std::left << std::setw(12)
                  << "0%" << std::setw(12)
                  << iteration << std::setw(12)
                  << t  << std::setw(12)
                  << dt_old << std::setw(12)
                  << vis_print << std::setw(12)
                  << relax_iter << std::setw(15)
                  << relax_res << std::setw(15)
                  << n_particles 
                  << "\r";
        std::cout.flush();
    }

    for (iteration = 1, vis_iteration = 1; !last; iteration++, vis_iteration++){

        //Update iteration parameters
        last = (t >= config.t_final - 1e-8*config.dt_init);
        dt = std::min(dt, config.t_final - t);

        //Perform a time step
        n_particles = system_evolve(config, dt, random, domain, phi, electric_field, particles);
        t += dt;

        //Update visualization steps
        vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

        //Calculate new electric field
        dt_old = dt;
        auto [relax_iter, relax_res] = relaxation(config, dt, domain, phi, electric_field);

        if (last || vis_steps <= vis_iteration){
            //Update parameters
            vis_iteration = 0;
            vis_print += 1;

            print(config, t, domain, phi, electric_field, particles, "results/data/data_"+std::to_string(vis_print));
        }

        if (config.verbose){
            std::cout << std::left << std::setw(12)
                      << std::to_string((int)(100*t/config.t_final))+"%" << std::setw(12)
                      << iteration << std::setw(12)
                      << t  << std::setw(12)
                      << dt_old  << std::setw(12)
                      << vis_print << std::setw(12)
                      << relax_iter << std::setw(15)
                      << relax_res << std::setw(15)
                      << n_particles 
                      << "\r";
            std::cout.flush();
        }

        if (n_particles == 0){
            if (config.verbose)
                std::cout << "No more particles.\n";
            break;
        }
    }

    auto t2 = high_resolution_clock::now();
    duration<double> ms_double = t2 - t1;
    std::cout << "\n\nExecution time: " << ms_double.count() << " s\n";

    return 0;
}
