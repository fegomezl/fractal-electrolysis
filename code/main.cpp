#include "header.h"

int main (int argc, char **argv){
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    Config config;
    std::vector<bool> domain(config.N);
    std::vector<double> particles;
    std::vector<double> phi(config.N);
    std::vector<double> electric_field_x(config.N,0.0);
    std::vector<double> electric_field_y(config.N,0.0);
    std::vector<std::vector<double>> electric_field = {electric_field_x, electric_field_y};

    initialization(config, domain, particles, phi, electric_field);

    benchmark(config, domain, particles, phi, electric_field);

    return 0;

    /*//La impresion se hace en un thread independiente del programa.
    //Por lo que el programa no necesita esperar 
    //a que la impresion termine para seguirl.
    //No necesariamente en un core aparte, 
    //pero en performace es casi equivalente a no imprimir nada.
    auto aux1 = std::async(std::launch::async, [&config, &domain, &phi, &electric_field]{print_fields(config, domain, phi, electric_field, "results/data/fields_"+std::to_string(0)+".dat");});
    auto aux2 = std::async(std::launch::async, [&config, &particles]{print_particles(config, particles, "results/data/particles_"+std::to_string(0)+".dat");});

    Crandom random(config.seed);

    std::cout << std::left << std::setw(12)
              << "------------------------------------------------\n"
              << std::left << std::setw(12)
              << "Step" << std::setw(12)
              << "Time" << std::setw(12)
              << "Printed" << std::setw(12)
              << "Progress"
              << std::left << std::setw(12)
              << "\n------------------------------------------------\n";

    auto t1 = high_resolution_clock::now();
    for (int ii = 1; ii <= config.iterations; ii++){

        double percentage = 100*ii/config.iterations;
        std::cout << std::left << std::setw(12)
             << ii << std::setw(12)
             << ii*config.dt  << std::setw(12)
             << ii/config.vis_iterations << std::setw(12)
             << std::to_string((int)percentage)+"%" << "\r";
        std::cout.flush();

        system_evolve(config, random, domain, particles, phi, electric_field);
        relaxation(config, domain, phi);
        get_electric_field(config, phi, electric_field);

        if (ii%config.vis_iterations == 0){
            aux1.get();
            aux2.get();
            aux1 = std::async(std::launch::async, [&config, &domain, &phi, &electric_field,ii]{print_fields(config, domain, phi, electric_field, "results/data/fields_"+std::to_string(ii/config.vis_iterations)+".dat");});
            aux2 = std::async(std::launch::async, [&config, &particles,ii]{print_particles(config, particles, "results/data/particles_"+std::to_string(ii/config.vis_iterations)+".dat");});
        }

        if (particles.size() == 0){
            std::cout << "No more particles.";
            break;
        }
    }
    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << "Execution time: " << ms_double.count()/1000.0 << "ms \n";

    //Esperar a la impresion del ultimo frame 
    aux1.get();
    aux2.get();
    return 0;*/
}
