#include "header.h"

void benchmark(const Config &config, std::vector<bool> &domain, std::vector<double> &particles, std::vector<double> &phi, std::vector<std::vector<double>> &electric_field)
{
    //cronometer parameters
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    //system_evolve(config, random, domain, particles, phi, electric_field);
    //relaxation(config, domain, phi);
    //get_electric_field(config, phi, electric_field);
    std::string name = "results/scaling.txt";
    std::ofstream fout;
    fout.open(name);

    int samples = 100;
    int max_cores = 8;

    for (int i = 1; i <  max_cores+1; ++i)
    {   
        double t=0;
        double t_2=0;

        for (int j = 0; j < samples; ++j)
        {
            initialization(config, domain, particles, phi, electric_field);

            auto t1 = high_resolution_clock::now();
            relaxation(config, domain, phi, i);
            //get_electric_field(config, phi, electric_field);
            auto t2 = high_resolution_clock::now();

            duration<double, std::milli> ms_double = t2 - t1;
            double a = ms_double.count();
            t+=a;
            t_2+=a*a;
        }
        double average = t/(samples*1.0);
        double average2 = t_2/(samples*1.0);
        std::cout << "cores:\t" << i << "\t<t>:\t" << average << "\tvar(t):\t" << average2 - average*average << "\tms\n";
        fout << i << "," << average << "," << average2 - average*average << "\n";
    }
    fout.close();
}