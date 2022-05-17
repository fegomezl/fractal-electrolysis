#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <string>
#include <omp.h>
#include <vector>

#pragma omp declare reduction (merge : std::vector<int> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction(vec_mult : std::vector<double> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::multiplies<double>())) initializer(omp_priv = omp_orig)

struct Config{
    Config();

    int nproc;
    bool verbose;

    int iterations;
    int vis_iterations;

    int n;
    int N;

    int seed;
    double relax_alpha;
    int relax_max_iter;
    double relax_res;

    double dt;
    double L;
    double l;
    double Rint;
    double Rext;

    double particle_proportion;

    double V;
    double sigma;
    double mu;
};

class Crandom{
	unsigned long long u,v,w;
	
	public:
		Crandom(unsigned long long j){
            v=4101842887655102017LL; w=1;
            u = j ^ v; int64();
            v = u; int64();
            w = v; int64();
        }
		unsigned long long int64(){
            u = u * 2862933555777941757LL + 7046029254386353087LL;
            v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
            w = 4294957665U*(w & 0xffffffff) + (w >> 32);
            unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
            return (x + v) ^ w;
        }
		double r(){
            return 5.42101086242752217E-20 * int64();
        }
		unsigned int int32(){
            return (unsigned int) int64();
        };
		double exponencial(float tau){
            return -tau*log(r());
        }
		double gauss(float mu,float sigma){
            return sigma*sqrt(-2*log(r()))*cos(2*M_PI*r())+mu;
        }
};

double initialization(const Config &config, std::vector<bool> &domain, std::vector<double> &particles, std::vector<double> &phi, std::vector<std::vector<double>> &electric_field);

double relaxation(const Config &config, const std::vector<bool> &domain, std::vector<double> &phi);
double relaxation(const Config &config, const std::vector<bool> &domain, std::vector<double> &phi, const bool verbose);
void get_electric_field(const Config &config, const std::vector<double> &phi, std::vector<std::vector<double>> &electric_field);

void print_fields(const Config &config, const std::vector<bool> &domain, const std::vector<double> &phi, std::vector<std::vector<double>> &electric_field, const std::string name = "results/data/fields.dat"); 
void print_particles(const Config &config, const std::vector<double> &particles, const std::string name = "results/data/particles.dat"); 

double system_evolve(const Config &config, Crandom &random, std::vector<bool> &domain, std::vector<double> &particles, std::vector<double> &phi, const std::vector<std::vector<double>> &electric_field);

void benchmark(const Config &config, std::vector<bool> &domain, std::vector<double> &particles, std::vector<double> &phi, std::vector<std::vector<double>> &electric_field);

