#pragma once
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <omp.h>
#include <future>

struct Config{
    Config();

    int nx;
    int ny;
    int N;

    double Lx;
    double Ly;
    double lx;
    double ly;

    double Rint;
    double Rext;

    double V;

    int max_iter_relax;
    double alpha_relax;
    double res_relax;

    double molar_volume;
    double molarity;
    double particle_proportion;

    int seed;
    
    double dt;
    int iterations;
    int vis_iterations;
    double diffusivity;
    int oxidation;
    double electro_boltzmann;
    double T;

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

void initialization(const Config &config, std::vector<bool> &domain, std::vector<double> &particles, std::vector<double> &phi, std::vector<std::vector<double>> &electric_field);

void relaxation(const Config &config, const std::vector<bool> &domain, std::vector<double> &phi);
void relaxation(const Config &config, const std::vector<bool> &domain, std::vector<double> &phi, const bool verbose);
void get_electric_field(const Config &config, const std::vector<double> &phi, std::vector<std::vector<double>> &electric_field);

void print_fields(const Config &config, const std::vector<bool> &domain, const std::vector<double> &phi, std::vector<std::vector<double>> &electric_field, const std::string name = "results/data/fields.dat"); 
void print_particles(const Config &config, const std::vector<double> &particles, const std::string name = "results/data/particles.dat"); 

void system_evolve(const Config &config, Crandom &random, std::vector<bool> &domain, std::vector<double> &particles, std::vector<double> &phi, const std::vector<std::vector<double>> &electric_field);
void benchmark(const Config &config, std::vector<bool> &domain, std::vector<double> &particles, std::vector<double> &phi, std::vector<std::vector<double>> &electric_field);