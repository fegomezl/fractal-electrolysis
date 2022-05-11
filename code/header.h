#pragma once
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

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
    double t;
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

void initialization(const Config &config, std::vector<bool> &boundary, std::vector<double> &phi, std::vector<double> &particles);
void relaxation(const Config &config, const std::vector<bool> &boundary, std::vector<double> &phi, const bool verbose=false);
std::vector<std::vector<double>> get_gradient(const Config &config, const std::vector<double> &phi);
