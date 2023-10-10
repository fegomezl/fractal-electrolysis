#pragma once
#include <chrono>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/IterativeLinearSolvers>

struct Config {
    Config();

    int nproc;
    bool verbose;

    double t_final;
    double dt_init;
    int vis_steps_max;

    int n;
    int N;

    int seed;
    double relax_alpha;
    int relax_max_iter;
    double relax_res;

    double L;
    double l;
    double Rint;
    double Rext;

    double particle_proportion;
    double diffusivity;
    double V_ref;

    double V;
    double sigma;
    double mu;

    double E_cte;
    int m;
};

class Crandom {
    unsigned long long u, v, w;

public:
    Crandom(unsigned long long j) {
        v = 4101842887655102017LL; w = 1;
        u = j ^ v; int64();
        v = u; int64();
        w = v; int64();
    }
    unsigned long long int64() {
        u = u * 2862933555777941757LL + 7046029254386353087LL;
        v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
        w = 4294957665U * (w & 0xffffffff) + (w >> 32);
        unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
        return (x + v) ^ w;
    }
    double r() {
        return 5.42101086242752217E-20 * int64();
    }
    unsigned int int32() {
        return (unsigned int) int64();
    };
    double exponencial(float tau) {
        return -tau * log(r());
    }
    double gauss(float mu, float sigma) {
        return sigma * sqrt(-2 * log(r())) * cos(2 * M_PI * r()) + mu;
    }
};

double initialization(const Config &config, Eigen::SparseMatrix<double, Eigen::RowMajor> &Domain, Eigen::VectorXd &phi, std::vector<Eigen::VectorXd> &electric_field, std::vector<double> &particles, Eigen::VectorXd &density);
void set_up_finite_laplacian(const Config &config, Eigen::SparseMatrix<double, Eigen::RowMajor> &Laplacian);
void set_up_finite_partial_dx(const Config &config, Eigen::SparseMatrix<double, Eigen::RowMajor> &Partial_dx);
void set_up_finite_partial_dy(const Config &config, Eigen::SparseMatrix<double, Eigen::RowMajor> &Partial_dy);

//std::tuple<int, double> relaxation(const Config &config, double &dt, const std::vector<bool> &domain, std::vector<double> &phi, std::vector<std::vector<double>> &electric_field, const std::vector<int> &density);
//std::tuple<int, double> relaxation(const Config &config, double &dt, const std::vector<bool> &domain, std::vector<double> &phi, std::vector<std::vector<double>> &electric_field, const std::vector<int> &density, const bool verbose);
void compute_electric_potential(Eigen::VectorXd &phi, const Eigen::SparseMatrix<double, Eigen::RowMajor> &Laplacian, const Eigen::SparseMatrix<double, Eigen::RowMajor> &Domain, const Eigen::SparseMatrix<double, Eigen::RowMajor> &I);
void compute_electric_field(std::vector<Eigen::VectorXd> &electric_field, const Eigen::VectorXd &phi, const Eigen::SparseMatrix<double, Eigen::RowMajor> &Partial_dx, const Eigen::SparseMatrix<double, Eigen::RowMajor> &Partial_dy);
 
void print(const Config &config, const double t, const Eigen::SparseMatrix<double, Eigen::RowMajor> &Domain, const Eigen::VectorXd &phi, const std::vector<Eigen::VectorXd> &electric_field, const Eigen::VectorXd &density, const std::string folder = "results/data/data_0");

double system_evolve(const Config &config, const double dt, Crandom &random, Eigen::SparseMatrix<double, Eigen::RowMajor> &Domain, Eigen::VectorXd &phi, const std::vector<Eigen::VectorXd> &electric_field, std::vector<double> &particles, Eigen::VectorXd &density);

