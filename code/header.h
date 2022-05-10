#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct Config{
    Config();
    int base;
    int refinements;

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
    double V_dis;

    int max_iter_relax;
    double alpha_relax;
    double res_relax;

    double molar_volume;
    double molarity;
    int M;
};

void initialization(const Config &config, std::vector<bool> &boundary, std::vector<double> &particles, std::vector<double> &phi);
void relaxation(const Config &config, const std::vector<bool> &boundary, std::vector<double> &phi, const bool verbose=false);
std::vector<std::vector<double>> get_gradient(const Config &config, const std::vector<double> &phi);
