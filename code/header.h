#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//using namespace std;

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
};

void initialization(Config config, std::vector<double> &phi, std::vector<int> &boundary, std::vector<int> &dissociation);
void relaxation(Config config, std::vector<double> &phi, const std::vector<int> &boundary, bool verbose=false);
std::vector<std::vector<double>> get_gradient(Config config, std::vector<double> &phi);