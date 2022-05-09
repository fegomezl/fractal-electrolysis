#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

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

extern void initialization(Config config, vector<double> &phi, vector<int> &boundary, vector<int> &dissociation);

extern void relaxation(Config config, vector<double> &phi, vector<int> boundary, bool verbose=false);
extern vector<vector<double>> get_gradient(Config config, vector<double> &phi);