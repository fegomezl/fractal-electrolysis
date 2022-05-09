#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "random64.h"

using namespace std;

int int_pow(int base, unsigned int exponent){
  if (exponent == 0) return 1;
  if (exponent == 1) return base;
  
  int results = int_pow(base, exponent/2);
  if (exponent%2 == 0) return results*results;
  else return base*results*results;
}

struct Config{

    Config(){
        nx = int_pow(base, refinements);
        ny = int_pow(base, refinements);
        N = nx*ny;
        lx = Lx/nx;
        ly = Ly/ny;
    }

    int base = 3;
    int refinements = 4;

    int nx;
    int ny;
    int N;

    double Lx = 11.;
    double Ly = 11.;
    double lx;
    double ly;

    double Rint = 1.;
    double Rext = 5.;

    double V = 1.;
    double V_dis = 0.5;
};
