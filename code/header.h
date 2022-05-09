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

class vector3D{
 private:
  double X, Y, Z;

 public:
  //Initialize the vector
  void load(double x0, double y0, double z0);
  //Get the components
  double x(void){return X;};
  double y(void){return Y;};
  double z(void){return Z;};
  //Show the vector
  void show(void);
  //-------------------------
  //Vectorial operators
  //-------------------------
  //Equal
  void operator= (vector3D v2);
  //Sum
  vector3D operator+ (vector3D v2);
  void operator+=(vector3D v2);
  //Substraction
  vector3D operator- (vector3D v2);
  void operator-=(vector3D v2);
  //Scalar multiplication
  vector3D operator* (double a);
  void operator*=(double a);
  friend vector3D operator* (double a,vector3D v1); 
  //Scalar division
  vector3D operator/ (double a);
  void operator/=(double a);
  //Dot product
  double operator* (vector3D v2);
  //Cross product
  vector3D operator^ (vector3D v2);
  //Norm operations
  double norm2(void);    
  double norm(void);
  //Angle between two vectors
  friend double angle(vector3D v1, vector3D v2);
};

extern void initialization(Config config, vector<double> &phi, vector<int> &boundary, vector<int> &dissociation);

extern void relaxation(Config config, vector<double> &phi, vector<int> boundary, bool verbose=false);
extern vector<vector3D> get_gradient(Config config, vector<double> &phi);