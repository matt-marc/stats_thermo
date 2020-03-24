#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;


struct density_function {
    double rho;
    double u;
    double p;
    density_function(double _rho, double _u, double _p) {
        rho = _rho;
        u = _u;
        p = _p;
    }
    double operator()(double v) {
        double expon = exp(-0.5 * rho / p * (v - u) * (v - u));
        double front = pow((rho / (2 * M_PI * p)), 0.5);
        return front * expon;
    }
};

class Domain {
   public:
    Domain()=default;
    Domain(double num, double min, double max);
    Domain(vector<double> dis, vector<double> vels);
    vector<double> num_dis;
    vector<double> vel_space;
    

    void setNumberDensity(density_function eq);

    double rho();
    double u();
    double p();
    double q();

    Domain operator+(const Domain& d);


    double min;
    double max;
    double dv;

};

#endif