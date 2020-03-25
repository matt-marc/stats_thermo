#ifndef DVS_H
#define DVS_H

#include <math.h>
#include <iostream>
#include <vector>

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
        double expon = exp(-0.5 * (rho / p) * (v - u) * (v - u));
        double front = pow((rho / (2 * M_PI * p)), 0.5);
        return rho * front * expon;
    }
};


class DiscreteVelocityScheme {
   public:
    DiscreteVelocityScheme(double cells, double lower_bound, double upper_bound);

    void setDensityInRange(double min, double max, density_function eq);

    void setDensityInRange(double min, double max);
    void setVelocitySpace(double number, double min, double max);

    void time_march_to(double df, double dt, double tau);

    friend std::ostream &operator<<(std::ostream &out, const DiscreteVelocityScheme &dsv);

    void writeF(double x, string filename);

    void testFuntions();


   private:
    double cells;
    double lower_bound;
    double upper_bound;
    double dx;
    vector<double> x_pos;

    vector<double> rho();
    vector<double> p();
    vector<double> q();
    vector<double> u();

    vector<double> F_minus_half(int index);
    vector<double> F_plus_half(int index);
    vector<double> F_flux(int index, double dx, double dt);
};

#endif