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
    void setVelocitySpace(double number, double min, double max);

    void setDensityInRange(double min, double max, density_function eq);

    void time_march_to(double df, double dt, double tau);

    void write_U(string filename);
    void write_F(string filename, double x);

    friend std::ostream &operator<<(std::ostream &out, const DiscreteVelocityScheme &dsv);

   private:
    double lower_bound;
    double upper_bound;
    double dx;
    double dv;

    vector<double> x_pos;
    vector<double> vel_space;

    //could probably of have used eigen.h
    vector<vector<double>> U;

    vector<double> rho();
    vector<double> p();
    vector<double> q();
    vector<double> u();

    double rho(int index);
    double u(int index);

    double rho(vector<double>);
    double u(vector<double>);
    double p(vector<double>);

    double F_flux(int index, int vel_index, double dx, double dt);

    void init_U(void);
};

#endif