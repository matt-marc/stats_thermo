#ifndef DVS_H
#define DVS_H

#include <iostream>
#include <vector>
#include <math.h>

#include "Domain.h"



class DiscreteVelocityScheme {
   public:
    DiscreteVelocityScheme(double cells, double lower_bound, double upper_bound);

    void setDensityInRange(double min, double max, density_function eq);

    void setDensityInRange(double min, double max);
    void setVelocitySpace(double number, double min, double max);

    

    void time_march_to(double df, double dt, double tau);

    friend std::ostream &operator<<(std::ostream &out, const DiscreteVelocityScheme &dsv);

    void writeF(double x, string filename);
    void printDomain(Domain d, string filename);
    void printDomain(double x, vector <Domain> domains, string filename);
    void printDomain(vector<Domain> d, string filename);

    void testFuntions();

   private:
    double cells;
    double lower_bound;
    double upper_bound;
    double dx;
    vector<double> x_pos;
    vector<Domain> U;

    vector<double> rho();
    vector<double> p();
    vector<double> q();
    vector<double> u();

    vector<double> F_minus_half(int index);
    vector<double> F_plus_half(int index);
    vector<double> F_flux(int index, double dx, double dt);

};

#endif