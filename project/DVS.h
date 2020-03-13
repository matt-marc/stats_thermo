#ifndef DVS_H
#define DVS_H

#include <iostream>
#include <vector>
#include "Domain.h"

class DiscreteVelocityScheme {
   public:
    DiscreteVelocityScheme(double cells, double lower_bound, double upper_bound);

    template <typename Dis_Eq>
    void setDensityInRange(double min, double max, Dis_Eq eq);

    void setVelocitySpace(double number, double min, double max);

    friend std::ostream &operator<<(std::ostream &out, const DiscreteVelocityScheme &dsv);

   private:
    double cells;
    double lower_bound;
    double upper_bound;
    double dx;
    vector<double> x_pos;
    vector<Domain> dis;

    vector<double> rho();
    vector<double> p();
    vector<double> q();
    vector<double> u();
};

#endif