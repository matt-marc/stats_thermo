#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <iostream>
#include <vector>

using namespace std;

class Domain {
   public:
    vector<double> num_dis;
    vector<double> vel_space;
    Domain(double num, double min, double max);

    template <typename Dis_Eq>
    void setNumberDensity(Dis_Eq eq);

    double rho();
    double u();
    double p();
    double q();

    double min;
    double max;
    double dv;

   
    

};

#endif