#ifndef DVS_H
#define DVS_H

#include <vector>
#include <iostream>

using namespace std;

class distribution {
   public:
    vector<double> num_dis;

   private:
    double num;
};

class DiscreteVelocityScheme {
   public:
    DiscreteVelocityScheme(double cells, double lower_bound, double upper_bound);

    friend std::ostream& operator<< (std::ostream &out, const DiscreteVelocityScheme &dsv);


   private:
    double cells;
    double lower_bound;
    double upper_bound;
    double dx;
    vector<double> x_pos;
    vector<distribution> dis;
};

#endif