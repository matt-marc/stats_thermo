#include "Domain.h"
using namespace std;

Domain::Domain(double num, double _min, double _max) {
    min = _min;
    max = _max;
    dv = (max - min) / (num - 1);
    vel_space.push_back(min);
    num_dis.push_back(0.0);

    //sets the velocity space to range betwen min and max
    // inital sets the number distribution to be zero
    for (int i = 1; i < num; ++i) {
        vel_space.push_back(vel_space[i - 1] + dv);
        num_dis.push_back(0.0);
    }
}

Domain::Domain(vector<double> dis, vector<double> vels) {
    for (size_t i = 0; i < dis.size(); ++i) {
        num_dis.push_back(dis[i]);
        vel_space.push_back(vels[i]);
    }
}

void Domain::setNumberDensity(density_function eq) {
    for (size_t i = 0; i < vel_space.size(); ++i) {
        num_dis[i] = eq(vel_space[i]) * dv;
    }
}

double Domain::rho() {
    double total = 0.0;

    for (size_t i = 0; i < num_dis.size(); ++i) {
        total += (num_dis[i]);
    }
    return total;
}

double Domain::u() {
    double total = 0;

    for (size_t i = 0; i < num_dis.size(); ++i) {
        total += vel_space[i] * num_dis[i];
    }
    return total;
}

double Domain::p() {
    double total = 0;
    double uu = u();
    double c;

    for (size_t i = 0; i < num_dis.size(); ++i) {
        c = vel_space[i] - uu;
        total += (c * c * num_dis[i]);
    }
    return total;
}

double Domain::q() {
    double total = 0;
    double uu = u();
    double c;

    for (size_t i = 0; i < num_dis.size(); ++i) {
        c = vel_space[i] - uu;
        total += c * c * c * num_dis[i];
    }
    return 0.5 * total;
}

Domain Domain::operator+(const Domain& d) {
    Domain new_d;
    for (size_t i = 0; i < d.vel_space.size(); ++i) {
        new_d.vel_space.push_back(d.vel_space[i] + vel_space[i]);
        new_d.num_dis.push_back(num_dis[i] + d.num_dis[i]);
    }

    return new_d;
}