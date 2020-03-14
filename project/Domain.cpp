#include "Domain.h"
using namespace std;

Domain::Domain(double num, double _min, double _max) {
    min = _min;
    max = _max;
    dv = (max - min) / (num - 1);
    vel_space.push_back(min);

    for (int i = 1; i < num; ++i) {
        vel_space.push_back(vel_space[i - 1] + dv);
    }

    for (int i = 0; i < num; ++i) {
        num_dis.push_back(0.0);
    }
}

void Domain::setNumberDensity(density_function eq) {
    for (size_t i = 0; i < vel_space.size(); ++i) {
        num_dis[i] = eq(vel_space[i]);
    }
}

double Domain::rho() {
    double total = 0;

    for (size_t i = 0; i < num_dis.size(); ++i) {
        total += num_dis[i];
    }
    return total;
}

double Domain::u() {
    double total = 0;
    double r = rho();

    for (size_t i = 0; i < num_dis.size(); ++i) {
        total += vel_space[i] * num_dis[i];
    }
    return total / r;
}

double Domain::p() {
    double total = 0;
    double uu = u();
    double c;

    for (size_t i = 0; i < num_dis.size(); ++i) {
        c = vel_space[i] - uu;
        total += c * c * num_dis[i];
    }
    return 0.5 * total;
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