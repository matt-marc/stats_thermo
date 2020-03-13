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
}

template <typename Dis_Eq>
void Domain::setNumberDensity(Dis_Eq eq) {
    for (auto v : vel_space) {
        num_dis = eq(v);
    }
}

double Domain::rho() {
    double total = 0;

    for (int i = 0; i < num_dis.size(); ++i) {
        total += num_dis[i];
    }
    return total;
}

double Domain::u() {
    double total = 0;
    double r = rho();

    for (int i = 0; i < num_dis.size(); ++i) {
        total += vel_space[i] * num_dis[i];
    }
    return total / r;
}

double Domain::p() {
    double total = 0;
    double uu = u();
    double c;

    for (int i = 0; i < num_dis.size(); ++i) {
        c = vel_space[i] - uu;
        total += c * c * num_dis[i];
    }
    return 0.5 * total;
}

double Domain::q() {
    double total = 0;
    double uu = u();
    double c;

    for (int i = 0; i < num_dis.size(); ++i) {
        c = vel_space[i] - uu;
        total += c * c * c * num_dis[i];
    }
    return 0.5 * total;
}