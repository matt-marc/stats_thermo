#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>

#include "DVS.h"

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
        double expon = exp(-0.5 * rho / p*(v-u)*(v-u));
        double front = pow((rho/(2*M_PI*p)),0.5);
        return front*expon;
    }
};

int main(void) {
    DiscreteVelocityScheme dsv(5, 0, 1);
    dsv.setVelocitySpace(5,-4,1);

    //dsv.set

    cout << dsv;

    return 0;
}