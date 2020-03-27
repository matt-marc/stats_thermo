#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "DVS.h"

using namespace std;

int main(void) {
    double rho = 1.225;
    double p = 101.325 * 1000;
    double gamma = 3.0;

    double mach = 2.0;

    double v = mach * sqrt(gamma * (p / rho));

    //creates a dvs that has 500 cells with x range [0m - 20m]
    DiscreteVelocityScheme dvs(300, 0, 20);

    //sets the velocity space to be 100 cells ranging from
    // [-1000m/s to 1000m/s]
    dvs.setVelocitySpace(300, -10000, 10000);

    //creates a density function with rho, u, p
    density_function left(rho, v, p);

    //creates a density function with rho, v, p
    density_function right(rho, 0.0, p);

    //sets the density based on the function to the range
    // of [0 - 5m] and [5m - 10m]
    dvs.setDensityInRange(0, 2, left);
    dvs.setDensityInRange(2, 20, right);

    dvs.write_U("initial_con");

    dvs.time_march_to(0.0006, 1E-6, 1E-7);

    dvs.write_U("ms6");

    return 0;
}