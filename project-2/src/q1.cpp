#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "DVS.h"

using namespace std;

int main(void) {
    //creates a dvs that has 500 cells with x range [0m - 10m]
    DiscreteVelocityScheme dvs(500, 0, 10);

    //sets the velocity space to be 500 cells ranging from
    // [-1000m/s to 1000m/s]
    dvs.setVelocitySpace(500, -1000, 1000);

    //creates a density function with rho, u, p
    density_function left(4.696, 0.0, 404.0 * 1000.0, dvs.dV());

    //creates a density function with rho, v, p
    density_function right(1.408, 0.0, 101.1 * 1000.0, dvs.dV());

    //sets the density based on the function to the range
    // of [0 - 5m] and [5m - 10m]
    dvs.setDensityInRange(0, 5, left);
    dvs.setDensityInRange(5, 10, right);

    //creates init F distribution graphs
    dvs.write_F("left_init", 0.5);
    dvs.write_F("right_init", 8.4);

    //prints out inital rho, u, p, q of domain
    dvs.write_U("initial_con");

    cout << "Domain set for Q1" << endl;
    cout << dvs;
    cout << "Begin time marching" << endl;

    //begin time marching to tf with dt and tau
    dvs.time_march_to(6.0E-3, 1.0E-5, 10.0);

    dvs.write_U("ms6");

    //prints out F function at given x location
    dvs.write_F("left_final", 3.0);
    dvs.write_F("right_final", 8.0);
    dvs.write_F("shock", 4.5);

    return 0;
}