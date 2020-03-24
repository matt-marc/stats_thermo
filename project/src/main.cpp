#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "DVS.h"

using namespace std;

int main(void) {
    //creates a dvs that has 500 cells with x range [0m - 10m]
    DiscreteVelocityScheme dsv(50, 0, 10);


    //sets the velocity space to be 100 cells ranging from 
    // [-1000m/s to 1000m/s]
    dsv.setVelocitySpace(100, -1000, 1000);


    //creates a density function with rho, v, p
    density_function left(4.696, 0.0, 404.0 * 1000.0);

    //creates a density function with rho, v, p
    density_function right(1.408, 0.0, 101.1 * 1000.0);


    //sets the density based on the function to the range
    // of [0 - 5m] and [5m - 10m]
    dsv.setDensityInRange(0, 5, left);
    dsv.setDensityInRange(5, 10, right);

    //dsv.writeF(0.5, "left");
    //dsv.writeF(8.4, "right");

    dsv.printDomain(dsv.U,"test");


    //dsv.time_march_to(0.06,1E-4,1.0);

    //cout << dsv;

    return 0;
}