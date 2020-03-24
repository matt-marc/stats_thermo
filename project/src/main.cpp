#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "DVS.h"

using namespace std;

int main(void) {
    DiscreteVelocityScheme dsv(500, 0, 10);

    dsv.setVelocitySpace(100, -1000, 1000);

    density_function left(4.696, 0.0, 404.0 * 1000.0);

    density_function right(1.408, 0.0, 101.1 * 1000.0);

    dsv.setDensityInRange(0, 5, left);
    dsv.setDensityInRange(5, 10, right);

    dsv.writeF(0.0, "left");
    dsv.writeF(8.4, "right");

    //dsv.set

    //dsv.testFuntions();

    dsv.time_march_to(0.06,1E-4,1.0);

    //cout << dsv;

    return 0;
}