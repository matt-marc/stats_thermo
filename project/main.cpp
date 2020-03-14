#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "DVS.h"

using namespace std;

int main(void) {
    DiscreteVelocityScheme dsv(5, 0, 10);

    dsv.setVelocitySpace(1000, -30, 30);

    density_function left(4.696, 0.0, 404.0);

    density_function right(1.408, 0.0, 101.1);

    dsv.setDensityInRange(0, 5, left);
    dsv.setDensityInRange(5, 10, right);

    dsv.writeF(0.0, "left");
    dsv.writeF(8.4, "right");

    //dsv.set

    cout << dsv;

    return 0;
}