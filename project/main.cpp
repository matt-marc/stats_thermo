#include <fstream>
#include <iostream>
#include <vector>

#include "DVS.h"


using namespace std;

int main(void) {
    DiscreteVelocityScheme dsv(20,0,1);

    cout <<dsv;

    return 0;
}
