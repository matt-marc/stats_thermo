#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "DVS.h"

using namespace std;

int main(void) {
    //sets upstream conditions
    double rho_l = 1.225;
    double p_l = 101.325 * 1000;
    double gamma = 3.0;

    //mach pre shock
    double mach = 4.0;

    //calculating downstream conditions right
    //see report for eq_1 and 2
    double eq_1 = (2.0 * gamma * mach * mach - (gamma - 1.0)) / (gamma + 1.0);
    double eq_2 = ((gamma + 1.0) * mach * mach) / ((gamma - 1.0) * mach * mach + 2.0);
    double mach_r = sqrt(((gamma - 1.0) * mach * mach + 2.0) / (2.0 * gamma * mach * mach - (gamma - 1.0)));

    //sets downstream rho and p
    double p_r = p_l * eq_1;
    double rho_r = rho_l * eq_2;

    //sets velocity for upstream and downstram based on computed values
    double v_l = mach * sqrt(gamma * (p_l / rho_l));
    double v_r = mach_r * sqrt(gamma * (p_r / rho_r));

    //creates a dvs that has 500 cells with x range [0m - 0.3m]
    DiscreteVelocityScheme dvs(50, 0, 0.5);

    //sets the velocity space to be 100 cells ranging from
    // [-20000/s to 20000/s]
    dvs.setVelocitySpace(1000, -5000, 5000);

    //creates a density function with rho, u, p
    density_function left(rho_l, v_l, p_l, dvs.dV());

    //creates a density function with rho, v, p
    density_function right(rho_r, v_r, p_r, dvs.dV());

    //sets the density based on the function to the range
    // of [0 - 5m] and [5m - 10m]
    dvs.setDensityInRange(0, 0.3, left);
    dvs.setDensityInRange(0.3, 0.5, right);

    dvs.write_U("initial_con");

    cout << dvs << endl;

    cout << "Begin shock time maching" << endl;
    cout << "UPSTEAM CONDITIONS" << endl;
    cout << "MACH:       " << mach << endl;
    cout << "Rho:        " << rho_l << "kg/m^3" << endl;
    cout << "Pressure:   " << p_l << "Pa" << endl;
    cout << "Velocity:   " << v_l << "m/s" << endl;

    cout << endl;

    cout << "DONWSTREAM CONDITIONS" << endl;
    cout << "Rho:        " << rho_r << "kg/m^3" << endl;
    cout << "Pressure:   " << p_r << "Pa" << endl;
    cout << "Velocity:   " << v_r << "m/s" << endl;
    cout << endl;

    dvs.time_march_to(6.0E-4, 1E-6, 1E-7);

    dvs.write_U("ms6");
    dvs.write_F("shock", 0.355);
    dvs.write_F("left_final", 0.1);

    return 0;
}