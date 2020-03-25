
#include <fstream>
#include <iostream>

#include "DVS.h"

using namespace std;

DiscreteVelocityScheme::DiscreteVelocityScheme(double cells, double lb, double ub) {
    lower_bound = lb;
    upper_bound = ub;

    // sets the x dir cell positions

    dx = (upper_bound - lower_bound) / cells;

    x_pos.push_back(lower_bound + (dx / 2.0));

    for (int i = 1; i < cells; ++i) {
        x_pos.push_back(x_pos[i - 1] + dx);
    }
}

void DiscreteVelocityScheme::setVelocitySpace(double num, double min, double max) {
    dv = (max - min) / (num-1);

    vel_space.push_back(min);

    for (int i = 1; i < num; ++i) {
        vel_space.push_back(vel_space[i-1] + dv);
    }
}

void DiscreteVelocityScheme::setDensityInRange(double min, double max, density_function eq) {
    for (size_t i = 0; i < x_pos.size(); ++i) {
        if (x_pos[i] >= min && x_pos[i] <= max) {
        }
    }
}

vector<double> DiscreteVelocityScheme::rho() {
    vector<double> rho;

    return rho;
}

vector<double> DiscreteVelocityScheme::p() {
    vector<double> p;

    return p;
}

vector<double> DiscreteVelocityScheme::u() {
    vector<double> u;

    return u;
}

vector<double> DiscreteVelocityScheme::q() {
    vector<double> q;

    return q;
}

void DiscreteVelocityScheme::time_march_to(double tf, double dt, double tau) {
    const double one_p_tau = (1.0 + (1.0 / tau));

    double c_t = 0.0;
}

vector<double> DiscreteVelocityScheme::F_minus_half(int index) {
    vector<double> f_i_mf;

    return f_i_mf;
}

vector<double> DiscreteVelocityScheme::F_plus_half(int index) {
    vector<double> f_i_pf;

    return f_i_pf;
}

vector<double> DiscreteVelocityScheme::F_flux(int index, double _dx, double dt) {
    vector<double> f_flux;

    return f_flux;
}

void DiscreteVelocityScheme::testFuntions() {
    auto a = F_minus_half(0);
}

ostream& operator<<(std::ostream& out, const DiscreteVelocityScheme(&dvs)) {
    out << "Current Domain for problem" << endl;
    out << "x [" << dvs.lower_bound << " " << dvs.upper_bound << "]" << endl;
    out << "dx " << dvs.dx << endl;

    out << "x pos cell middle cell valls" << endl;
    out << dvs.x_pos[0] << endl;
    out << dvs.x_pos[1] << endl;
    out << "..." << endl;
    out << dvs.x_pos[dvs.x_pos.size() - 2] << endl;
    out << dvs.x_pos[dvs.x_pos.size() - 1] << endl;

    out << endl;
    out << "Velocity space set for each cell" << endl;
    out << "v [" << dvs.vel_space[0] << " " << dvs.vel_space[dvs.vel_space.size() - 1]<< "]" << endl;
    out << "dv " << dvs.dv << endl;
    out << dvs.vel_space[0] << endl;
    out << dvs.vel_space[1] << endl;
    out << "..." << endl;
    out << dvs.vel_space[dvs.vel_space.size() - 2] << endl;
    out << dvs.vel_space[dvs.vel_space.size() - 1] << endl;

    return out;
}

void DiscreteVelocityScheme::writeF(double x, string filename) {
    int index = 0;
    for (size_t i = 0; i < x_pos.size(); ++i) {
        if (x_pos[i] >= x) {
            index = i;
            break;
        }
    }

    ofstream outfile;
    outfile.open(filename + ".dat");
    outfile << "# F distribution at position " << x_pos[index] << endl;
    outfile << "vel num_dis" << endl;
    outfile.close();

    outfile.open(filename + ".dat", ios_base::app);
}