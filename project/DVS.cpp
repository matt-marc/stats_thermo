
#include <fstream>
#include <iostream>

#include "DVS.h"
#include "Domain.h"

using namespace std;

DiscreteVelocityScheme::DiscreteVelocityScheme(double _cells, double lb, double ub) {
    cells = _cells;
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
    for (size_t i = 0; i < x_pos.size(); ++i) {
        U.push_back(Domain(num, min, max));
    }
}

void DiscreteVelocityScheme::setDensityInRange(double min, double max, density_function eq) {
    for (size_t i = 0; i < x_pos.size(); ++i) {
        if (x_pos[i] >= min && x_pos[i] <= max) {
            U[i].setNumberDensity(eq);
        }
    }
}

vector<double> DiscreteVelocityScheme::rho() {
    vector<double> rho;
    for (size_t i = 0; i < U.size(); ++i) {
        rho.push_back(U[i].rho());
    }

    return rho;
}

vector<double> DiscreteVelocityScheme::p() {
    vector<double> p;
    for (size_t i = 0; i < U.size(); ++i) {
        p.push_back(U[i].p());
    }

    return p;
}

vector<double> DiscreteVelocityScheme::u() {
    vector<double> u;
    for (size_t i = 0; i < U.size(); ++i) {
        u.push_back(U[i].u());
    }

    return u;
}

vector<double> DiscreteVelocityScheme::q() {
    vector<double> q;
    for (size_t i = 0; i < U.size(); ++i) {
        q.push_back(U[i].q());
    }
    return q;
}

void DiscreteVelocityScheme::time_march_to(double tf, double dt, double tau){

}




ostream& operator<<(std::ostream& out, const DiscreteVelocityScheme(&dsv)) {
    out << "Current Domain for problem" << endl;
    out << "x [" << dsv.lower_bound << " " << dsv.upper_bound << "]" << endl;
    out << "dx " << dsv.dx << endl;
    for (auto x : dsv.x_pos) {
        out << x << endl;
    }
    out << endl;
    out << "Velocity space set for each cell" << endl;
    out << "v [" << dsv.U[0].min << " " << dsv.U[0].max << "]" << endl;
    out << "dx " << dsv.U[0].dv << endl;

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
    outfile.open(filename+".dat");
    outfile << "# F distribution at position " << x_pos[index] << endl;
    outfile << "vel num_dis " << x_pos[index] << endl;
    outfile.close();

    outfile.open(filename+".dat", ios_base::app);

    Domain d = U[index];

    for (size_t i = 0; i < d.vel_space.size(); ++i) {
        outfile << d.vel_space[i] << " " << d.num_dis[i] << endl;
    }
}