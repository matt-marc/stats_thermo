
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

void DiscreteVelocityScheme::time_march_to(double tf, double dt, double tau) {
    const double one_p_tau = (1.0 + (1.0 / tau));

    vector<Domain> uf;

    double c_t = 0.0;
    while (c_t < tf) {
        vector<Domain> u_hat;
        vector<double> f_flux;
        vector<Domain> f_flux_domain;
        if (c_t + dt > tf) {
            dt = tf - c_t;
        }
        // probably not the most efficent thing in the world
        for (size_t i = 0; i < U.size(); ++i) {
            f_flux = F_flux(i, dx, dt);
            Domain d(f_flux, U[0].vel_space);
            u_hat.push_back(U[i]+d); 
        }

        c_t += dt;
        uf = u_hat;
    }
    //printDomain(5.0, uf, "test");
    printDomain(uf, "ms6");
}

vector<double> DiscreteVelocityScheme::F_minus_half(int index) {
    vector<double> f_i_mf;
    vector<double> f_i_m1;

    vector<double> f_i = U[index].num_dis;
    vector<double> v_i = U[index].vel_space;

    //edge case, if we are at the left most cell
    // then F_(i-1) woudl be zero
    if (index == 0) {
        f_i_m1.resize(v_i.size(), 0.0);
    } else {
        f_i_m1 = U[index - 1].num_dis;
    }

    for (size_t i = 0; i < f_i.size(); ++i) {
        //if moving to left use F_i
        if (v_i[i] < 0) {
            f_i_mf.push_back(v_i[i] * f_i[i]);
        } else {
            f_i_mf.push_back(v_i[i] * f_i_m1[i]);
        }
    }

    return f_i_mf;
}

vector<double> DiscreteVelocityScheme::F_plus_half(int index) {
    vector<double> f_i_pf;
    vector<double> f_i_p1;

    vector<double> f_i = U[index].num_dis;
    vector<double> v_i = U[index].vel_space;

    //egde case if we are at the rightmost cell
    // then F_(i+1) = 0
    if (index == (int)x_pos.size() - 1) {
        f_i_p1.resize(v_i.size(), 0.0);
    } else {
        f_i_p1 = U[index + 1].num_dis;
    }

    for (size_t i = 0; i < f_i.size(); ++i) {
        //if moving to right use F_i
        if (v_i[i] > 0) {
            f_i_pf.push_back(v_i[i] * f_i[i]);
        } else {
            f_i_pf.push_back(v_i[i] * f_i_p1[i]);
        }
    }

    return f_i_pf;
}

vector<double> DiscreteVelocityScheme::F_flux(int index, double _dx, double dt) {
    vector<double> f_flux;
    vector<double> f_mf = F_minus_half(index);
    vector<double> f_pf = F_plus_half(index);

    for (size_t i = 0; i < f_mf.size(); ++i) {
        f_flux.push_back(dt * (f_mf[i] - f_pf[i]) / _dx);
    }

    return f_flux;
}

void DiscreteVelocityScheme::testFuntions() {
    auto a = F_minus_half(0);
}

void DiscreteVelocityScheme::printDomain(double x, vector<Domain> d, string filename) {
    int index = 0;
    for (size_t i = 0; i < x_pos.size(); ++i) {
        if (x_pos[i] >= x) {
            index = i;
            break;
        }
    }
    printDomain(d[index], filename);
}

void DiscreteVelocityScheme::printDomain(Domain d, string filename) {
    filename += ".dat";
    ofstream outfile;
    outfile.open(filename);
    outfile << "# F distribution" << endl;
    outfile << "vel num_dis" << endl;
    outfile.close();

    outfile.open(filename, ios_base::app);

    for (size_t i = 0; i < d.vel_space.size(); ++i) {
        outfile << d.vel_space[i] << " " << d.num_dis[i] << endl;
    }
}

void DiscreteVelocityScheme::printDomain(vector<Domain> d, string filename) {
    filename += ".dat";
    ofstream outfile;
    outfile.open(filename);
    outfile << "# domain properties" << endl;
    outfile << "x rho u p q" << endl;
    outfile.close();

    outfile.open(filename, ios_base::app);

    for (size_t i = 0; i < x_pos.size(); ++i) {
        outfile << x_pos[i] << " " << d[i].rho() << " " << d[i].u() << endl;
    }
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
    outfile.open(filename + ".dat");
    outfile << "# F distribution at position " << x_pos[index] << endl;
    outfile << "vel num_dis" << x_pos[index] << endl;
    outfile.close();

    outfile.open(filename + ".dat", ios_base::app);

    Domain d = U[index];

    for (size_t i = 0; i < d.vel_space.size(); ++i) {
        outfile << d.vel_space[i] << " " << d.num_dis[i] << endl;
    }
}