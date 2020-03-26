
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
    dv = (max - min) / (num - 1);

    vel_space.push_back(min);

    for (int i = 1; i < num; ++i) {
        vel_space.push_back(vel_space[i - 1] + dv);
    }
}

void DiscreteVelocityScheme::init_U() {
    U = vector<vector<double>>(x_pos.size(), vector<double>(vel_space.size(), 0.0));
}

void DiscreteVelocityScheme::setDensityInRange(double min, double max, density_function eq) {
    for (size_t i = 0; i < x_pos.size(); ++i) {
        if (x_pos[i] >= min && x_pos[i] <= max) {
            for (size_t j = 0; j < vel_space.size(); ++j) {
                U[i][j] = eq(vel_space[j]) * dv;
            }
        }
    }
}

vector<double> DiscreteVelocityScheme::rho() {
    vector<double> rho;
    double total;

    for (size_t i = 0; i < x_pos.size(); ++i) {
        total = 0.0;
        for (size_t j = 0; j < vel_space.size(); ++j) {
            total += U[i][j];
        }
        rho.push_back(total);
    }
    return rho;
}

double DiscreteVelocityScheme::rho(int index) {
    double total = 0.0;

    for (size_t j = 0; j < vel_space.size(); ++j) {
        total += U[index][j];
    }
    return total;
}

vector<double> DiscreteVelocityScheme::u() {
    vector<double> u;

    double total;

    for (size_t i = 0; i < x_pos.size(); ++i) {
        total = 0.0;
        for (size_t j = 0; j < vel_space.size(); ++j) {
            total += U[i][j] * vel_space[j];
        }
        u.push_back(total / rho(i));
    }

    return u;
}

double DiscreteVelocityScheme::u(int index) {
    double total = 0.0;

    for (size_t j = 0; j < vel_space.size(); ++j) {
        total += U[index][j] * vel_space[j];
    }
    return total / rho(index);
}

vector<double> DiscreteVelocityScheme::p() {
    vector<double> p;

    double total;
    double c;
    double uu;

    for (size_t i = 0; i < x_pos.size(); ++i) {
        total = 0.0;
        uu = u(i);
        for (size_t j = 0; j < vel_space.size(); ++j) {
            c = vel_space[j] - uu;
            total += U[i][j] * c * c;
        }
        p.push_back(total);
    }
    return p;
}

vector<double> DiscreteVelocityScheme::q() {
    vector<double> q;

    double total;
    double c;
    double uu;

    for (size_t i = 0; i < x_pos.size(); ++i) {
        total = 0.0;
        uu = u(i);
        for (size_t j = 0; j < vel_space.size(); ++j) {
            c = vel_space[j] - uu;
            total += U[i][j] * c * c * c;
        }
        q.push_back(total / 2.0);
    }
    return q;
}

void DiscreteVelocityScheme::time_march_to(double tf, double dt, double tau) {
    const double one_p_tau = (1.0 + (1.0 / tau));
    double c_t = 0.0;

    vector<vector<double>> u_hat;

    while (c_t < tf) {
        if (c_t + dt > tf) {
            dt = tf - c_t;
        }

        for (size_t i = 0; i < x_pos.size(); ++i) {
            auto flux = F_flux(i, dx, dt);
            vector<double> u_i;
            for (size_t j = 0; j < vel_space.size(); ++j) {
                u_i.push_back(U[i][j] + flux[j]);
            }
            u_hat.push_back(u_i);
        }
        //cout << c_t << endl;

        U = u_hat;
        c_t += dt;
    }
}

vector<double> DiscreteVelocityScheme::F_flux(int index, double _dx, double dt) {
    vector<double> f_flux;

    double v_i;
    double f_i_mf;
    double f_i_pf;

    double f_i;
    double f_i_m1;
    double f_i_p1;

    for (size_t i = 0; i < vel_space.size(); ++i) {
        if (index == 0) {
            f_i = U[index][i];
            f_i_m1 = U[index][i];
            f_i_p1 = U[index + 1][i];
        } else if (index == x_pos.size() - 1) {
            f_i = U[index][i];
            f_i_m1 = U[index - 1][i];
            f_i_p1 = U[index][i];
        } else {
            f_i = U[index][i];
            f_i_m1 = U[index - 1][i];
            f_i_p1 = U[index + 1][i];
        }

        v_i = vel_space[i];

        if (v_i >= 0) {
            f_i_mf = v_i * f_i;
            f_i_pf = v_i * f_i;

        } else {
            f_i_mf = v_i * f_i_m1;
            f_i_pf = v_i * f_i_p1;
        }

        double flux = dt * (f_i_mf - f_i_pf) / _dx;
        f_flux.push_back(flux);
    }

    return f_flux;
}

void DiscreteVelocityScheme::testFuntions() {
    string filename = "test";
    ofstream outfile;
    outfile.open(filename + ".dat");
    outfile << "# F distribution at position" << endl;
    outfile.close();

    outfile.open(filename + ".dat", ios_base::app);

    auto r = rho();

    auto test = U[0];
    for (auto t : test) {
        // cout << t << endl;
    }

    for (size_t i = 0; i < x_pos.size(); ++i) {
        outfile << x_pos[i] << " " << r[i] << endl;
    }
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
    out << "v [" << dvs.vel_space[0] << " " << dvs.vel_space[dvs.vel_space.size() - 1] << "]" << endl;
    out << "dv " << dvs.dv << endl;
    out << dvs.vel_space[0] << endl;
    out << dvs.vel_space[1] << endl;
    out << "..." << endl;
    out << dvs.vel_space[dvs.vel_space.size() - 2] << endl;
    out << dvs.vel_space[dvs.vel_space.size() - 1] << endl;

    return out;
}

void DiscreteVelocityScheme::write_U(string filename) {
    ofstream outfile;
    outfile.open(filename + ".dat");
    outfile << "# x rho u p q" << endl;
    outfile.close();

    outfile.open(filename + ".dat", ios_base::app);

    auto r = rho();
    auto uu = u();
    auto pp = p();
    auto qq = q();

    for (size_t i = 0; i < x_pos.size(); ++i) {
        outfile << x_pos[i] << " " << r[i] << " " << uu[i] << " " << pp[i] << " " << qq[i] << endl;
    }
}