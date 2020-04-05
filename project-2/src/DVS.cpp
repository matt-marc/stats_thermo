
#include <fstream>
#include <iostream>

#include "DVS.h"

using namespace std;

//Constructor that sets the x cells
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

//sets the velocity space vector
void DiscreteVelocityScheme::setVelocitySpace(double num, double min, double max) {
    dv = (max - min) / (num - 1);

    vel_space.push_back(min);

    for (int i = 1; i < num; ++i) {
        vel_space.push_back(vel_space[i - 1] + dv);
    }
    //sets the size of the U matrix, filled with zeros
    init_U();
}

//createds matrix U of size [x_pos][vel_space]
void DiscreteVelocityScheme::init_U() {
    U = vector<vector<double>>(x_pos.size(), vector<double>(vel_space.size(), 0.0));
}

//sets density function of given x-pos range [min - max]
void DiscreteVelocityScheme::setDensityInRange(double min, double max, density_function eq) {
    for (size_t i = 0; i < x_pos.size(); ++i) {
        if (x_pos[i] >= min && x_pos[i] <= max) {
            for (size_t j = 0; j < vel_space.size(); ++j) {
                U[i][j] = eq(vel_space[j]);
            }
        }
    }
}

//Computes rho for U matrix
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

//computes rho of matrix U at given index
double DiscreteVelocityScheme::rho(int index) {
    double total = 0.0;

    for (size_t j = 0; j < vel_space.size(); ++j) {
        total += U[index][j];
    }
    return total;
}

//computes rho for given distribution u
double DiscreteVelocityScheme::rho(vector<double> u) {
    double total = 0.0;
    for (size_t j = 0; j < vel_space.size(); ++j) {
        total += u[j];
    }
    return total;
}

//Computes u for U matrix
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

//computes u for matrix U at index
double DiscreteVelocityScheme::u(int index) {
    double total = 0.0;

    for (size_t j = 0; j < vel_space.size(); ++j) {
        total += U[index][j] * vel_space[j];
    }
    return total / rho(index);
}

//computes u for given distribution u
double DiscreteVelocityScheme::u(vector<double> u) {
    double total = 0.0;

    for (size_t j = 0; j < vel_space.size(); ++j) {
        total += u[j] * vel_space[j];
    }
    return total / rho(u);
}

//computes pressure for U matrix
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

//computes pressure for given distribution f
double DiscreteVelocityScheme::p(vector<double> f) {
    vector<double> p;

    double total;
    double c;
    double uu;

    total = 0.0;
    uu = u(f);
    for (size_t j = 0; j < vel_space.size(); ++j) {
        c = vel_space[j] - uu;
        total += f[j] * c * c;
    }

    return total;
}

//computes q for U matrix
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
    double c_t = 0.0;

    while (c_t < tf) {
        if (c_t + dt > tf) {
            dt = tf - c_t;
        }

        double one_p_tau = (1.0 + (dt / tau));
        vector<vector<double>> u_n_p1;
        vector<vector<double>> u_hat;

        for (size_t i = 0; i < x_pos.size(); ++i) {
            vector<double> u_i_hat;
            for (size_t j = 0; j < vel_space.size(); ++j) {
                u_i_hat.push_back(U[i][j] + F_flux(i, j, dx, dt));
            }
            u_hat.push_back(u_i_hat);

            double m_rho = rho(u_hat[i]);
            double m_u = u(u_hat[i]);
            double m_p = p(u_hat[i]);

            vector<double> m_i;
            vector<double> u_n_i;

            density_function m(m_rho, m_u, m_p, dv);

            for (size_t j = 0; j < vel_space.size(); ++j) {
                m_i.push_back(dt * m(vel_space[j]) / tau);

                u_n_i.push_back((u_hat[i][j] + m_i[j]) / one_p_tau);
            }
            u_n_p1.push_back(u_n_i);
        }

        U = u_n_p1;
        c_t += dt;
    }
}

//Function computes the flux at given index x posiition and point in vel_space
double DiscreteVelocityScheme::F_flux(int index, int vel_index, double _dx, double dt) {
    double f_i_mf;
    double f_i_pf;

    double f_i;
    double f_i_m1;
    double f_i_p1;

    double v_i = vel_space[vel_index];

    if (index == 0) {
        f_i = U[index][vel_index];
        f_i_m1 = U[index][vel_index];
        f_i_p1 = U[index + 1][vel_index];
    } else if (index == (int)x_pos.size() - 1) {
        f_i = U[index][vel_index];
        f_i_m1 = U[index - 1][vel_index];
        f_i_p1 = U[index][vel_index];
    } else {
        f_i = U[index][vel_index];
        f_i_m1 = U[index - 1][vel_index];
        f_i_p1 = U[index + 1][vel_index];
    }

    if (v_i < 0) {
        f_i_mf = v_i * f_i;
        f_i_pf = v_i * f_i_p1;

    } else {
        f_i_mf = v_i * f_i_m1;
        f_i_pf = v_i * f_i;
    }

    double flux = dt * (f_i_mf - f_i_pf) / _dx;

    return flux;
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

void DiscreteVelocityScheme::write_F(string filename, double x) {
    int index = 0;

    for (size_t i = 0; i < x_pos.size(); ++i) {
        if (x_pos[i] >= x) {
            index = i;
            break;
        }
    }

    double r = rho(index);
    double uu = u(index);

    ofstream outfile;
    outfile.open(filename + ".dat");
    outfile << "# At position " << x_pos[index] << endl;
    outfile << "# values are following" << endl;
    outfile << "# rho " << r << endl;
    outfile << "# u " << uu << endl;

    outfile.close();

    outfile.open(filename + ".dat", ios_base::app);

    for (size_t i = 0; i < vel_space.size(); ++i) {
        outfile << vel_space[i] << " " << U[index][i] << endl;
    }
}