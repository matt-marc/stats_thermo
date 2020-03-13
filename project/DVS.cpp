#include "DVS.h"

using namespace std;

DiscreteVelocityScheme::DiscreteVelocityScheme(double _cells, double lb, double ub) {
    cells = _cells;
    lower_bound = lb;
    upper_bound = ub;

    // sets the x dir cell positions

    dx = (upper_bound - lower_bound) / (cells - 1);

    x_pos.push_back(lower_bound + (dx / 2));

    for (int i = 1; i < cells; ++i) {
        x_pos.push_back(x_pos[i - 1] + dx);
    }
}

ostream& operator<<(std::ostream& out, const DiscreteVelocityScheme(&dsv)) {
    out << "Current Domain for problem" << endl;
    out << "x [" << dsv.lower_bound << " " << dsv.upper_bound << "]" << endl;
    for (auto x : dsv.x_pos) {
        out << x;
    }
    out << endl;
    return out;
}