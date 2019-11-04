#include "trjtools/trajectory.hpp"
#include <iostream>

namespace trj = trjtools::trajectory;

int main(int argc, char**argv)
{
    double abar = 5.;
    double vbar = 1.5;
    double xbar = 1.;
    auto [tx, tv, ta] = trj::axis1d::second_order(xbar,vbar,abar);
    auto time_points = std::vector<double>{
        0.2, 0.4, 0.6, 0.8
    };
    auto [xs, vs, as] = trj::axis1d::second_order_profile(
        xbar, vbar, abar, tx, tv, ta, time_points
    );
    assert(xs.size()==vs.size() && vs.size()==as.size());
    for (size_t i=0; i<xs.size(); i++)
        std::cout << "[" << xs[i] << "; " << vs[i] << "; " << as[i] << "]" << std::endl;
    return 0;
}