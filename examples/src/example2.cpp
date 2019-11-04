#include "trjtools/trajectory.hpp"
#include <iostream>
#include <iomanip>

namespace trj = trjtools::trajectory;

int main(int argc, char**argv)
{
    double dbar = 1000.;
    double jbar = 50.;
    double abar = 5.;
    double vbar = 1.5;
    double xbar = 1.;
    auto [tv, ta, tj, td] = trj::axis1d::fourth_order(xbar,vbar,abar,jbar,dbar);
    std::cout << "[" << tv << "," << ta << "," << tj << "," << td << "]" << std::endl;
    auto time_pts = trj::axis1d::fourth_order_time_points(tv,ta,tj,td);
    for (const auto tp : time_pts) 
        std::cout << "[" << tp << "]" << std::endl;
    auto [xp, vp, ap, jp, dp] = trj::axis1d::fourth_order_pivot(xbar, vbar, abar, jbar, dbar, time_pts);
    std::cout << std::fixed << std::setprecision(4);
    for (size_t i=0; i<std::size(xp); i++)
        std::cout << i << "\t: [" << xp[i] << ",\t" << vp[i] << ",\t" << ap[i] << ", \t" << jp[i] << ", \t" << dp[i] << "]" << std::endl;
    return 0;
}