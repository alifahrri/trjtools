#include "trjtools/trajectory.hpp"
#include "trjtools/utility.hpp"
#include <iostream>
#include <iomanip>
#include <thread>
#include <chrono>

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
    auto profile_time_pts = std::vector<double>{
        0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1
    };
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        auto [xfn, vfn, afn, jfn, dfn] = trj::axis1d::helper::fourth_order_profile_fn(dbar, time_pts);
        auto [xp, vp, ap, jp, dp] = trj::axis1d::fourth_order_profile<50>(xbar, vbar, abar, jbar, dbar, tv, ta, tj, td, profile_time_pts);
        auto t1 = std::chrono::high_resolution_clock::now();
        std::cout << "completed in " << (t1-t0).count() << "ns"  << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        for (size_t i=0; i<std::size(xp); i++){
            std::cout << "d " << dfn(profile_time_pts[i]) << std::endl;
            std::cout << i << "\t: [" << xp[i] << ",\t" << vp[i] << ",\t" << ap[i] << ", \t" << jp[i] << ", \t" << dp[i] << "]" << std::endl;
        }
    }
    {
        constexpr size_t time_segment = 25;
        constexpr size_t integral_segment = 25;
        auto t0 = std::chrono::high_resolution_clock::now();
        auto [xfn, vfn, afn, jfn, dfn] = trj::axis1d::helper::fourth_order_profile_fn<time_segment,integral_segment>(dbar, time_pts);
        std::cout << std::fixed << std::setprecision(4);
        std::vector<double> ds, js, as, vs, xs;
        std::thread dthread = std::thread([&](){
            for (const auto t : profile_time_pts)
                ds.push_back(dfn(t));
        });
        std::thread jthread = std::thread([&](){ 
            for (const auto t : profile_time_pts)
                js.push_back(jfn(t));
        });
        std::thread athread = std::thread([&](){
            for (const auto t : profile_time_pts)
                as.push_back(afn(t));
        });
        std::thread vthread = std::thread([&](){
            for (const auto t : profile_time_pts)
                vs.push_back(vfn(t));
        });
        std::thread xthread = std::thread([&](){
            for (const auto t : profile_time_pts)
                xs.push_back(xfn(t));
        });
        dthread.join();
        jthread.join();
        athread.join();
        vthread.join();
        xthread.join();
        auto t1 = std::chrono::high_resolution_clock::now();
        std::cout << "completed in " << (t1-t0).count() << "ns"  << std::endl;
        for (size_t i=0; i<std::size(profile_time_pts); i++){
            std::cout << "time : " << profile_time_pts[i] << std::endl; 
            std::cout << "d " << ds[i] << std::endl;
            std::cout << "j " << js[i] << std::endl;
            std::cout << "a " << as[i] << std::endl;
            std::cout << "v " << vs[i] << std::endl;
            std::cout << "x " << xs[i] << std::endl;
        }
    }
    return 0;
}