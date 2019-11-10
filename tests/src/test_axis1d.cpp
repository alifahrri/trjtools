#include "trjtools/trajectory.hpp"
#include <gtest/gtest.h>
#include <vector>

namespace trj = trjtools::trajectory;

TEST(axis1d, second_order)
{
    double abar = 5.;
    double vbar = 1.5;
    double xbar = 1.;
    auto [tx, tv, ta] = trj::axis1d::second_order(xbar,vbar,abar);
    ASSERT_NEAR(tx,0.97,1e-2);
    ASSERT_NEAR(tv,0.36,1e-2);
    ASSERT_NEAR(ta,0.3,1e-2);
}

TEST(axis1d, constexpr_second_order)
{
    constexpr double abar = 5.;
    constexpr double vbar = 1.5;
    constexpr double xbar = 1.;
    constexpr auto time = trj::axis1d::second_order(xbar,vbar,abar);
    constexpr auto tx = std::get<0>(time);
    constexpr auto tv = std::get<1>(time);
    constexpr auto ta = std::get<2>(time);
    static_assert(
        fabs(tx-0.97) < 1e-2
    );
    static_assert(
        fabs(tv-0.36) < 1e-2
    );
    static_assert(
        fabs(ta-0.3) < 1e-2
    );
}

TEST(axis1d, solver_4th_3rd_poly)
{
    double t_dbar = 0.04;
    double x = 1.;
    double d = 1000.;
    auto [root1, root2, root3] = trj::axis1d::helper::solver_4th_3rd_poly(t_dbar, x, d);
    static_assert(std::is_same_v<std::decay_t<decltype(root1)>,double>);
    static_assert(std::is_same_v<std::decay_t<decltype(root2)>,std::complex<double>>);
    static_assert(std::is_same_v<std::decay_t<decltype(root3)>,std::complex<double>>);
    EXPECT_NEAR(root1, 0.166208032359266, 1e-6);
    EXPECT_NEAR(root2.real(), -0.183104016179633, 1e-6);
    EXPECT_NEAR(root2.imag(), 0.200348785250823, 1e-6);
    EXPECT_NEAR(root3.real(), -0.183104016179633, 1e-6);
    EXPECT_NEAR(root3.imag(), -0.200348785250823, 1e-6);
}

TEST(axis1d, constexpr_solver_4th_3rd_poly)
{
    constexpr double t_dbar = 0.04;
    constexpr double x = 1.;
    constexpr double d = 1000.;
    constexpr auto roots = trj::axis1d::helper::solver_4th_3rd_poly(t_dbar, x, d);
    constexpr auto root1 = std::get<0>(roots);
    constexpr auto root2 = std::get<1>(roots);
    constexpr auto root3 = std::get<2>(roots);
    static_assert(
        fabs(root1-0.166208032359266) < 1e-6
    );
    static_assert(
        fabs(root2.real()-(-0.183104016179633)) < 1e-6
    );
    static_assert(
        fabs(root2.imag()-0.200348785250823) < 1e-6
    );
    static_assert(
        fabs(root3.real()-(-0.183104016179633)) < 1e-6
    );
    static_assert(
        fabs(root3.imag()-(-0.200348785250823)) < 1e-6
    );
}

TEST(axis1d, solver_4th_2nd_poly)
{
    double t_dbar = 0.04;
    double v = .25;
    double d = 1000.;
    auto [root1, root2] = trj::axis1d::helper::solver_4th_2nd_poly(t_dbar, v, d);
    static_assert(std::is_same_v<std::decay_t<decltype(root1)>,double>);
    static_assert(std::is_same_v<std::decay_t<decltype(root2)>,double>);
    EXPECT_NEAR(root1, 0.0215475321515004, 1e-6);
    EXPECT_NEAR(root2, -0.141547532151500, 1e-6);
}

TEST(axis1d, constexpr_solver_4th_2nd_poly)
{
    constexpr double t_dbar = 0.04;
    constexpr double v = .25;
    constexpr double d = 1000.;
    constexpr auto roots = trj::axis1d::helper::solver_4th_2nd_poly(t_dbar, v, d);
    constexpr auto root1 = std::get<0>(roots);
    constexpr auto root2 = std::get<1>(roots);
    static_assert( fabs(root1 - (0.0215475321515004)) < 1e-6);
    static_assert( fabs(root2 - (-0.141547532151500)) < 1e-6);
}

TEST(axis1d, solver_4th_2nd_poly_ta)
{
    double t_dbar = 0.04;
    double t_jbar = 0.08;
    double x = 1.;
    double d = 1000.;
    auto [root1, root2] = trj::axis1d::helper::solver_4th_2nd_poly_ta(t_jbar, t_dbar, x, d);
    static_assert(std::is_same_v<std::decay_t<decltype(root1)>,double>);
    static_assert(std::is_same_v<std::decay_t<decltype(root2)>,double>);
    EXPECT_NEAR(root1, 0.267473196040797, 1e-6);
    EXPECT_NEAR(root2, -0.587473196040797, 1e-6);
}

TEST(axis1d, constexpr_solver_4th_2nd_poly_ta)
{
    constexpr double t_dbar = 0.04;
    constexpr double t_jbar = 0.08;
    constexpr double x = 1.;
    constexpr double d = 1000.;
    constexpr auto roots = trj::axis1d::helper::solver_4th_2nd_poly_ta(t_jbar, t_dbar, x, d);
    constexpr auto root1 = std::get<0>(roots);
    constexpr auto root2 = std::get<1>(roots);
    static_assert( fabs(root1- (0.267473196040797)  ) < 1e-6);
    static_assert( fabs(root2- (-0.587473196040797) ) < 1e-6);
}

TEST(axis1d, fourth_order)
{
    double dbar = 1000.;
    double jbar = 50.;
    double abar = 5.;
    double vbar = 1.5;
    double xbar = 1.;
    auto [tv, ta, tj, td] = trj::axis1d::fourth_order(xbar,vbar,abar,jbar,dbar);
    auto time_pts = trj::axis1d::fourth_order_time_points(tv,ta,tj,td);
    EXPECT_NEAR(time_pts.back(), 1.108, 1e-3);
}

TEST(axis1d, fourth_order_djerk_profile)
{
    double dbar = 1000.;
    double jbar = 50.;
    double abar = 5.;
    double vbar = 1.5;
    double xbar = 1.;
    auto [tv, ta, tj, td] = trj::axis1d::fourth_order(xbar,vbar,abar,jbar,dbar);
    auto time_pts = trj::axis1d::fourth_order_time_points(tv,ta,tj,td);
    auto dj_fn = trj::axis1d::helper::fourth_order_djerk_profile(dbar, time_pts);
    EXPECT_NEAR( dj_fn(0.01), dbar, 1e-6);
    EXPECT_NEAR( dj_fn(0.09), 0., 1e-6);
    EXPECT_NEAR( dj_fn(0.11), -dbar, 1e-6);
    EXPECT_NEAR( dj_fn(0.2), 0., 1e-6);
    EXPECT_NEAR( dj_fn(0.33), -dbar, 1e-6);
    EXPECT_NEAR( dj_fn(0.39), 0., 1e-6);
    EXPECT_NEAR( dj_fn(0.4), dbar, 1e-6);
    EXPECT_NEAR( dj_fn(0.6), 0., 1e-6);
    EXPECT_NEAR( dj_fn(0.7), -dbar, 1e-6);
    EXPECT_NEAR( dj_fn(0.75), 0, 1e-6);
    EXPECT_NEAR( dj_fn(0.8), dbar, 1e-6);
    EXPECT_NEAR( dj_fn(0.9), 0, 1e-6);
    EXPECT_NEAR( dj_fn(1.0), dbar, 1e-6);
    EXPECT_NEAR( dj_fn(1.02), 0, 1e-6);
    EXPECT_NEAR( dj_fn(1.08), -dbar, 1e-6);
    EXPECT_NEAR( dj_fn(1.2), 0, 1e-6);
}

TEST(axis1d, constexpr_fourth_order_djerk_profile)
{
    constexpr double dbar = 1000.;
    constexpr double jbar = 50.;
    constexpr double abar = 5.;
    constexpr double vbar = 1.5;
    constexpr double xbar = 1.;
    constexpr auto time_profile = trj::axis1d::fourth_order(xbar,vbar,abar,jbar,dbar);
    constexpr auto tv = std::get<0>(time_profile);
    constexpr auto ta = std::get<1>(time_profile);
    constexpr auto tj = std::get<2>(time_profile);
    constexpr auto td = std::get<3>(time_profile);
    constexpr auto time_pts = trj::axis1d::fourth_order_time_points(tv,ta,tj,td);
    constexpr auto dj_fn = trj::axis1d::helper::fourth_order_djerk_profile(dbar, time_pts);
    static_assert( fabs(dj_fn(0.01) - (dbar))   < 1e-6);
    static_assert( fabs(dj_fn(0.09) - (0))      < 1e-6);
    static_assert( fabs(dj_fn(0.11) - (-dbar))  < 1e-6);
    static_assert( fabs(dj_fn(0.2)  - (0.))     < 1e-6);
    static_assert( fabs(dj_fn(0.33) - (-dbar))  < 1e-6);
    static_assert( fabs(dj_fn(0.39) - (0))      < 1e-6);
    static_assert( fabs(dj_fn(0.4)  - (dbar))   < 1e-6);
    static_assert( fabs(dj_fn(0.6)  - (0))      < 1e-6);
    static_assert( fabs(dj_fn(0.7)  - (-dbar))  < 1e-6);
    static_assert( fabs(dj_fn(0.75) - (0))      < 1e-6);
    static_assert( fabs(dj_fn(0.8)  - (dbar))   < 1e-6);
    static_assert( fabs(dj_fn(0.9)  - (0))      < 1e-6);
    static_assert( fabs(dj_fn(1.0)  - (dbar))   < 1e-6);
    static_assert( fabs(dj_fn(1.02) - (0))      < 1e-6);
    static_assert( fabs(dj_fn(1.08) - (-dbar))  < 1e-6);
    static_assert( fabs(dj_fn(1.2)  - (0))      < 1e-6);
}

TEST(axis1d, constexpr_fourth_order_profile_fn)
{
    /* 
        note computing on compile time with high precision 
        using gcc will blow-up your memory 
    */
    constexpr double dbar = 1000.;
    constexpr double jbar = 50.;
    constexpr double abar = 5.;
    constexpr double vbar = 1.5;
    constexpr double xbar = 1.;
    constexpr auto time_profile = trj::axis1d::fourth_order(xbar,vbar,abar,jbar,dbar);
    constexpr auto tv = std::get<0>(time_profile);
    constexpr auto ta = std::get<1>(time_profile);
    constexpr auto tj = std::get<2>(time_profile);
    constexpr auto td = std::get<3>(time_profile);
    constexpr auto time_pts = trj::axis1d::fourth_order_time_points(tv,ta,tj,td);
    constexpr auto trajectory = trj::axis1d::helper::fourth_order_profile_fn<25,25>(dbar, time_pts);
    constexpr auto x_fn = std::get<0>(trajectory);
    constexpr auto v_fn = std::get<1>(trajectory);
    constexpr auto a_fn = std::get<2>(trajectory);
    constexpr auto j_fn = std::get<3>(trajectory);
    constexpr auto d_fn = std::get<4>(trajectory);
    static_assert( fabs(d_fn(0.01) - (dbar))   < 1e-6);
    static_assert( fabs(d_fn(0.09) - (0))      < 1e-6);
    static_assert( fabs(d_fn(0.11) - (-dbar))  < 1e-6);
    static_assert( fabs(d_fn(0.2)  - (0.))     < 1e-6);
    static_assert( fabs(d_fn(0.33) - (-dbar))  < 1e-6);
    static_assert( fabs(d_fn(0.39) - (0))      < 1e-6);
    static_assert( fabs(d_fn(0.4)  - (dbar))   < 1e-6);
    static_assert( fabs(d_fn(0.6)  - (0))      < 1e-6);
    static_assert( fabs(d_fn(0.7)  - (-dbar))  < 1e-6);
    static_assert( fabs(d_fn(0.75) - (0))      < 1e-6);
    static_assert( fabs(d_fn(0.8)  - (dbar))   < 1e-6);
    static_assert( fabs(d_fn(0.9)  - (0))      < 1e-6);
    static_assert( fabs(d_fn(1.0)  - (dbar))   < 1e-6);
    static_assert( fabs(d_fn(1.02) - (0))      < 1e-6);
    static_assert( fabs(d_fn(1.08) - (-dbar))  < 1e-6);
    static_assert( fabs(d_fn(1.2)  - (0))      < 1e-6);

    static_assert( fabs(j_fn(0.05) - (49))   < 1e-0);
    static_assert( fabs(j_fn(0.2)  - (22))   < 1e-0);
    static_assert( fabs(j_fn(0.35) - (-21))  < 1e-0);
    static_assert( fabs(j_fn(0.6)  - (21))   < 1e-0);
    static_assert( fabs(j_fn(0.75) - (-19))  < 1e-0);
    static_assert( fabs(j_fn(0.9)  - (26))   < 1e-0);
    static_assert( fabs(j_fn(1.05) - (55))   < 1e-0);

    static_assert( fabs(a_fn(0.2)  - (7.86))   < 1e-0);
    static_assert( fabs(a_fn(0.6)  - (10.9))   < 1e-0);
    static_assert( fabs(a_fn(0.9)  - (11.1))   < 1e-0);

    static_assert( fabs(v_fn(0.6)  - (4.4))   < 1e-0);

    static_assert( fabs(x_fn(1.1)  - (4.7)) < 1e-0);
}

TEST(axis1d, fourth_order_profile_fn)
{
    double dbar = 1000.;
    double jbar = 50.;
    double abar = 5.;
    double vbar = 1.5;
    double xbar = 1.;
    auto [tv, ta, tj, td] = trj::axis1d::fourth_order(xbar,vbar,abar,jbar,dbar);
    auto time_pts = trj::axis1d::fourth_order_time_points(tv,ta,tj,td);
    // auto [xfn, vfn, afn, jfn, dfn] = trj::axis1d::helper::fourth_order_profile_fn<75,50>(dbar, time_pts);
    // EXPECT_NEAR( jfn(0.06), jbar, 1e-1);
    // EXPECT_NEAR( jfn(0.2), 0., 1e-1);
}

TEST(axis1d, constexpr_fourth_order)
{
    constexpr double dbar = 1000.;
    constexpr double jbar = 50.;
    constexpr double abar = 5.;
    constexpr double vbar = 1.5;
    constexpr double xbar = 1.;
    constexpr auto time = trj::axis1d::fourth_order(xbar,vbar,abar,jbar,dbar);
    constexpr auto time_pts = trj::axis1d::fourth_order_time_points(
        std::get<0>(time), std::get<1>(time),
        std::get<2>(time), std::get<3>(time)
    );
    static_assert(
        fabs(time_pts.back() - 1.108) < 1e-3
    );
}

TEST(axis1d, second_order_profile)
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
    ASSERT_EQ(xs.size(), 4);
    ASSERT_EQ(vs.size(), 4);
    ASSERT_EQ(as.size(), 4);
    ASSERT_NEAR(vs.at(0), 1., 1e-3);
    ASSERT_NEAR(as.at(0), 5., 1e-3);
    ASSERT_NEAR(vs.at(1), 1.5, 1e-3);
    ASSERT_NEAR(as.at(1), 0., 1e-3);
    ASSERT_NEAR(vs.at(2), 1.5, 1e-3);
    ASSERT_NEAR(as.at(2), 0., 1e-3);
    ASSERT_NEAR(vs.at(3), 0.83, 1e-2);
    ASSERT_NEAR(as.at(3), -5., 1e-3);
}

TEST(axis1d, constexpr_second_order_profile)
{
    constexpr double abar = 5.;
    constexpr double vbar = 1.5;
    constexpr double xbar = 1.;
    constexpr auto time_profile = trj::axis1d::second_order(xbar,vbar,abar);
    constexpr auto time_points = std::array<double,4>{
        0.2, 0.4, 0.6, 0.8
    };
    constexpr auto tx = std::get<0>(time_profile);
    constexpr auto tv = std::get<1>(time_profile);
    constexpr auto ta = std::get<2>(time_profile);
    constexpr auto trajectory = trj::axis1d::second_order_profile(
        xbar, vbar, abar, tx, tv, ta, time_points
    );
    constexpr auto xs = std::get<0>(trajectory);
    constexpr auto vs = std::get<1>(trajectory);
    constexpr auto as = std::get<2>(trajectory);
    static_assert(xs.size()==4);
    static_assert(vs.size()==4);
    static_assert(as.size()==4);
    static_assert( fabs(vs.at(0)-1.) < 1e-3);
    static_assert( fabs(as.at(0)-5.) < 1e-3);
    static_assert( fabs(vs.at(1)-1.5) < 1e-3);
    static_assert( fabs(as.at(1)-0.) < 1e-3);
    static_assert( fabs(vs.at(2)-1.5) < 1e-3);
    static_assert( fabs(as.at(2)-0.) < 1e-3);
    static_assert( fabs(vs.at(3)-0.83) < 1e-2);
    static_assert( fabs(as.at(3)-(-5.)) < 1e-3);
}