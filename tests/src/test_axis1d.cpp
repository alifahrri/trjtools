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

TEST(axis1d, fourth_order)
{
    double dbar = 1000.;
    double jbar = 50.;
    double abar = 5.;
    double vbar = 1.5;
    double xbar = 1.;
    auto [tv, ta, tj, td] = trj::axis1d::fourth_order(xbar,vbar,abar,jbar,dbar);
    // EXPECT_NEAR(tv,0.3,1e-2);
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