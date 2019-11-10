#include <gtest/gtest.h>
#include "trjtools/trajectory.hpp"

namespace trj = trjtools::trajectory;

TEST(axis3d, second_order)
{
    double abar = 5.;
    double vbar = 1.5;
    {
        double xbar = 1.;
        double ybar = 1.;
        double zbar = 1.;
        auto [tx, ty, tz] = trj::axis3d::second_order(xbar,ybar,zbar,vbar,abar);
        auto [tpx, tvx, tax] = tx;
        auto [tpy, tvy, tay] = ty;
        auto [tpz, tvz, taz] = tz;
        EXPECT_NEAR(tpx,tpy,1e-6);
        EXPECT_NEAR(tpx,tpz,1e-6);
    }
    {
        double xbar = .3;
        double ybar = 1.;
        double zbar = 1.;
        auto [tx, ty, tz] = trj::axis3d::second_order(xbar,ybar,zbar,vbar,abar);
        auto [tpx, tvx, tax] = tx;
        auto [tpy, tvy, tay] = ty;
        auto [tpz, tvz, taz] = tz;
        EXPECT_NEAR(tpx,tpy,1e-6);
        EXPECT_NEAR(tpx,tpz,1e-6);
    }
    {
        double xbar = 1.;
        double ybar = .3;
        double zbar = 1.;
        auto [tx, ty, tz] = trj::axis3d::second_order(xbar,ybar,zbar,vbar,abar);
        auto [tpx, tvx, tax] = tx;
        auto [tpy, tvy, tay] = ty;
        auto [tpz, tvz, taz] = tz;
        EXPECT_NEAR(tpx,tpy,1e-6);
        EXPECT_NEAR(tpx,tpz,1e-6);
    }
    {
        double xbar = 1.;
        double ybar = 1.;
        double zbar = .3;
        auto [tx, ty, tz] = trj::axis3d::second_order(xbar,ybar,zbar,vbar,abar);
        auto [tpx, tvx, tax] = tx;
        auto [tpy, tvy, tay] = ty;
        auto [tpz, tvz, taz] = tz;
        EXPECT_NEAR(tpx,tpy,1e-6);
        EXPECT_NEAR(tpx,tpz,1e-6);
    }
    {
        double xbar = .3;
        double ybar = .3;
        double zbar = 1.;
        auto [tx, ty, tz] = trj::axis3d::second_order(xbar,ybar,zbar,vbar,abar);
        auto [tpx, tvx, tax] = tx;
        auto [tpy, tvy, tay] = ty;
        auto [tpz, tvz, taz] = tz;
        EXPECT_NEAR(tpx,tpy,1e-6);
        EXPECT_NEAR(tpx,tpz,1e-6);
    }
    {
        double xbar = 1.;
        double ybar = .3;
        double zbar = .3;
        auto [tx, ty, tz] = trj::axis3d::second_order(xbar,ybar,zbar,vbar,abar);
        auto [tpx, tvx, tax] = tx;
        auto [tpy, tvy, tay] = ty;
        auto [tpz, tvz, taz] = tz;
        EXPECT_NEAR(tpx,tpy,1e-6);
        EXPECT_NEAR(tpx,tpz,1e-6);
    }
}

TEST(axis3d, constexpr_second_order)
{
    constexpr double abar = 5.;
    constexpr double vbar = 1.5;
    {
        constexpr double xbar = 1.;
        constexpr double ybar = 1.;
        constexpr double zbar = 1.;
        constexpr auto time = trj::axis3d::second_order(xbar,ybar,zbar,vbar,abar);
        constexpr auto tx = std::get<0>(time);
        constexpr auto ty = std::get<1>(time);
        constexpr auto tz = std::get<2>(time);
        constexpr auto tpx = std::get<0>(tx);
        constexpr auto tpy = std::get<0>(ty);
        constexpr auto tpz = std::get<0>(tz);
        static_assert( fabs(tpx-tpy) < 1e-6);
        static_assert( fabs(tpx-tpz) < 1e-6);
    }
    {
        constexpr double xbar = .3;
        constexpr double ybar = 1.;
        constexpr double zbar = 1.;
        constexpr auto time = trj::axis3d::second_order(xbar,ybar,zbar,vbar,abar);
        constexpr auto tx = std::get<0>(time);
        constexpr auto ty = std::get<1>(time);
        constexpr auto tz = std::get<2>(time);
        constexpr auto tpx = std::get<0>(tx);
        constexpr auto tpy = std::get<0>(ty);
        constexpr auto tpz = std::get<0>(tz);
        static_assert( fabs(tpx-tpy) < 1e-6);
        static_assert( fabs(tpx-tpz) < 1e-6);
    }
    {
        constexpr double xbar = 1.;
        constexpr double ybar = .3;
        constexpr double zbar = 1.;
        constexpr auto time = trj::axis3d::second_order(xbar,ybar,zbar,vbar,abar);
        constexpr auto tx = std::get<0>(time);
        constexpr auto ty = std::get<1>(time);
        constexpr auto tz = std::get<2>(time);
        constexpr auto tpx = std::get<0>(tx);
        constexpr auto tpy = std::get<0>(ty);
        constexpr auto tpz = std::get<0>(tz);
        static_assert( fabs(tpx-tpy) < 1e-6);
        static_assert( fabs(tpx-tpz) < 1e-6);
    }
    {
        constexpr double xbar = 1.;
        constexpr double ybar = 1.;
        constexpr double zbar = .3;
        constexpr auto time = trj::axis3d::second_order(xbar,ybar,zbar,vbar,abar);
        constexpr auto tx = std::get<0>(time);
        constexpr auto ty = std::get<1>(time);
        constexpr auto tz = std::get<2>(time);
        constexpr auto tpx = std::get<0>(tx);
        constexpr auto tpy = std::get<0>(ty);
        constexpr auto tpz = std::get<0>(tz);
        static_assert( fabs(tpx-tpy) < 1e-6);
        static_assert( fabs(tpx-tpz) < 1e-6);
    }
    {
        constexpr double xbar = .3;
        constexpr double ybar = .3;
        constexpr double zbar = 1.;
        constexpr auto time = trj::axis3d::second_order(xbar,ybar,zbar,vbar,abar);
        constexpr auto tx = std::get<0>(time);
        constexpr auto ty = std::get<1>(time);
        constexpr auto tz = std::get<2>(time);
        constexpr auto tpx = std::get<0>(tx);
        constexpr auto tpy = std::get<0>(ty);
        constexpr auto tpz = std::get<0>(tz);
        static_assert( fabs(tpx-tpy) < 1e-6);
        static_assert( fabs(tpx-tpz) < 1e-6);
    }
    {
        constexpr double xbar = 1.;
        constexpr double ybar = .3;
        constexpr double zbar = .3;
        constexpr auto time = trj::axis3d::second_order(xbar,ybar,zbar,vbar,abar);
        constexpr auto tx = std::get<0>(time);
        constexpr auto ty = std::get<1>(time);
        constexpr auto tz = std::get<2>(time);
        constexpr auto tpx = std::get<0>(tx);
        constexpr auto tpy = std::get<0>(ty);
        constexpr auto tpz = std::get<0>(tz);
        static_assert( fabs(tpx-tpy) < 1e-6);
        static_assert( fabs(tpx-tpz) < 1e-6);
    }
}

TEST(axis3d, fourth_order)
{
    double dbar = 1000.;
    double jbar = 50.;
    double abar = 5.;
    double vbar = 1.5;
    {
        double xbar = 1.;
        double ybar = 1.;
        double zbar = 1.;
        auto [tx, ty, tz] = trj::axis3d::fourth_order(xbar,ybar,zbar,vbar,abar,jbar, dbar);
        auto [tvx, tax, tjx, tdx] = tx;
        auto [tvy, tay, tjy, tdy] = ty;
        auto [tvz, taz, tjz, tdz] = tz;
        auto tpx = trj::axis1d::fourth_order_time_points(tvx, tax, tjx, tdx);
        auto tpy = trj::axis1d::fourth_order_time_points(tvy, tay, tjy, tdy);
        auto tpz = trj::axis1d::fourth_order_time_points(tvz, taz, tjz, tdz);
        EXPECT_NEAR(tpx.back(),tpy.back(),1e-6);
        EXPECT_NEAR(tpx.back(),tpz.back(),1e-6);
    }
    {
        double xbar = 1.;
        double ybar = .3;
        double zbar = 1.;
        auto [tx, ty, tz] = trj::axis3d::fourth_order(xbar,ybar,zbar,vbar,abar,jbar, dbar);
        auto [tvx, tax, tjx, tdx] = tx;
        auto [tvy, tay, tjy, tdy] = ty;
        auto [tvz, taz, tjz, tdz] = tz;
        auto tpx = trj::axis1d::fourth_order_time_points(tvx, tax, tjx, tdx);
        auto tpy = trj::axis1d::fourth_order_time_points(tvy, tay, tjy, tdy);
        auto tpz = trj::axis1d::fourth_order_time_points(tvz, taz, tjz, tdz);
        EXPECT_NEAR(tpx.back(),tpy.back(),1e-6);
        EXPECT_NEAR(tpx.back(),tpz.back(),1e-6);
    }
    {
        double xbar = 1.;
        double ybar = 1.;
        double zbar = 3.;
        auto [tx, ty, tz] = trj::axis3d::fourth_order(xbar,ybar,zbar,vbar,abar,jbar, dbar);
        auto [tvx, tax, tjx, tdx] = tx;
        auto [tvy, tay, tjy, tdy] = ty;
        auto [tvz, taz, tjz, tdz] = tz;
        auto tpx = trj::axis1d::fourth_order_time_points(tvx, tax, tjx, tdx);
        auto tpy = trj::axis1d::fourth_order_time_points(tvy, tay, tjy, tdy);
        auto tpz = trj::axis1d::fourth_order_time_points(tvz, taz, tjz, tdz);
        EXPECT_NEAR(tpx.back(),tpy.back(),1e-6);
        EXPECT_NEAR(tpx.back(),tpz.back(),1e-6);
    }
    {
        double xbar = 3.;
        double ybar = .3;
        double zbar = 1.;
        auto [tx, ty, tz] = trj::axis3d::fourth_order(xbar,ybar,zbar,vbar,abar,jbar, dbar);
        auto [tvx, tax, tjx, tdx] = tx;
        auto [tvy, tay, tjy, tdy] = ty;
        auto [tvz, taz, tjz, tdz] = tz;
        auto tpx = trj::axis1d::fourth_order_time_points(tvx, tax, tjx, tdx);
        auto tpy = trj::axis1d::fourth_order_time_points(tvy, tay, tjy, tdy);
        auto tpz = trj::axis1d::fourth_order_time_points(tvz, taz, tjz, tdz);
        EXPECT_NEAR(tpx.back(),tpy.back(),1e-6);
        EXPECT_NEAR(tpx.back(),tpz.back(),1e-6);
    }
    {
        double xbar = 1.;
        double ybar = .3;
        double zbar = 3.;
        auto [tx, ty, tz] = trj::axis3d::fourth_order(xbar,ybar,zbar,vbar,abar,jbar, dbar);
        auto [tvx, tax, tjx, tdx] = tx;
        auto [tvy, tay, tjy, tdy] = ty;
        auto [tvz, taz, tjz, tdz] = tz;
        auto tpx = trj::axis1d::fourth_order_time_points(tvx, tax, tjx, tdx);
        auto tpy = trj::axis1d::fourth_order_time_points(tvy, tay, tjy, tdy);
        auto tpz = trj::axis1d::fourth_order_time_points(tvz, taz, tjz, tdz);
        EXPECT_NEAR(tpx.back(),tpy.back(),1e-6);
        EXPECT_NEAR(tpx.back(),tpz.back(),1e-6);
    }
}

TEST(axis3d, constexpr_fourth_order)
{
    constexpr double dbar = 1000.;
    constexpr double jbar = 50.;
    constexpr double abar = 5.;
    constexpr double vbar = 1.5;
    {
        constexpr double xbar = 1.;
        constexpr double ybar = 1.;
        constexpr double zbar = 1.;
        constexpr auto time = trj::axis3d::fourth_order(xbar,ybar,zbar,vbar,abar,jbar, dbar);
        constexpr auto tx = std::get<0>(time);
        constexpr auto ty = std::get<1>(time);
        constexpr auto tz = std::get<2>(time);
        constexpr auto tvx = std::get<0>(tx);
        constexpr auto tax = std::get<1>(tx);
        constexpr auto tjx = std::get<2>(tx);
        constexpr auto tdx = std::get<3>(tx);
        constexpr auto tvy = std::get<0>(ty);
        constexpr auto tay = std::get<1>(ty);
        constexpr auto tjy = std::get<2>(ty);
        constexpr auto tdy = std::get<3>(ty);
        constexpr auto tvz = std::get<0>(tz);
        constexpr auto taz = std::get<1>(tz);
        constexpr auto tjz = std::get<2>(tz);
        constexpr auto tdz = std::get<3>(tz);
        constexpr auto tpx = trj::axis1d::fourth_order_time_points(tvx, tax, tjx, tdx);
        constexpr auto tpy = trj::axis1d::fourth_order_time_points(tvy, tay, tjy, tdy);
        constexpr auto tpz = trj::axis1d::fourth_order_time_points(tvz, taz, tjz, tdz);
        static_assert( fabs(tpx.back()-tpy.back()) < 1e-6);
        static_assert( fabs(tpx.back()-tpz.back()) < 1e-6);
    }
    {
        constexpr double xbar = 1.;
        constexpr double ybar = .3;
        constexpr double zbar = 1.;
        constexpr auto time = trj::axis3d::fourth_order(xbar,ybar,zbar,vbar,abar,jbar, dbar);
        constexpr auto tx = std::get<0>(time);
        constexpr auto ty = std::get<1>(time);
        constexpr auto tz = std::get<2>(time);
        constexpr auto tvx = std::get<0>(tx);
        constexpr auto tax = std::get<1>(tx);
        constexpr auto tjx = std::get<2>(tx);
        constexpr auto tdx = std::get<3>(tx);
        constexpr auto tvy = std::get<0>(ty);
        constexpr auto tay = std::get<1>(ty);
        constexpr auto tjy = std::get<2>(ty);
        constexpr auto tdy = std::get<3>(ty);
        constexpr auto tvz = std::get<0>(tz);
        constexpr auto taz = std::get<1>(tz);
        constexpr auto tjz = std::get<2>(tz);
        constexpr auto tdz = std::get<3>(tz);
        constexpr auto tpx = trj::axis1d::fourth_order_time_points(tvx, tax, tjx, tdx);
        constexpr auto tpy = trj::axis1d::fourth_order_time_points(tvy, tay, tjy, tdy);
        constexpr auto tpz = trj::axis1d::fourth_order_time_points(tvz, taz, tjz, tdz);
        static_assert( fabs(tpx.back()-tpy.back()) < 1e-6);
        static_assert( fabs(tpx.back()-tpz.back()) < 1e-6);
    }
    {
        constexpr double xbar = .3;
        constexpr double ybar = .3;
        constexpr double zbar = 1.;
        constexpr auto time = trj::axis3d::fourth_order(xbar,ybar,zbar,vbar,abar,jbar, dbar);
        constexpr auto tx = std::get<0>(time);
        constexpr auto ty = std::get<1>(time);
        constexpr auto tz = std::get<2>(time);
        constexpr auto tvx = std::get<0>(tx);
        constexpr auto tax = std::get<1>(tx);
        constexpr auto tjx = std::get<2>(tx);
        constexpr auto tdx = std::get<3>(tx);
        constexpr auto tvy = std::get<0>(ty);
        constexpr auto tay = std::get<1>(ty);
        constexpr auto tjy = std::get<2>(ty);
        constexpr auto tdy = std::get<3>(ty);
        constexpr auto tvz = std::get<0>(tz);
        constexpr auto taz = std::get<1>(tz);
        constexpr auto tjz = std::get<2>(tz);
        constexpr auto tdz = std::get<3>(tz);
        constexpr auto tpx = trj::axis1d::fourth_order_time_points(tvx, tax, tjx, tdx);
        constexpr auto tpy = trj::axis1d::fourth_order_time_points(tvy, tay, tjy, tdy);
        constexpr auto tpz = trj::axis1d::fourth_order_time_points(tvz, taz, tjz, tdz);
        static_assert( fabs(tpx.back()-tpy.back()) < 1e-6);
        static_assert( fabs(tpx.back()-tpz.back()) < 1e-6);
    }
    {
        constexpr double xbar = 1.;
        constexpr double ybar = .3;
        constexpr double zbar = .3;
        constexpr auto time = trj::axis3d::fourth_order(xbar,ybar,zbar,vbar,abar,jbar, dbar);
        constexpr auto tx = std::get<0>(time);
        constexpr auto ty = std::get<1>(time);
        constexpr auto tz = std::get<2>(time);
        constexpr auto tvx = std::get<0>(tx);
        constexpr auto tax = std::get<1>(tx);
        constexpr auto tjx = std::get<2>(tx);
        constexpr auto tdx = std::get<3>(tx);
        constexpr auto tvy = std::get<0>(ty);
        constexpr auto tay = std::get<1>(ty);
        constexpr auto tjy = std::get<2>(ty);
        constexpr auto tdy = std::get<3>(ty);
        constexpr auto tvz = std::get<0>(tz);
        constexpr auto taz = std::get<1>(tz);
        constexpr auto tjz = std::get<2>(tz);
        constexpr auto tdz = std::get<3>(tz);
        constexpr auto tpx = trj::axis1d::fourth_order_time_points(tvx, tax, tjx, tdx);
        constexpr auto tpy = trj::axis1d::fourth_order_time_points(tvy, tay, tjy, tdy);
        constexpr auto tpz = trj::axis1d::fourth_order_time_points(tvz, taz, tjz, tdz);
        static_assert( fabs(tpx.back()-tpy.back()) < 1e-6);
        static_assert( fabs(tpx.back()-tpz.back()) < 1e-6);
    }
    {
        constexpr double xbar = 1.;
        constexpr double ybar = 1.;
        constexpr double zbar = .3;
        constexpr auto time = trj::axis3d::fourth_order(xbar,ybar,zbar,vbar,abar,jbar, dbar);
        constexpr auto tx = std::get<0>(time);
        constexpr auto ty = std::get<1>(time);
        constexpr auto tz = std::get<2>(time);
        constexpr auto tvx = std::get<0>(tx);
        constexpr auto tax = std::get<1>(tx);
        constexpr auto tjx = std::get<2>(tx);
        constexpr auto tdx = std::get<3>(tx);
        constexpr auto tvy = std::get<0>(ty);
        constexpr auto tay = std::get<1>(ty);
        constexpr auto tjy = std::get<2>(ty);
        constexpr auto tdy = std::get<3>(ty);
        constexpr auto tvz = std::get<0>(tz);
        constexpr auto taz = std::get<1>(tz);
        constexpr auto tjz = std::get<2>(tz);
        constexpr auto tdz = std::get<3>(tz);
        constexpr auto tpx = trj::axis1d::fourth_order_time_points(tvx, tax, tjx, tdx);
        constexpr auto tpy = trj::axis1d::fourth_order_time_points(tvy, tay, tjy, tdy);
        constexpr auto tpz = trj::axis1d::fourth_order_time_points(tvz, taz, tjz, tdz);
        static_assert( fabs(tpx.back()-tpy.back()) < 1e-6);
        static_assert( fabs(tpx.back()-tpz.back()) < 1e-6);
    }
}