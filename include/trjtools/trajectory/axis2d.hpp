#ifndef TRJTOOLS_TRAJECTORY_AXIS2D_HPP
#define TRJTOOLS_TRAJECTORY_AXIS2D_HPP

#include "trjtools/trajectory/axis1d.hpp"

#include <type_traits>
#include <vector>
#include <array>
#include <tuple>
#include <cmath>
#include <cassert>

namespace trjtools {
    namespace trajectory {
        namespace axis2d {

            template <typename ScalarPos, typename ScalarVel, typename ScalarAcc>
            constexpr auto second_order(const ScalarPos &xbar, const ScalarPos &ybar, const ScalarVel &vbar, const ScalarAcc &abar)
            {
                auto angle = atan2(ybar, xbar);
                auto xvbar = cos(angle) * vbar;
                auto yvbar = sin(angle) * vbar;
                auto xabar = cos(angle) * abar;
                auto yabar = sin(angle) * abar;
                return std::make_tuple(
                    axis1d::second_order(xbar,xvbar,xabar),
                    axis1d::second_order(ybar,yvbar,yabar)
                );
            } // second_order

            template <typename ScalarPos, typename ScalarVel, typename ScalarAcc,
                typename ScalarJerk, typename ScalarDJerk>
            constexpr auto fourth_order(const ScalarPos &xbar, const ScalarPos &ybar, const ScalarVel &vbar, 
                const ScalarAcc &abar, const ScalarJerk &jbar, const ScalarDJerk &dbar)
            {
                auto angle = atan2(ybar, xbar);
                auto xvbar = cos(angle) * vbar;
                auto yvbar = sin(angle) * vbar;
                auto xabar = cos(angle) * abar;
                auto yabar = sin(angle) * abar;
                auto xjbar = cos(angle) * jbar;
                auto yjbar = sin(angle) * jbar;
                auto xdbar = cos(angle) * dbar;
                auto ydbar = sin(angle) * dbar;
                return std::make_tuple(
                    axis1d::fourth_order(xbar,xvbar,xabar,xjbar,xdbar),
                    axis1d::fourth_order(ybar,yvbar,yabar,yjbar,ydbar)
                );
            } // fourth_order

        } // namespace axis2d
    } // namespace trajectory
} // namespace trjtools
#endif