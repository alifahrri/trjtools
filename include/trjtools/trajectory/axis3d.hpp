#ifndef TRJTOOLS_TRAJECTORY_AXIS3D_HPP
#define TRJTOOLS_TRAJECTORY_AXIS3D_HPP

namespace trjtools {
    namespace trajectory {
        namespace axis3d {
            template <typename ScalarPos, typename ScalarVel, typename ScalarAcc>
            constexpr auto second_order(const ScalarPos &xbar, const ScalarPos &ybar, 
                const ScalarPos &zbar, const ScalarVel &vbar, const ScalarAcc &abar)
            {
                auto xy_angle = atan2(ybar, xbar);
                auto zx_angle = atan2(zbar, hypot(ybar,xbar));
                auto xy_vbar = cos(zx_angle) * vbar;
                auto zvbar = sin(zx_angle) * vbar;
                auto xy_abar = cos(zx_angle) * abar;
                auto zabar = sin(zx_angle) * abar;
                auto xvbar = cos(xy_angle) * xy_vbar;
                auto yvbar = sin(xy_angle) * xy_vbar;
                auto xabar = cos(xy_angle) * xy_abar;
                auto yabar = sin(xy_angle) * xy_abar;
                return std::make_tuple(
                    axis1d::second_order(xbar,xvbar,xabar),
                    axis1d::second_order(ybar,yvbar,yabar),
                    axis1d::second_order(zbar,zvbar,zabar)
                );
            } // second_order

            template <typename ScalarPos, typename ScalarVel, typename ScalarAcc,
                typename ScalarJerk, typename ScalarDJerk>
            constexpr auto fourth_order(const ScalarPos &xbar, const ScalarPos &ybar, const ScalarPos &zbar, 
                const ScalarVel &vbar, const ScalarAcc &abar, const ScalarJerk &jbar, const ScalarDJerk &dbar)
            {
                auto xy_angle = atan2(ybar, xbar);
                auto zx_angle = atan2(zbar, hypot(ybar,xbar));
                auto xy_vbar = cos(zx_angle) * vbar;
                auto zvbar = sin(zx_angle) * vbar;
                auto xy_abar = cos(zx_angle) * abar;
                auto zabar = sin(zx_angle) * abar;
                auto xy_jbar = cos(zx_angle) * jbar;
                auto zjbar = sin(zx_angle) * jbar;
                auto xy_dbar = cos(zx_angle) * dbar;
                auto zdbar = sin(zx_angle) * dbar;
                auto xvbar = cos(xy_angle) * xy_vbar;
                auto yvbar = sin(xy_angle) * xy_vbar;
                auto xabar = cos(xy_angle) * xy_abar;
                auto yabar = sin(xy_angle) * xy_abar;
                auto xjbar = cos(xy_angle) * xy_jbar;
                auto yjbar = sin(xy_angle) * xy_jbar;
                auto xdbar = cos(xy_angle) * xy_dbar;
                auto ydbar = sin(xy_angle) * xy_dbar;
                return std::make_tuple(
                    axis1d::fourth_order(xbar,xvbar,xabar,xjbar,xdbar),
                    axis1d::fourth_order(ybar,yvbar,yabar,yjbar,ydbar),
                    axis1d::fourth_order(zbar,zvbar,zabar,zjbar,zdbar)
                );
            } // fourth_order
        } // namespace axis3d
    } // namespace trajectory
} // namespace trjtools

#endif