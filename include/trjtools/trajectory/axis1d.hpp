#ifndef TRJTOOLS_TRAJECTORY_AXIS1D_HPP
#define TRJTOOLS_TRAJECTORY_AXIS1D_HPP

#include "trjtools/utility.hpp"
#include "nmtools/integration/trapezoidal.hpp"

#include <type_traits>
#include <vector>
#include <array>
#include <tuple>
#include <cmath>
#include <cassert>
#include <complex>

namespace trjtools {
    namespace trajectory {
        namespace axis1d {

            /**
             * @brief compute travel time (tx), constant vel period (tv), and acceleration period (ta)
             *  for second order trajectory profile
             * @note this routine assume the initial position and derivatives are zero
             */
            template <typename ScalarPos, typename ScalarVel, typename ScalarAcc>
            constexpr auto second_order(const ScalarPos &xbar, const ScalarVel &vbar, const ScalarAcc &abar)
            {
                assert(xbar > 0. && vbar > 0. && abar > 0.);
                auto t_abar = sqrt(xbar / abar);
                auto vhat = abar * t_abar;
                if (vhat > vbar) {
                    t_abar = vbar / abar;
                    vhat = vbar;
                }
                auto x_abar = 2 * ( 0.5 * abar * pow(t_abar, 2));
                auto t_vbar = (xbar - x_abar) / vbar;
                auto xbar_ = abar * pow(t_abar, 2) + vbar * t_vbar;
                auto t_xbar = 2 * t_abar + t_vbar;
                return std::make_tuple(t_xbar, t_vbar, t_abar);
            } // second_order

            /**
             * @brief compute second order profile from computed time coeff.
             * @note this routine assumes : 
             *  - the time is computed from 0
             *  - initial velocity is 0
             *  - initial position is 0
             *  - final velocity is 0
             */
            template <typename ScalarPos, typename ScalarVel, typename ScalarAcc, 
                typename Scalar, typename Iterable>
            constexpr auto second_order_profile(
                const ScalarPos &xbar, const ScalarVel &vbar, const ScalarAcc &abar, 
                const Scalar &tx, const Scalar &tv, const Scalar &ta, const Iterable &time_points)
            {
                static_assert(utility::traits::is_std_array_or_vector<Iterable>::value,
                    "currently only support std::array or std::vector for Iterable (time_points)"
                );
                utility::mpl::copy_std_container_t<Iterable,ScalarPos> xs{};
                utility::mpl::copy_std_container_t<Iterable,ScalarVel> vs{};
                utility::mpl::copy_std_container_t<Iterable,ScalarAcc> as{};
                auto tv0 = ta;
                auto tv1 = ta + tv;
                auto x_tv0 = 0.5 * abar * tv0 * tv0;
                auto x_tv1 = x_tv0 + vbar * tv;
                assert(tv < tx && ta < tx);
                if constexpr (utility::traits::is_std_vector<Iterable>::value) {
                    xs.resize(time_points.size());
                    vs.resize(time_points.size());
                    as.resize(time_points.size());
                }
                for (size_t i=0; i<std::size(time_points); i++) {
                    auto tp = time_points[i];
                    assert(tp < tx);
                    if (tp < tv0) {
                        as[i] = abar;
                        vs[i] = abar * tp;
                        xs[i] = 0.5 * vs[i] * tp;
                    } else if (tp < tv1) {
                        as[i] = 0.;
                        vs[i] = vbar;
                        xs[i] = x_tv0 + vbar * (tp - tv0);
                    } else {
                        as[i] = -abar;
                        vs[i] = vbar - abar * (tp - tv1);
                        xs[i] = x_tv1 + 0.5 * vs[i] * (tp - tv1);
                    }
                }
                return std::make_tuple(xs,vs,as);
            } // second_order_profile

            namespace helper {
                constexpr auto solver_4th_3rd_poly(auto t_dbar, auto xbar, auto dbar) {
                    using namespace std::complex_literals;
                    auto I = 1i;
                    auto root1 = -1.0/3.0*pow(t_dbar, 2)/cbrt(-55*pow(t_dbar, 3) 
                        + (1.0/2.0)*sqrt(-4*pow(t_dbar, 6) + pow(-110*pow(t_dbar, 3) 
                        + (27.0/2.0)*(8*dbar*pow(t_dbar, 4) - xbar)/(dbar*t_dbar), 2)) 
                        + (27.0/4.0)*(8*dbar*pow(t_dbar, 4) - xbar)/(dbar*t_dbar)) 
                        - 5.0/3.0*t_dbar - 1.0/3.0*cbrt(-55*pow(t_dbar, 3) 
                        + (1.0/2.0)*sqrt(-4*pow(t_dbar, 6) + pow(-110*pow(t_dbar, 3) 
                        + (27.0/2.0)*(8*dbar*pow(t_dbar, 4) - xbar)/(dbar*t_dbar), 2)) 
                        + (27.0/4.0)*(8*dbar*pow(t_dbar, 4) - xbar)/(dbar*t_dbar));
                    auto root2 = -1.0/3.0*pow(t_dbar, 2)/((-1.0/2.0 - 1.0/2.0*sqrt(3)*I)*cbrt(-55*pow(t_dbar, 3) 
                        + (1.0/2.0)*sqrt(-4*pow(t_dbar, 6) + pow(-110*pow(t_dbar, 3) 
                        + (27.0/2.0)*(8*dbar*pow(t_dbar, 4) - xbar)/(dbar*t_dbar), 2)) 
                        + (27.0/4.0)*(8*dbar*pow(t_dbar, 4) - xbar)/(dbar*t_dbar))) - 5.0/3.0*t_dbar 
                        - 1.0/3.0*(-1.0/2.0 - 1.0/2.0*sqrt(3)*I)*cbrt(-55*pow(t_dbar, 3) 
                        + (1.0/2.0)*sqrt(-4*pow(t_dbar, 6) + pow(-110*pow(t_dbar, 3) 
                        + (27.0/2.0)*(8*dbar*pow(t_dbar, 4) - xbar)/(dbar*t_dbar), 2)) 
                        + (27.0/4.0)*(8*dbar*pow(t_dbar, 4) - xbar)/(dbar*t_dbar));
                    auto root3 = -1.0/3.0*pow(t_dbar, 2)/((-1.0/2.0 
                        + (1.0/2.0)*sqrt(3)*I)*cbrt(-55*pow(t_dbar, 3) 
                        + (1.0/2.0)*sqrt(-4*pow(t_dbar, 6) + pow(-110*pow(t_dbar, 3) 
                        + (27.0/2.0)*(8*dbar*pow(t_dbar, 4) - xbar)/(dbar*t_dbar), 2)) 
                        + (27.0/4.0)*(8*dbar*pow(t_dbar, 4) - xbar)/(dbar*t_dbar))) - 5.0/3.0*t_dbar 
                        - 1.0/3.0*(-1.0/2.0 + (1.0/2.0)*sqrt(3)*I)*cbrt(-55*pow(t_dbar, 3) 
                        + (1.0/2.0)*sqrt(-4*pow(t_dbar, 6) + pow(-110*pow(t_dbar, 3) 
                        + (27.0/2.0)*(8*dbar*pow(t_dbar, 4) - xbar)/(dbar*t_dbar), 2)) 
                        + (27.0/4.0)*(8*dbar*pow(t_dbar, 4) - xbar)/(dbar*t_dbar));
                    return std::make_tuple(root1, std::complex<decltype(xbar)>{root2}, std::complex<decltype(xbar)>{root3});
                } // solver_4th_3rd_poly

                constexpr auto solver_4th_2nd_poly(auto t_dbar, auto vbar, auto dbar)
                {
                    auto root1 = (1.0/2.0)*(-3*dbar*pow(t_dbar, 2) + sqrt(dbar*t_dbar*(dbar*pow(t_dbar, 3) + 4*vbar)))/(dbar*t_dbar);
                    auto root2 = -1.0/2.0*(3*dbar*pow(t_dbar, 2) + sqrt(dbar*t_dbar*(dbar*pow(t_dbar, 3) + 4*vbar)))/(dbar*t_dbar);
                    return std::make_tuple(root1, root2);
                } // solver_4th_2nd_poly

                constexpr auto solver_4th_2nd_poly_ta(auto t_jbar, auto t_dbar, auto xbar, auto dbar)
                {
                    auto root1 = (1.0/2.0)*(-3*dbar*pow(t_dbar, 2)*(2*t_dbar + 3*t_jbar) + sqrt(dbar*t_dbar*(4*dbar*pow(t_dbar, 5) 
                        + 12*dbar*pow(t_dbar, 4)*t_jbar - 23*dbar*pow(t_dbar, 3)*pow(t_jbar, 2) - 48*dbar*pow(t_dbar, 2)*pow(t_jbar, 3) 
                        - 8*dbar*t_dbar*pow(t_jbar, 4) + 4*t_dbar*xbar + 4*t_jbar*xbar)))/(dbar*t_dbar*(t_dbar + t_jbar));
                    auto root2 = -1.0/2.0*(6*dbar*pow(t_dbar, 3) + 9*dbar*pow(t_dbar, 2)*t_jbar + sqrt(dbar*t_dbar*(4*dbar*pow(t_dbar, 5) 
                        + 12*dbar*pow(t_dbar, 4)*t_jbar - 23*dbar*pow(t_dbar, 3)*pow(t_jbar, 2) - 48*dbar*pow(t_dbar, 2)*pow(t_jbar, 3) 
                        - 8*dbar*t_dbar*pow(t_jbar, 4) + 4*t_dbar*xbar + 4*t_jbar*xbar)))/(dbar*t_dbar*(t_dbar + t_jbar));
                    return std::make_tuple(root1, root2);
                } // solver_4th_2nd_poly_ta

                /**
                 * @brief construct lambda defining the trajectory of the deriv. of the jerk
                 * 
                 * @tparam ScalarDJerk 
                 * @tparam Iterable 
                 * @param dbar constraint of djerk
                 * @param time_profile solved time profile, should have 16 size
                 * @return auto callable lambda
                 */
                template <typename ScalarDJerk, typename Iterable>
                constexpr auto fourth_order_djerk_profile(const ScalarDJerk &dbar, const Iterable &time_profile)
                {
                    /* define function of the derivative of the jerk */
                    assert(std::size(time_profile)==16);
                    return [=](const auto &t) {
                        ScalarDJerk d{0.};
                        size_t idx{0};
                        if ((t < time_profile[0]) || (t > time_profile[15]))
                            idx = 0;
                        else {
                            for (size_t i=1; i<std::size(time_profile); i++)
                                if (t <= time_profile[i]) {
                                    idx = i;
                                    break;
                                }
                        }
                        switch (idx)
                        {
                        case 3: case 5: case 9: case 15:
                            d = -dbar;
                            break;
                        case 1: case 7: case 11: case 13:
                            d = dbar;
                            break;
                        default:
                            d = 0.;
                            break;
                        }
                        return d;
                    };
                } // fourth_order_djerk_profile

                /**
                 * @brief naive implementation of 4th order trajectory profiles
                 *      returns callables that define the trajectories
                 * 
                 * @tparam integral_segment=100 
                 * @tparam ScalarDJerk 
                 * @tparam Iterable 
                 * @param dbar maximum value of the derivative of the jerk
                 * @param time_profile the computed time profile of the trajectory
                 * @return auto tuple of size 4 of callable functions
                 */
                template <size_t integral_segment=100, typename ScalarDJerk, typename Iterable>
                auto fourth_order_profile_fn(const ScalarDJerk &dbar, const Iterable &time_profile)
                {
                    auto dj_fn = fourth_order_djerk_profile(dbar, time_profile);
                    /* define function for j, a, v, x */
                    auto j_fn = utility::integrate<integral_segment>(dj_fn);
                    auto a_fn = utility::integrate<integral_segment>(j_fn);
                    auto v_fn = utility::integrate<integral_segment>(a_fn);
                    auto x_fn = utility::integrate<integral_segment>(v_fn);
                    return std::make_tuple(x_fn,v_fn,a_fn,j_fn,dj_fn);
                } // fourth_order_profile_fn

                /**
                 * @brief template-overloaded version of 4th order trajectory func
                 *      returns callables that define the trajectories
                 * 
                 * @tparam time_segment 
                 * @tparam integral_segment 
                 * @tparam ScalarDJerk 
                 * @tparam Iterable 
                 * @param dbar maximum value of the derivative of the jerk
                 * @param time_profile the computed time profile of the trajectory
                 * @return auto tuple of size 4 of callable functions
                 */
                template <size_t time_segment, size_t integral_segment, typename ScalarDJerk, typename Iterable>
                constexpr auto fourth_order_profile_fn(const ScalarDJerk &dbar, const Iterable &time_profile)
                {
                    using split_container_t = utility::mpl::copy_std_container_t<
                        Iterable,typename Iterable::value_type, time_segment+1>;
                    auto dj_fn = fourth_order_djerk_profile(dbar, time_profile);
                    auto time_pts = utility::split<split_container_t>(time_profile.back(),time_segment);
                    /* define function for j, a, v, x */
                    auto dj_int = utility::interpolate(time_pts, 
                        utility::sample(dj_fn, time_pts));
                    auto j_fn   = utility::integrate<integral_segment>(dj_int);
                    auto j_int  = utility::interpolate(time_pts,
                        utility::sample(j_fn, time_pts));
                    auto a_fn   = utility::integrate<integral_segment>(j_int);
                    auto a_int  = utility::interpolate(time_pts,
                        utility::sample(a_fn, time_pts));
                    auto v_fn   = utility::integrate<integral_segment>(a_int);
                    auto v_int  = utility::interpolate(time_pts,
                        utility::sample(v_fn,time_pts));
                    auto x_fn   = utility::integrate<integral_segment>(v_int);
                    return std::make_tuple(x_fn,v_fn,a_fn,j_fn,dj_fn);
                } // fourth_order_profile_fn
            } // namespace helper

            /**
             * @brief compute the time coefficient that define the trajectory while 
             *      satisfying the constraints
             * 
             * @tparam ScalarPos 
             * @tparam ScalarVel 
             * @tparam ScalarAcc 
             * @tparam ScalarJerk 
             * @tparam ScalarDJerk 
             * @param xbar target position 
             * @param vbar allowed maximum value of velocity
             * @param abar allowed maximum value of acceleration
             * @param jbar allowed maximum value of jerk
             * @param dbar allowed maximum value of the derivative of the jerk
             * @return constexpr auto 
             */
            template <typename ScalarPos, typename ScalarVel, typename ScalarAcc,
                typename ScalarJerk, typename ScalarDJerk>
            constexpr auto fourth_order(const ScalarPos &xbar, const ScalarVel &vbar, 
                const ScalarAcc &abar, const ScalarJerk &jbar, const ScalarDJerk &dbar)
            {
                auto t_dbar = sqrt(sqrt(xbar/(8*dbar)));
                auto v_hat = 2 * dbar * pow(t_dbar,3);
                if (v_hat > vbar) {
                    t_dbar = cbrt(vbar / (2*dbar));
                }
                auto a_hat = dbar * pow(t_dbar,2);
                if (a_hat > abar) {
                    t_dbar = sqrt(abar / dbar);
                }
                auto j_hat = dbar * t_dbar;
                if (j_hat > jbar) {
                    t_dbar = jbar / dbar;
                }
                auto [tjr0, tjr1, tjr2] = helper::solver_4th_3rd_poly(t_dbar, xbar, dbar);
                assert(tjr0 > 0. || tjr1.real() > 0. || tjr2.real() > 0.);
                auto t_jbar = tjr0;
                if (t_jbar < 0.)
                    t_jbar = (tjr1.real() > 0. ? tjr1.real() : tjr2.real());
                v_hat = 2*dbar*pow(t_dbar,3) + 3*dbar*pow(t_dbar,2)*t_jbar + dbar*t_dbar*pow(t_jbar,2);
                if (v_hat > vbar) {
                    auto [tjr0, tjr1] = helper::solver_4th_2nd_poly(t_dbar, vbar, dbar);
                    assert(tjr0 > 0. || tjr1 > 0.);
                    t_jbar = (tjr0 > 0. ? tjr0 : tjr1);
                }
                a_hat = dbar*pow(t_dbar,2) + dbar*t_dbar*t_jbar;
                if (a_hat > abar) 
                    t_jbar = (abar / jbar) - t_dbar;
                auto [tar0, tar1] = helper::solver_4th_2nd_poly_ta(t_jbar, t_dbar, xbar, dbar);
                assert(tar0 > 0. || tar1 > 0.);
                auto t_abar = (tar0 > 0. ? tar0 : tar1);
                v_hat = 2*dbar*pow(t_dbar,3) + 3*dbar*pow(t_dbar,2)*t_jbar + dbar*t_dbar*pow(t_jbar,2) + dbar*pow(t_dbar,2)*t_abar + dbar*t_dbar*t_jbar*t_abar;
                if (v_hat > vbar)
                    t_abar = (vbar-2*dbar*pow(t_dbar,3)-3*dbar*pow(t_dbar,2)*t_jbar-dbar*t_dbar*pow(t_jbar,2)) 
                            / (dbar*pow(t_dbar,2) + dbar*t_dbar*t_jbar);
                auto x_abar = (8*pow(t_dbar,4)+16*pow(t_dbar,3)*t_jbar+10*pow(t_dbar,2)*pow(t_jbar,2)+2*t_dbar*2*pow(t_jbar,3)
                                + pow(t_dbar,2)*pow(t_abar,2) + t_dbar*t_jbar*pow(t_abar,2) + 6*pow(t_dbar,3)*t_abar 
                                + 9*pow(t_dbar,2)*t_jbar*t_abar + 3*t_dbar*pow(t_jbar,2)*t_abar) * dbar;
                auto t_vbar = (xbar - x_abar) / vbar;
                return std::make_tuple(t_vbar, t_abar, t_jbar, t_dbar);
            } // fourth_order

            /**
             * @brief construct array of time profile defining the trajectory
             * 
             * @tparam Scalar 
             * @param t_vbar computed value of constant velocity
             * @param t_abar computed value of constant acceleration
             * @param t_jbar computed value of constant jerk
             * @param t_dbar computed value of time swicth of the derivative of the jerk
             * @return constexpr auto array of Scalar with size of 16
             */
            template <typename Scalar>
            constexpr auto fourth_order_time_points(const Scalar &t_vbar, 
                const Scalar &t_abar, const Scalar &t_jbar, const Scalar &t_dbar)
            {
                std::array<Scalar,16> t{};
                t[0]  = 0.;
                t[1]  = t[0]  + t_dbar;
                t[2]  = t[1]  + t_jbar;
                t[3]  = t[2]  + t_dbar;
                t[4]  = t[3]  + t_abar;
                t[5]  = t[4]  + t_dbar;
                t[6]  = t[5]  + t_jbar;
                t[7]  = t[6]  + t_dbar;
                t[8]  = t[7]  + t_vbar;
                t[9]  = t[8]  + t_dbar;
                t[10] = t[9]  + t_jbar;
                t[11] = t[10] + t_dbar;
                t[12] = t[11] + t_abar;
                t[13] = t[12] + t_dbar;
                t[14] = t[13] + t_jbar;
                t[15] = t[14] + t_dbar;
                return t;
            }

            /**
             * @brief compute the 4th order trajectory given constraints and time profile
             * 
             * @tparam integral_segment=100 
             * @tparam ScalarPos 
             * @tparam ScalarVel 
             * @tparam ScalarAcc 
             * @tparam ScalarJerk 
             * @tparam ScalarDJerk 
             * @tparam Scalar 
             * @tparam Iterable 
             * @param xbar target position
             * @param vbar maximum allowed velocity
             * @param abar maximum allowed acceleration
             * @param jbar maximum allowed jerk
             * @param dbar maximum allowed derivative of the jerk
             * @param t_vbar computed constant velocity time profile
             * @param t_abar computed constant acceleration time profile
             * @param t_jbar computed constant jerk time profile
             * @param t_dbar computed time switch for djerk
             * @param time_points 
             * @return auto tuple of size 5 of container defining the trajectory
             */
            template <size_t integral_segment=100, typename ScalarPos, typename ScalarVel, typename ScalarAcc,
                typename ScalarJerk, typename ScalarDJerk, typename Scalar, typename Iterable>
            auto fourth_order_profile(const ScalarPos &xbar, const ScalarVel &vbar, 
                const ScalarAcc &abar, const ScalarJerk &jbar, const ScalarDJerk &dbar, 
                const Scalar &t_vbar, const Scalar &t_abar, const Scalar &t_jbar, 
                const Scalar &t_dbar, const Iterable &time_points)
            {
                static_assert(utility::traits::is_std_array_or_vector<Iterable>::value,
                    "currently only support std::array or std::vector for Iterable (time_points)"
                );
                utility::mpl::copy_std_container_t<Iterable,ScalarPos>   xs;
                utility::mpl::copy_std_container_t<Iterable,ScalarVel>   vs;
                utility::mpl::copy_std_container_t<Iterable,ScalarAcc>   as;
                utility::mpl::copy_std_container_t<Iterable,ScalarJerk>  js;
                utility::mpl::copy_std_container_t<Iterable,ScalarDJerk> ds;
                auto tp = fourth_order_time_points(t_vbar, t_abar, t_jbar, t_dbar);
                if constexpr (utility::traits::is_std_vector<Iterable>::value) {
                    xs.resize(time_points.size());
                    vs.resize(time_points.size());
                    as.resize(time_points.size());
                    js.resize(time_points.size());
                    ds.resize(time_points.size());
                }
                auto [x_fn, v_fn, a_fn, j_fn, dj_fn] = helper::fourth_order_profile_fn<integral_segment>(dbar, tp);
                // auto dj_fn = helper::fourth_order_djerk_profile(dbar, tp);

                for (size_t i=0; i<std::size(time_points); i++) {
                    auto t = time_points[i];
                    assert(t > tp[0] && t < tp[15]);
                    ds[i] = dj_fn(t);
                    js[i] = j_fn(t);
                    as[i] = a_fn(t);
                    vs[i] = v_fn(t);
                    xs[i] = x_fn(t);
                }
                return std::make_tuple(xs,vs,as,js,ds);
            }
        } // namespace axis1d
    } // namespace trajectory
} // namespace trjtools

#endif // TRJTOOLS_TRAJECTORY_AXIS1D_HPP