#ifndef TRJTOOLS_TRAJECTORY_AXIS1D_HPP
#define TRJTOOLS_TRAJECTORY_AXIS1D_HPP

#include "trjtools/utility.hpp"

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
            auto second_order(const ScalarPos &xbar, const ScalarVel &vbar, const ScalarAcc &abar)
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
            auto second_order_profile(
                const ScalarPos &xbar, const ScalarVel &vbar, const ScalarAcc &abar, 
                const Scalar &tx, const Scalar &tv, const Scalar &ta, const Iterable &time_points)
            {
                static_assert(utility::traits::is_std_array_or_vector<Iterable>::value,
                    "currently only support std::array or std::vector for Iterable (time_points)"
                );
                utility::mpl::copy_std_container_t<Iterable,ScalarPos> xs;
                utility::mpl::copy_std_container_t<Iterable,ScalarVel> vs;
                utility::mpl::copy_std_container_t<Iterable,ScalarAcc> as;
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
                auto solver_4th_3rd_poly(auto t_dbar, auto xbar, auto dbar) {
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

                auto solver_4th_2nd_poly(auto t_dbar, auto vbar, auto dbar)
                {
                    auto root1 = (1.0/2.0)*(-3*dbar*pow(t_dbar, 2) + sqrt(dbar*t_dbar*(dbar*pow(t_dbar, 3) + 4*vbar)))/(dbar*t_dbar);
                    auto root2 = -1.0/2.0*(3*dbar*pow(t_dbar, 2) + sqrt(dbar*t_dbar*(dbar*pow(t_dbar, 3) + 4*vbar)))/(dbar*t_dbar);
                    return std::make_tuple(root1, root2);
                }

                auto solver_4th_2nd_poly_ta(auto t_jbar, auto t_dbar, auto xbar, auto dbar)
                {
                    auto root1 = (1.0/2.0)*(-3*dbar*pow(t_dbar, 2)*(2*t_dbar + 3*t_jbar) + sqrt(dbar*t_dbar*(4*dbar*pow(t_dbar, 5) 
                        + 12*dbar*pow(t_dbar, 4)*t_jbar - 23*dbar*pow(t_dbar, 3)*pow(t_jbar, 2) - 48*dbar*pow(t_dbar, 2)*pow(t_jbar, 3) 
                        - 8*dbar*t_dbar*pow(t_jbar, 4) + 4*t_dbar*xbar + 4*t_jbar*xbar)))/(dbar*t_dbar*(t_dbar + t_jbar));
                    auto root2 = -1.0/2.0*(6*dbar*pow(t_dbar, 3) + 9*dbar*pow(t_dbar, 2)*t_jbar + sqrt(dbar*t_dbar*(4*dbar*pow(t_dbar, 5) 
                        + 12*dbar*pow(t_dbar, 4)*t_jbar - 23*dbar*pow(t_dbar, 3)*pow(t_jbar, 2) - 48*dbar*pow(t_dbar, 2)*pow(t_jbar, 3) 
                        - 8*dbar*t_dbar*pow(t_jbar, 4) + 4*t_dbar*xbar + 4*t_jbar*xbar)))/(dbar*t_dbar*(t_dbar + t_jbar));
                    return std::make_tuple(root1, root2);
                }
            }

            template <typename ScalarPos, typename ScalarVel, typename ScalarAcc,
                typename ScalarJerk, typename ScalarDJerk>
            auto fourth_order(const ScalarPos &xbar, const ScalarVel &vbar, 
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

            template <typename Scalar>
            auto fourth_order_time_points(const Scalar &t_vbar, 
                const Scalar &t_abar, const Scalar &t_jbar, const Scalar &t_dbar)
            {
                std::array<Scalar,16> t;
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

            template <typename ScalarPos, typename ScalarVel, typename ScalarAcc,
                typename ScalarJerk, typename ScalarDJerk, typename Scalar>
            auto fourth_order_pivot(const ScalarPos &xbar, const ScalarVel &vbar, 
                const ScalarAcc &abar, const ScalarJerk &jbar, const ScalarDJerk &dbar, 
                const std::array<Scalar,16> t)
            {
                std::array<ScalarPos,16>   xs;
                std::array<ScalarVel,16>   vs;
                std::array<ScalarAcc,16>   as;
                std::array<ScalarJerk,16>  js;
                std::array<ScalarDJerk,16> ds;
                ds[0] = 0.;     js[0] = 0.;     as[0] = 0.;
                ds[1] = dbar;   js[1] = jbar;   as[1] = as[0] + (1./2.) * jbar * (t[1]-t[0]);
                ds[2] = 0.;     js[2] = jbar;   as[2] = as[1] + jbar * (t[2]-t[1]);
                ds[3] = -dbar;  js[3] = 0.;     as[3] = abar;
                ds[4] = 0.;     js[4] = 0.;     as[4] = abar;
                ds[5] = -dbar;  js[5] = -jbar;  as[5] = as[4] - (1./2.) * jbar * (t[5]-t[4]);
                ds[6] = 0.;     js[6] = -jbar;  as[6] = as[5] - jbar * (t[6]-t[5]);
                ds[7] = dbar;   js[7] = 0.;     as[7] = 0.;
                ds[8] = 0.;     js[8] = 0.;     as[8] = 0.;
                ds[9] = -dbar;  js[9] = -jbar;  as[9] = as[8] - (1./2.) * jbar * (t[9]-t[8]);
                ds[10] = 0.;    js[10] = -jbar; as[10] = as[9] - jbar * (t[10]-t[9]);
                ds[11] = dbar;  js[11] = 0.;    as[11] = -abar;
                ds[12] = 0.;    js[12] = 0.;    as[12] = -abar;
                ds[13] = dbar;  js[13] = jbar;  as[13] = as[12] + (1./2.) * jbar * (t[13]-t[12]);
                ds[14] = 0.;    js[14] = jbar;  as[14] = as[13] + jbar * (t[14]-t[13]);
                ds[15] = -dbar; js[15] = 0.;    as[15] = 0.;
                vs[0] = 0.;
                vs[1] = vs[0] + (1./6.) * dbar * pow(t[1]-t[0],3);
                vs[2] = vs[1] + (1./2.) * jbar * pow(t[2]-t[1],2);
                vs[3] = vs[2] + (1./4.) * jbar * pow(t[3]-t[2],3);
                vs[4] = vs[3] + abar * (t[4]-t[3]);
                vs[5] = vs[4] + (1./6.) * dbar * pow(t[5]-t[4],3);
                vs[6] = vs[5] + (1./2.) * jbar * pow(t[6]-t[5],2);
                vs[7] = vs[6] + (1./6.) * dbar * pow(t[7]-t[6],3);
                vs[8] = vbar;
                vs[9] = vs[8] - (1./6.) * dbar * pow(t[9]-t[8],3);
                vs[10] = vs[9] - (1./2.) * jbar * pow(t[10]-t[9],2);
                vs[11] = vs[10] - (1./6.) * dbar * pow(t[11]-t[10],3);
                vs[12] = vs[11] - abar * (t[12]-t[11]);
                vs[13] = vs[12] - (1./6.) * dbar * pow(t[13]-t[12],3);
                vs[14] = vs[13] - (1./2.) * jbar * pow(t[14]-t[13],2);
                vs[15] = vs[14] - (1./6.) * dbar * pow(t[15]-t[14],3);
                xs[0] = 0.;
                xs[1] = xs[0] + (1./24.) * dbar * pow(t[1]-t[0],4);
                xs[2] = xs[1] + (1./6.) * jbar * pow(t[2]-t[1],3);
                xs[3] = xs[2] + (1./24.) * dbar * pow(t[3]-t[2],4);
                xs[4] = xs[3] + (1./2.) * abar * pow(t[4]-t[3],2);
                xs[5] = xs[4] + (1./24.) * dbar * pow(t[5]-t[4],4);
                xs[6] = xs[5] + (1./6.) * jbar * pow(t[6]-t[5],3);
                xs[7] = xs[6] + (1./24.) * dbar * pow(t[7]-t[6],4);
                xs[8] = xs[7] + vbar * (t[8]-t[7]);
                xs[9] = xs[8] + (1./24.) * dbar * pow(t[9]-t[8],4);
                xs[10] = xs[9] + (1./6.) * jbar * pow(t[10]-t[9],3);
                xs[11] = xs[10] + (1./24.) * dbar * pow(t[11]-t[10],4);
                xs[12] = xs[11] + (1./2.) * abar * pow(t[12]-t[11],2);
                xs[13] = xs[12] + (1./24.) * dbar * pow(t[13]-t[12],4);
                xs[14] = xs[13] + (1./6.) * jbar * pow(t[14]-t[13],3);
                xs[15] = xs[14] + (1./24.) * dbar * pow(t[15]-t[14],4);
                return std::make_tuple(xs, vs, as, js, ds);
            }

            template <typename ScalarPos, typename ScalarVel, typename ScalarAcc,
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
                auto [xp, vp, ap, jp, dp] = fourth_order_pivot(xbar, vbar, abar, jbar, dbar, tp);
                for (size_t i=0; i<std::size(time_points); i++) {
                    auto t = time_points[i];
                    assert(t > tp[0] && t < tp[15]);
                    if (t < tp[1]) {
                        ds[i] = dbar;
                        js[i] = 0. + dbar * t;
                        as[i] = 0. + (1./2.) * dbar * pow(t,2);
                        vs[i] = 0. + (1./6.) * dbar * pow(t,3);
                        as[i] = 0. + (1./24.) * dbar * pow(t,4);
                    } else if (t < tp[2]) {
                        ds[i] = 0;
                    }
                }
            }
        } // namespace axis1d
    } // namespace trajectory
} // namespace trjtools

#endif // TRJTOOLS_TRAJECTORY_AXIS1D_HPP