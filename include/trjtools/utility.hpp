#ifndef TRJTOOLS_UTILITY_HPP
#define TRJTOOLS_UTILITY_HPP

#include "nmtools/integration/trapezoidal.hpp"

#include <iostream>
#include <array>
#include <tuple>
#include <vector>
#include <iterator>
#include <type_traits>

#define TRJTOOLS_USE_EXCEPTION
#ifdef TRJTOOLS_USE_EXCEPTION
#include <exception>
#define TRJTOOLS_ASSERT1(cond) if (!(cond)) throw std::runtime_error(#cond);
#define TRJTOOLS_ASSERT2(cond, msg) if (!(cond)) throw std::runtime_error(msg);
#else
#include <cassert>
#define TRJTOOLS_ASSERT1(cond) assert(cond);
#define TRJTOOLS_ASSERT2(cond, msg) assert(cond);
#endif

#define GET_MACRO(_1,_2,NAME,...) NAME
#define TRJTOOLS_ASSERT(...) GET_MACRO(__VA_ARGS__, TRJTOOLS_ASSERT2, TRJTOOLS_ASSERT1)(__VA_ARGS__)

namespace trjtools {
    namespace utility {
        namespace traits {

            template <typename T, typename = void>
            struct is_std_array : std::false_type {};
            template <typename T>
            struct is_std_array<T,std::enable_if_t<
                std::is_same_v<std::array<typename T::value_type, std::tuple_size<T>::value>, T>
            > > : std::true_type {};

            template <typename T, typename = void>
            struct is_std_vector : std::false_type {};
            template <typename T>
            struct is_std_vector<T,std::enable_if_t<
                std::is_same_v<std::vector<typename T::value_type>, T>
            > > : std::true_type {};


            template <typename T, typename ...Args>
            struct is_callable {
            private:
                template <typename F>
                static constexpr auto test(int) 
                    -> decltype(std::declval<F>()(std::declval<Args>()...), bool())
                { return true; }
                template <typename F>
                static constexpr auto test(char) 
                { return false; }
            public:
                constexpr static bool value = test<T>(int());
            };

            /* TODO : move (?) */
            using std::begin;
            using std::end;

            template <typename T>
            struct is_iterable {
            private:
                template <typename It>
                static constexpr auto test(int)
                    -> decltype(begin(std::declval<It>()), end(std::declval<It>()), bool())
                { return true; }
                template <typename It>
                static constexpr auto test(char)
                { return false; }
            public:
                constexpr static bool value = test<T>(int());
            };

            template <typename T, typename = void>
            struct is_indexable : std::false_type {};
            template <typename T>
            struct is_indexable<T, std::void_t<decltype(std::declval<T>()[size_t{}])> > : std::true_type {};

            template <typename T, typename = void>
            struct has_push_back : std::false_type {};
            template <typename T>
            struct has_push_back<T, 
                std::void_t<decltype(std::declval<T>().push_back(std::declval<typename T::value_type>()))> > 
            : std::true_type {};

            template <typename T, typename = void>
            struct has_resize : std::false_type {};
            template <typename T>
            struct has_resize<T, 
                std::void_t<decltype(std::declval<T>().resize(std::declval<typename T::size_type>()))> > 
            : std::true_type {};

            template <typename T, typename = void>
            struct is_resizeable : std::false_type {};
            template <typename T>
            struct is_resizeable<T, 
                std::void_t<decltype(std::declval<T>().resize(std::declval<typename T::size_type>()))> > 
            : std::true_type {};

            template <typename T>
            using is_std_array_or_vector = std::disjunction<is_std_array<std::decay_t<T>>,is_std_vector<std::decay_t<T>>>;

        } // namespace traits

        namespace mpl {
            template <typename T, typename V, size_t, typename=void>
            struct copy_std_container {
                using type = std::enable_if_t<
                    traits::is_std_array_or_vector<T>::value, std::vector<std::decay_t<V>>
                >;
            };
            template <typename T, typename V>
            struct copy_std_container<T,V,0,std::enable_if_t<traits::is_std_array<T>::value>> {
                using type = std::enable_if_t<
                    traits::is_std_array_or_vector<T>::value, std::array<std::decay_t<V>,std::tuple_size<T>::value>
                >;
            };
            template <typename T, typename V, size_t n>
            struct copy_std_container<T,V,n,std::enable_if_t<traits::is_std_array<T>::value>> {
                using type = std::enable_if_t<
                    traits::is_std_array_or_vector<T>::value, std::array<std::decay_t<V>,n>
                >;
            };
            template <typename T, typename V, size_t n=0>
            using copy_std_container_t = typename copy_std_container<T,V,n>::type;
        } // namespace mpl

        /* TODO : move to nmtools */
        /**
         * @brief 
         * 
         * @param y1 
         * @param y0 
         * @param x1 
         * @param x0 
         * @param x 
         * @return auto 
         */
        constexpr auto linear_interpolation(auto y1, auto y0, auto x1, auto x0, auto x)
        {
            TRJTOOLS_ASSERT(x1>x0);
            auto m = (y1-y0) / (x1-x0);
            return y0 + (x-x0) * m;
        }

        /**
         * @brief overloaded version of interpolation
         * 
         * @tparam IterableX 
         * @tparam IterableY 
         * @param x 
         * @param y 
         * @return auto 
         */
        template <typename IterableX, typename IterableY>
        constexpr auto interpolate(const IterableX &x, const IterableY &y)
        {
            static_assert(traits::is_indexable<IterableX>::value && traits::is_indexable<IterableY>::value,
                "Iterable class should be indexable, e.g. i[0] is valid!"
            );
            if constexpr(traits::has_resize<IterableX>::value && 
                traits::has_resize<IterableY>::value) {
                    TRJTOOLS_ASSERT(std::size(x)==std::size(y));
            } else {
                static_assert(std::tuple_size<IterableX>::value==std::tuple_size<IterableY>::value, 
                    "expect tuple_size(x)==tuple_size(y)"
                );
            }
            return [=](const typename IterableX::value_type &t){
                TRJTOOLS_ASSERT((t>=x.front() && t<=x.back()), 
                    std::to_string(t)
                );
                auto x0 = x[0]; auto x1 = x[1]; 
                auto y0 = y[0]; auto y1 = y[1]; 
                for (size_t i=1; i<std::size(x); i++) {
                    auto _x0 = x[i-1]; auto _y0 = y[i-1];
                    auto _x1 = x[i]; auto _y1 = y[i];
                    if ( t<_x1 && t>=_x0 ) {
                        x0 = _x0; x1 = _x1;
                        y0 = _y0; y1 = _y1;
                        break;
                    }
                }
                return linear_interpolation(
                    y1, y0, x1, x0, t
                );
            };
        }

        /**
         * @brief convert callable to iterable by sampling fn from tp
         * 
         * @tparam Callable 
         * @tparam Iterable 
         * @param fn 
         * @param tp 
         * @return auto 
         */
        template <typename Callable, typename Iterable>
        constexpr auto sample(const Callable &fn, const Iterable &tp)
        {
            using result_value_type = decltype(fn(tp[0]));
            mpl::copy_std_container_t<Iterable,result_value_type> result{};
            if constexpr (traits::has_resize<Iterable>::value)
                result.resize(tp.size());
            for (size_t i=0;  i<std::size(tp); i++)
                result[i] = fn(tp[i]);
            return result;
        }

        /**
         * @brief convert callable to callable integral
         * 
         * @tparam integral_segment=100 
         * @tparam Callable 
         * @param fn 
         * @return auto 
         */
        template <size_t integral_segment=100, typename Callable>
        constexpr auto integrate(const Callable &fn)
        {
            return [fn](const auto &t) {
                return numeric::integration::trapezoid<integral_segment>(fn, 0., t);
            };
        } // integrate

        /**
         * @brief split s to n segment, returning (n+1) points
         * 
         * @tparam Container target container type
         * @tparam Scalar 
         * @param s 
         * @param n 
         * @return constexpr auto container with size (n+1)
         */
        template <typename Container, typename Scalar>
        constexpr auto split(const Scalar &s, size_t n)
        {
            static_assert(
                std::is_same_v<typename Container::value_type, Scalar>
            );
            Container c{};
            if constexpr (traits::is_resizeable<Container>::value)
                c.resize(n+1);
            TRJTOOLS_ASSERT(std::size(c)==(n+1), 
                std::to_string(std::size(c))
            );
            auto ds = s / n;
            for (size_t i=0; i<=n; i++)
                c[i] = i * ds;
            return c;
        }

        /**
         * @brief overloaded version of split, container type defaulted to std::vector
         * 
         * @tparam Scalar 
         * @param s 
         * @param n 
         * @return auto std::vector<Scalar> of size (n+1)
         * @note this function can't be constexpr
         */
        template <typename Scalar>
        auto split(const Scalar &s, size_t n)
        {
            return split<std::vector<Scalar>>(s,n);
        }
    } // namespace utility
} // namespace trjtools

#endif // TRJTOOLS_UTILITY_HPP