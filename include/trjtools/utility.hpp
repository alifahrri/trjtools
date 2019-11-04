#ifndef TRJTOOLS_UTILITY_HPP
#define TRJTOOLS_UTILITY_HPP

#include <array>
#include <vector>
#include <tuple>
#include <type_traits>

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

            template <typename T>
            using is_std_array_or_vector = std::disjunction<is_std_array<std::decay_t<T>>,is_std_vector<std::decay_t<T>>>;
        } // namespace traits

        namespace mpl {
            template <typename T, typename V, typename=void>
            struct copy_std_container {
                using type = std::enable_if_t<
                    traits::is_std_array_or_vector<T>::value, std::vector<std::decay_t<V>>
                >;
            };
            template <typename T, typename V>
            struct copy_std_container<T,V,std::enable_if_t<traits::is_std_array<T>::value>> {
                using type = std::enable_if_t<
                    traits::is_std_array_or_vector<T>::value, std::array<std::decay_t<V>,std::tuple_size<T>::value>
                >;
            };
            template <typename T, typename V>
            using copy_std_container_t = typename copy_std_container<T,V>::type;
        } // namespace mpl
    } // namespace utility
} // namespace trjtools

#endif // TRJTOOLS_UTILITY_HPP