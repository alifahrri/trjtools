#include "trjtools/utility.hpp"
#include <gtest/gtest.h>
#include <vector>
#include <array>

namespace trj = trjtools;

TEST(traits, is_std_array_or_vector_uni_ref)
{
    static_assert(
        trj::utility::traits::is_std_array_or_vector<std::vector<double>&&>::value
    );
}

TEST(traits, is_std_array_or_vector_const)
{
    static_assert(
        trj::utility::traits::is_std_array_or_vector<const std::vector<double>>::value
    );
}

TEST(traits, is_std_array_or_vector)
{
    static_assert(
        trj::utility::traits::is_std_array_or_vector<std::vector<double>>::value
    );
}

TEST(mpl, copy_std_container_vector)
{
    using Iterable = std::vector<double>;
    using ScalarPos = float;
    static_assert(
        std::is_same_v<
            trj::utility::mpl::copy_std_container_t<Iterable,ScalarPos>,
            std::vector<float>
        >
    );
}

TEST(mpl, copy_std_container_array)
{
    using Iterable = std::array<double,3>;
    using ScalarPos = float;
    static_assert(
        std::is_same_v<
            trj::utility::mpl::copy_std_container_t<Iterable,ScalarPos>,
            std::array<float,3>
        >
    );
}