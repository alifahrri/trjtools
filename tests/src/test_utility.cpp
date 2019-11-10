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

TEST(traits, is_callable)
{
    auto callable_no_args = [](){};
    static_assert(
        trj::utility::traits::is_callable<decltype(callable_no_args)>::value
    );
    static_assert(
        !trj::utility::traits::is_callable<decltype(callable_no_args),double>::value
    );
    static_assert(
        !trj::utility::traits::is_callable<decltype(callable_no_args),double,double>::value
    );
    auto callable_single_args = [](auto){};
    static_assert(
        !trj::utility::traits::is_callable<decltype(callable_single_args)>::value
    );
    static_assert(
        trj::utility::traits::is_callable<decltype(callable_single_args),double>::value
    );
    static_assert(
        !trj::utility::traits::is_callable<decltype(callable_single_args),double,double>::value
    );
    auto callable_double_args = [](auto,auto){};
    static_assert(
        !trj::utility::traits::is_callable<decltype(callable_double_args)>::value
    );
    static_assert(
        !trj::utility::traits::is_callable<decltype(callable_double_args),double>::value
    );
    static_assert(
        trj::utility::traits::is_callable<decltype(callable_double_args),double,double>::value
    );
}

TEST(traits, is_iterable)
{
    static_assert(
        trj::utility::traits::is_iterable<std::vector<double>>::value
    );
    static_assert(
        trj::utility::traits::is_iterable<std::array<double,10>>::value
    );
    static_assert(
        !trj::utility::traits::is_iterable<double>::value
    );
}

TEST(traits, is_indexable)
{
    static_assert(
        trj::utility::traits::is_indexable<std::vector<double>>::value
    );
    static_assert(
        trj::utility::traits::is_indexable<std::array<double,10>>::value
    );
    static_assert(
        !trj::utility::traits::is_indexable<double>::value
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

TEST(mpl, copy_std_container_vector_different_size)
{
    using Iterable = std::vector<double>;
    using ScalarPos = float;
    constexpr size_t new_size = 6;
    static_assert(
        std::is_same_v<
            trj::utility::mpl::copy_std_container_t<Iterable,ScalarPos,new_size>,
            std::vector<float>
        >
    );
}

TEST(mpl, copy_std_container_array_different_size)
{
    using Iterable = std::array<double,3>;
    using ScalarPos = float;
    constexpr size_t new_size = 6;
    static_assert(
        std::is_same_v<
            trj::utility::mpl::copy_std_container_t<Iterable,ScalarPos,new_size>,
            std::array<float,6>
        >
    );
}

TEST(utility, interpolate)
{
    std::array<double,3> t{
        0.0, 1.0, 2.0
    };
    std::array<double,3> v{
        1.0, 2.0, 3.0
    };

    auto interpolated_fn = trj::utility::interpolate(t,v);
    auto p1 = interpolated_fn(0.5);
    auto p2 = interpolated_fn(1.5);
    EXPECT_NEAR(p1,1.5,1e-9);
    EXPECT_NEAR(p2,2.5,1e-9);
}

TEST(utility, sample)
{
    auto fn = [](double x){
        double x0 = 0.;
        double y0 = 1.;
        double m  = 1.;
        return y0 + m * (x-x0);
    };

    auto quad_fn = [](double x){
        return x*x;
    };

    std::array<double,3> tp_array{
        0.0, 1.0, 2.0
    };

    std::vector<double> tp_vector{
        0.0, 1.0, 2.0
    };

    auto sampled = trj::utility::sample(fn,tp_array);
    using sampled_type = decltype(sampled);
    static_assert(std::tuple_size<sampled_type>::value==3);
    EXPECT_NEAR(sampled[0],1.0,1e-9);
    EXPECT_NEAR(sampled[1],2.0,1e-9);
    EXPECT_NEAR(sampled[2],3.0,1e-9);

    auto sampled_vector = trj::utility::sample(fn,tp_vector);
    ASSERT_EQ(sampled_vector.size(),3);
    EXPECT_NEAR(sampled_vector[0],1.0,1e-9);
    EXPECT_NEAR(sampled_vector[1],2.0,1e-9);
    EXPECT_NEAR(sampled_vector[2],3.0,1e-9);
}

TEST(utility, sample_interpolate)
{
    auto fn = [](double x){
        double x0 = 0.;
        double y0 = 1.;
        double m  = 1.;
        return y0 + m * (x-x0);
    };

    auto quad_fn = [](double x){
        return x*x;
    };

    std::array<double,3> tp_array{
        0.0, 1.0, 2.0
    };

    std::vector<double> tp_vector{
        0.0, 1.0, 2.0
    };

    auto sampled = trj::utility::sample(fn,tp_array);
    auto interpolated = trj::utility::interpolate(sampled, tp_array);
    using sampled_type = decltype(sampled);
    using interpolated_type = decltype(interpolated);
    static_assert(
        trj::utility::traits::is_indexable<sampled_type>::value
    );
    static_assert(
        !trj::utility::traits::is_indexable<interpolated_type>::value
    );
    static_assert(
        !trj::utility::traits::is_callable<sampled_type,double>::value
    );
    static_assert(
        trj::utility::traits::is_callable<interpolated_type,double>::value
    );
}

TEST(utility, interpolate_integrate)
{
    std::array<double,3> x{
        0.0, 1.0, 2.0
    };

    std::array<double,3> y{
        1.0, 2.0, 3.0
    };

    auto interpolated_fn = trj::utility::interpolate(x,y);
    auto integrated_fn = trj::utility::integrate(interpolated_fn);

    using interpolated_type = decltype(interpolated_fn);
    using integrated_type = decltype(integrated_fn);
    static_assert(
        !trj::utility::traits::is_indexable<interpolated_type>::value
    );
    static_assert(
        trj::utility::traits::is_callable<interpolated_type,double>::value
    );
    static_assert(
        !trj::utility::traits::is_indexable<integrated_type>::value
    );
    static_assert(
        trj::utility::traits::is_callable<integrated_type,double>::value
    );
    EXPECT_NEAR(interpolated_fn(0.),1.,1e-9);
    EXPECT_NEAR(interpolated_fn(1.),2.,1e-9);
    EXPECT_NEAR(interpolated_fn(2.),3.,1e-9);
    EXPECT_NEAR(integrated_fn(1.),1.5,1e-9);
    EXPECT_NEAR(integrated_fn(2.),4.0,1e-9);
}

TEST(utility, split)
{
    double s = 1.0;
    size_t n = 20;
    auto splitted = trj::utility::split<std::vector<double>>(s,n);
    EXPECT_EQ(splitted.size(),21);
    EXPECT_NEAR(splitted.at(0),0.,1e-6);
    EXPECT_NEAR(splitted.at(19),0.95,1e-6);
}

TEST(utility, split_interpolate)
{
    double s = 1.0;
    size_t n = 20;
    auto splitted = trj::utility::split<std::vector<double>>(s,n);
    auto interpolated = trj::utility::interpolate(splitted, splitted);
    
    EXPECT_EQ(splitted.size(),21);
    EXPECT_NEAR(splitted.at(0),0.,1e-6);
    EXPECT_NEAR(splitted.at(19),0.95,1e-6);

    EXPECT_NEAR( interpolated(0.), 0., 1e-6);
    EXPECT_NEAR( interpolated(0.3), 0.3, 1e-6);
    EXPECT_NEAR( interpolated(0.5), 0.5, 1e-6);
    EXPECT_NEAR( interpolated(0.9), 0.9, 1e-6);
}

TEST(utility, split_interpolate_sample)
{
    double s = 1.0;
    size_t n = 20;
    auto splitted = trj::utility::split<std::vector<double>>(s,n);
    auto interpolated = trj::utility::interpolate(splitted, splitted);
    auto sampled = trj::utility::sample(interpolated, splitted);
    
    EXPECT_EQ(splitted.size(),21);
    EXPECT_EQ(sampled.size(),21);
    EXPECT_NEAR(splitted.at(0),0.,1e-6);
    EXPECT_NEAR(splitted.at(19),0.95,1e-6);
    EXPECT_NEAR(sampled.at(0),0.,1e-6);
    EXPECT_NEAR(sampled.at(19),0.95,1e-6);

    EXPECT_NEAR( interpolated(0.), 0., 1e-6);
    EXPECT_NEAR( interpolated(0.3), 0.3, 1e-6);
    EXPECT_NEAR( interpolated(0.5), 0.5, 1e-6);
    EXPECT_NEAR( interpolated(0.9), 0.9, 1e-6);
}