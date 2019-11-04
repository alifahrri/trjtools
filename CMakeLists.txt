cmake_minimum_required(VERSION 3.0)
project(trjtools)

include_directories(include)

option(BUILD_TESTS "Build test programs" ON)
option(BUILD_EXAMPLES "Build example programs" OFF)

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif(BUILD_EXAMPLES)

if(BUILD_TESTS)
    add_subdirectory(tests)
endif(BUILD_TESTS)