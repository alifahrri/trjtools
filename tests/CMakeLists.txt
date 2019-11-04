find_package(GTest REQUIRED)
add_executable(
    ${PROJECT_NAME}-test
    src/test_axis1d.cpp
    src/test_utility.cpp
)
target_link_libraries(${PROJECT_NAME}-test ${GTEST_BOTH_LIBRARIES} pthread)
target_compile_options(${PROJECT_NAME}-test PRIVATE -g)
add_test(all-${PROJECT_NAME}-test ${PROJECT_NAME}-test)

set_target_properties(
    ${PROJECT_NAME}-test
    PROPERTIES
    CXX_STANDARD 17
)