add_executable(
    ${PROJECT_NAME}-example1
    src/example1.cpp
)

add_executable(
    ${PROJECT_NAME}-example2
    src/example2.cpp
)

set_target_properties(
    ${PROJECT_NAME}-example1
    PROPERTIES
    CXX_STANDARD 17
)

set_target_properties(
    ${PROJECT_NAME}-example2
    PROPERTIES
    CXX_STANDARD 17
)


target_link_libraries(
    ${PROJECT_NAME}-example2
    pthread
)