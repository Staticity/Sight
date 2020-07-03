cmake_minimum_required(VERSION 3.10)

project(tests)

add_executable(tests
    test.cpp
    geometric/homography_test.cpp
    geometric/triangulation_test.cpp
    features/fast_test.cpp
    features/orb_test.cpp
    image/image_gen_test.cpp
    image/image_ops_test.cpp
    lie/so3_test.cpp
    lie/se3_test.cpp
    lie/vec_test.cpp
    lie/exp_utils_test.cpp
)

# Include Eigen
find_package(Eigen3 3.3 REQUIRED)
target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC
  Eigen3::Eigen
)

#include Catch
find_package(Catch2 REQUIRED)
target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    Catch2::Catch2)

# Add custom libraries
target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC
  features
  geometric
  image
  lie
  linear
  utils
)