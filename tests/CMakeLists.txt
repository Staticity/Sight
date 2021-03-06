cmake_minimum_required(VERSION 3.12)

project(tests)

add_executable(tests
    test.cpp
    calibration/camera_test.cpp
    geometric/homography_test.cpp
    # geometric/polar_rectify_test.cpp # Tested functions moved to anonymous namespace
    geometric/triangulation_test.cpp
    features/fast_test.cpp
    features/orb_test.cpp
    image/image_gen_test.cpp
    image/image_ops_test.cpp
    lie/so3_test.cpp
    lie/se3_test.cpp
    lie/vec_test.cpp
    lie/exp_utils_test.cpp
    optimization/camera_test.cpp
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
  io
  lie
  linear
  optimization
  utils
)