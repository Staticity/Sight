cmake_minimum_required(VERSION 3.10)

project(features)

add_library(${PROJECT_NAME} INTERFACE)

# Include Eigen
find_package(Eigen3 3.3 REQUIRED)
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)

target_include_directories(
    ${PROJECT_NAME}
    INTERFACE
    ${PROJECT_SOURCE_DIR}/inc)

target_link_libraries(
    ${PROJECT_NAME}
    INTERFACE
    linear
    image
    utils)