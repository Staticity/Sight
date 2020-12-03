cmake_minimum_required(VERSION 3.10)

project(utils)

list(APPEND sources
  src/string_helpers.cpp
)

add_library(${PROJECT_NAME} ${sources})

# Include Eigen
find_package(Eigen3 3.3 REQUIRED)
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
    ${PROJECT_SOURCE_DIR}/inc)

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    linear)

