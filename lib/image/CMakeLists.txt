cmake_minimum_required(VERSION 3.10)

project(image)

list(APPEND sources
  src/image.cpp
  src/image_ops.cpp
)

add_library(${PROJECT_NAME} ${sources})

# Include OpenCV
find_package(OpenCV REQUIRED)

target_include_directories(
  ${PROJECT_NAME}
  PUBLIC
  ${OpenCV_INCLUDE_DIRS})

target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC
  ${OpenCV_LIBS}
)

# Include Eigen
find_package(Eigen3 3.3 REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
    ${PROJECT_SOURCE_DIR}/inc)

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    utils)

