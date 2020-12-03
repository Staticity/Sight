cmake_minimum_required(VERSION 3.12)

project(io)

list(APPEND sources
  src/video_device.cpp
  src/video_sources.cpp
)

add_library(${PROJECT_NAME} ${sources})

# Include OpenCV (for image reading & windows)
find_package(OpenCV REQUIRED)
target_include_directories(
  ${PROJECT_NAME}
  PRIVATE
  ${OpenCV_INCLUDE_DIRS})
target_link_libraries(
  ${PROJECT_NAME}
  PRIVATE
  ${OpenCV_LIBS}
)

# Include Eigen
find_package(Eigen3 3.3 REQUIRED)
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)

# Include yaml-cpp
find_package(yaml-cpp REQUIRED)
target_link_libraries(
    ${PROJECT_NAME}
    PRIVATE
    yaml-cpp)

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
    ${PROJECT_SOURCE_DIR}/inc)

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    calibration
    image)

