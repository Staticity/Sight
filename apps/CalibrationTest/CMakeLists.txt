cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(CalibrationTest)

add_executable(${PROJECT_NAME} main.cpp)

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

find_package(yaml-cpp REQUIRED)
target_link_libraries(
  ${PROJECT_NAME}
  PRIVATE
  yaml-cpp)

# Add custom libraries
target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC
  calibration
  features
  geometric
  image
  io
  lie
  linear
  optimization
  utils
)
