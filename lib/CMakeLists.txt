cmake_minimum_required(VERSION 3.10)

project(SightLib)

# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_subdirectory(calibration)
add_subdirectory(estimation)
add_subdirectory(features)
add_subdirectory(geometric)
add_subdirectory(image)
add_subdirectory(io)
add_subdirectory(lie)
add_subdirectory(linear)
add_subdirectory(optimization)
add_subdirectory(utils)
