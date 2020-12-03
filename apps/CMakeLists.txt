cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(apps)

add_subdirectory(CalibrationTest)
add_subdirectory(OpenCVCalibration)
add_subdirectory(FindSO3ExpThreshold)
add_subdirectory(LucasKanadeVisual)
add_subdirectory(PlanarTracking)
add_subdirectory(Sandbox)
add_subdirectory(StereoDepth)
