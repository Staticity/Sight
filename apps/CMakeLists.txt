cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(apps)

add_subdirectory(FindSO3ExpThreshold)
add_subdirectory(LucasKanadeVisual)
add_subdirectory(Sandbox)
