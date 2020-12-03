cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set("OpenCV_DIR" "C:\\Users\\Jaime\\source\\repos\\opencv\\bin")

project(Sight)

add_subdirectory(apps)
add_subdirectory(lib)
add_subdirectory(tests)
