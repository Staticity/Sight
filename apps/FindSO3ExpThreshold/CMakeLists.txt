cmake_minimum_required(VERSION 3.10)

project(FindSO3ExpThreshold)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(
    ${PROJECT_NAME}
    PRIVATE
    lie)