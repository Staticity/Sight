cmake_minimum_required(VERSION 3.10)

project(calibration)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(
    ${PROJECT_NAME}
    INTERFACE
    ${PROJECT_SOURCE_DIR}/inc)

# target_link_libraries(
#     ${PROJECT_NAME}
#     INTERFACE
#     linear
#     image
#     utils)