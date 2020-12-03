cmake_minimum_required(VERSION 3.10)

project(algorithm)

add_library(
    ${PROJECT_NAME}
    src/union_find.cpp
)

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
    ${PROJECT_SOURCE_DIR}/inc)

# target_link_libraries(
#     ${PROJECT_NAME}
#     INTERFACE
#     linear
#     image
#     utils)