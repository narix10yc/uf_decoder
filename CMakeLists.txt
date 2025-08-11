cmake_minimum_required(VERSION 3.10)

project(uf_decoder LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
include_directories(include)
add_compile_options(-Wall)

add_subdirectory(src)
add_subdirectory(tools)