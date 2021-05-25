cmake_minimum_required(VERSION 3.0)

PROJECT(microstructureGenerator CXX)

add_executable(microstructureGenerator microstructureGenerator.cpp)

include_directories(/usr/local/include)
include_directories(../../header_symlinks)

set(CMAKE_BUILD_TYPE Debug)

set_property(TARGET microstructureGenerator PROPERTY CXX_STANDARD 17)


