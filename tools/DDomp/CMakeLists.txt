cmake_minimum_required(VERSION 3.0)

PROJECT(DDomp CXX)

add_executable(DDomp main.cpp)

include_directories(/usr/local/include)
include_directories(../../header_symlinks)

set(CMAKE_BUILD_TYPE Debug)

set_property(TARGET DDomp PROPERTY CXX_STANDARD 17)


