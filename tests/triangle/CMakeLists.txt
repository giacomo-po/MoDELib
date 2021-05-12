cmake_minimum_required(VERSION 3.0)

project(tri CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall")
set(CMAKE_CXX_FLAGS_DUBUG "-O3 -Wall")

include_directories(../../header_symlinks)
include_directories(/usr/local/include)
#add_definitions(-DTRILIBRARY)
#add_definitions(-DANSI_DECLARATORS)

add_executable(tri main.cpp ${src_files})
