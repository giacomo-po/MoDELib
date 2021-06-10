cmake_minimum_required(VERSION 3.0)

PROJECT(test CXX)
set(MoDELibRootDir "/Users/giacomo/Documents/MODELIB2/tests/templates")

include_directories(${MoDELibRootDir}/include/Module1/)
add_subdirectory(${MoDELibRootDir}/src/Module1 Module1)

add_executable(test main.cpp)
target_link_libraries(test Module1)

set(CMAKE_BUILD_TYPE Release)
set_property(TARGET test PROPERTY CXX_STANDARD 20)

