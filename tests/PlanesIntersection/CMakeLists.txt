# on Windows use: cd build; cmake -G "NMake Makefiles" ..
cmake_minimum_required(VERSION 3.1.0)
project(test VERSION 1.0.0 LANGUAGES CXX)

set(EigenRootDir "/usr/local/include")
set(MoDELibRootDir "${CMAKE_CURRENT_SOURCE_DIR}/../..")

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#message("CMAKE_SYSTEM_NAME IS: ${CMAKE_SYSTEM_NAME}")
#if(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
#	set(CMAKE_CXX_COMPILER "/opt/local/bin/g++")
#	set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -Ofast")
#elseif(${CMAKE_SYSTEM_NAME} MATCHES Linux)
#	set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=native -O3 -Ofast")
#endif()
set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=native -Ofast")


add_executable(test main.cpp)

include_directories(${EigenRootDir})
include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/Math/)
include_directories(${MoDELibRootDir}/include/Geometry/)
include_directories(${MoDELibRootDir}/include/IO/)

add_subdirectory(${MoDELibRootDir}/src/Geometry Geometry)


