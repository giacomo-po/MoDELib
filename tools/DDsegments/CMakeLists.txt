# on Windows use: cd build; cmake -G "NMake Makefiles" ..
cmake_minimum_required(VERSION 3.1.0)
project(DDsegments VERSION 1.0.0 LANGUAGES CXX)

set(EigenRootDir "/usr/local/include")
set(MoDELibRootDir "${CMAKE_CURRENT_SOURCE_DIR}/../..")

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message("CMAKE_SYSTEM_NAME IS: ${CMAKE_SYSTEM_NAME}")
if(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
	set(CMAKE_CXX_COMPILER "/opt/local/bin/g++-mp-11")
#	set(CMAKE_CXX_FLAGS "-Wall -Wextra -fopenmp -O3 -Ofast")
elseif(${CMAKE_SYSTEM_NAME} MATCHES Linux)
#	set(CMAKE_CXX_FLAGS "-Wall -Wextra -fopenmp -march=native -O3 -Ofast")
endif()
set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=native -fopenmp -Ofast")



add_executable(DDsegments DDsegments.cpp)

include_directories(${EigenRootDir})
include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/IO/)
include_directories(${MoDELibRootDir}/include/Math/)
include_directories(${MoDELibRootDir}/include/Geometry/)
include_directories(${MoDELibRootDir}/include/Mesh/)
include_directories(${MoDELibRootDir}/include/Lattices/)
include_directories(${MoDELibRootDir}/include/FEM/)
include_directories(${MoDELibRootDir}/include/Quadrature/)
include_directories(${MoDELibRootDir}/include/GlidePlanes/)
include_directories(${MoDELibRootDir}/include/PolycrystallineMaterials/)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsIO/)
include_directories(${MoDELibRootDir}/include/DislocationMobilities/)
include_directories(${MoDELibRootDir}/include/DislocationMicrostructure/)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsBase/)


add_subdirectory(${MoDELibRootDir}/src/Geometry Geometry)
add_subdirectory(${MoDELibRootDir}/src/Lattices Lattices)
add_subdirectory(${MoDELibRootDir}/src/Mesh Mesh)
add_subdirectory(${MoDELibRootDir}/src/GlidePlanes GlidePlanes)
add_subdirectory(${MoDELibRootDir}/src/DislocationMobilities DDDmobilities)
add_subdirectory(${MoDELibRootDir}/src/DislocationDynamicsIO DDDIO)
add_subdirectory(${MoDELibRootDir}/src/PolycrystallineMaterials PolyMaterials)

add_subdirectory(${MoDELibRootDir}/src/DislocationMicrostructure Microstructure)
