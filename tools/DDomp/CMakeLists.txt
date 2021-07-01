cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Release)

PROJECT(DDomp CXX)
set(CMAKE_CXX_STANDARD 17)
set(MoDELibRootDir "/Users/giacomo/Documents/MODELIB2")

message("CMAKE_SYSTEM_NAME IS: ${CMAKE_SYSTEM_NAME}")
if(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
#	set(CMAKE_CXX_FLAGS "-fopenmp -O3 -Ofast")
elseif(${CMAKE_SYSTEM_NAME} MATCHES Linux)
	set(CMAKE_CXX_FLAGS "-fopenmp -march=native")
endif()


#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(DDomp DDomp.cpp)


include_directories(/usr/local/include)

include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/IO/)
include_directories(${MoDELibRootDir}/include/Math/)
include_directories(${MoDELibRootDir}/include/PolycrystallineMaterials/)
include_directories(${MoDELibRootDir}/include/FEM/)
include_directories(${MoDELibRootDir}/include/Quadrature/)
include_directories(${MoDELibRootDir}/include/LoopNetwork/)
include_directories(${MoDELibRootDir}/include/ParticleInteraction/)
include_directories(${MoDELibRootDir}/include/SpaceDecomposition/)
include_directories(${MoDELibRootDir}/include/DiscreteCrackMechanics/)


include_directories(${MoDELibRootDir}/include/Geometry/)
add_subdirectory(${MoDELibRootDir}/src/Geometry Geometry)
target_link_libraries(DDomp MoDELibGeometry)

include_directories(${MoDELibRootDir}/include/Mesh/)
add_subdirectory(${MoDELibRootDir}/src/Mesh Mesh)
target_link_libraries(DDomp MoDELibMesh)

include_directories(${MoDELibRootDir}/include/Lattices/)
add_subdirectory(${MoDELibRootDir}/src/Lattices Lattices)
target_link_libraries(DDomp MoDELibLattices)

include_directories(${MoDELibRootDir}/include/DislocationDynamics/)
#add_subdirectory(${MoDELibRootDir}/src/DislocationDynamics DislocationDynamics)
#target_link_libraries(microstructureGenerator MoDELibDDD)





