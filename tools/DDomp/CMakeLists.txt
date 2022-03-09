cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Debug)

PROJECT(DDomp CXX)
set(CMAKE_CXX_STANDARD 17)
set(MoDELibRootDir "/home/ypachaur/Collaboration_Model/Model_CMAKE")

message("CMAKE_SYSTEM_NAME IS: ${CMAKE_SYSTEM_NAME}")
if(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
#	set(CMAKE_CXX_FLAGS "-fopenmp -O3 -Ofast")
elseif(${CMAKE_SYSTEM_NAME} MATCHES Linux)
	set(CMAKE_CXX_FLAGS "-fopenmp -march=native")
endif()


set(CMAKE_CXX_FLAGS_DEBUG "-O3")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#set(CMAKE_CXX_FLAGS_RELEASE "-H")


add_executable(DDomp DDomp.cpp)


include_directories(/home/ypachaur/ModelEssentials/Eigen)

include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/MPI/)
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
include_directories(${MoDELibRootDir}/include/Mesh/)
include_directories(${MoDELibRootDir}/include/Lattices/)
include_directories(${MoDELibRootDir}/include/DislocationDynamics/)


add_subdirectory(${MoDELibRootDir}/src/DislocationDynamics DislocationDynamics)
target_link_libraries(DDomp MoDELibDDD)

add_subdirectory(${MoDELibRootDir}/src/Lattices Lattices)
target_link_libraries(DDomp MoDELibLattices)

add_subdirectory(${MoDELibRootDir}/src/Mesh Mesh)
target_link_libraries(DDomp MoDELibMesh)

add_subdirectory(${MoDELibRootDir}/src/Geometry Geometry)
target_link_libraries(DDomp MoDELibGeometry)



