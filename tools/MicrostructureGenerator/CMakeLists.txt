cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Debug)

PROJECT(microstructureGenerator CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_PREFIX_PATH "/Applications/Qt/6.1.3/macos/lib/cmake/")
set(MoDELibRootDir "${CMAKE_CURRENT_SOURCE_DIR}/../..")
set(EigenRootDir "/usr/local/include")

include_directories(${EigenRootDir})


set(CMAKE_CXX_FLAGS_DEBUG "-O3")

add_executable(microstructureGenerator microstructureGenerator.cpp)



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
target_link_libraries(microstructureGenerator MoDELibDDD)


add_subdirectory(${MoDELibRootDir}/src/Lattices Lattices)
target_link_libraries(microstructureGenerator MoDELibLattices)

add_subdirectory(${MoDELibRootDir}/src/Mesh Mesh)
target_link_libraries(microstructureGenerator MoDELibMesh)

add_subdirectory(${MoDELibRootDir}/src/Geometry Geometry)
target_link_libraries(microstructureGenerator MoDELibGeometry)







