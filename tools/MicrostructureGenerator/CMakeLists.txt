cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Release)

PROJECT(microstructureGenerator CXX)
set(CMAKE_CXX_STANDARD 17)
set(MoDELibRootDir "/Users/giacomo/Documents/MODELIB2")

add_executable(microstructureGenerator microstructureGenerator.cpp)

#set_property(TARGET microstructureGenerator PROPERTY CXX_STANDARD 17)



include_directories(/usr/local/include)

include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/IO/)
include_directories(${MoDELibRootDir}/include/Math/)
include_directories(${MoDELibRootDir}/include/PolycrystallineMaterials/)
include_directories(${MoDELibRootDir}/include/FEM/)

include_directories(${MoDELibRootDir}/include/Geometry/)
add_subdirectory(${MoDELibRootDir}/src/Geometry Geometry)
target_link_libraries(microstructureGenerator MoDELibGeometry)

include_directories(${MoDELibRootDir}/include/Mesh/)
add_subdirectory(${MoDELibRootDir}/src/Mesh Mesh)
target_link_libraries(microstructureGenerator MoDELibMesh)

include_directories(${MoDELibRootDir}/include/Lattices/)
add_subdirectory(${MoDELibRootDir}/src/Lattices Lattices)
target_link_libraries(microstructureGenerator MoDELibLattices)

include_directories(${MoDELibRootDir}/include/DislocationDynamics/)
#add_subdirectory(${MoDELibRootDir}/src/DislocationDynamics DislocationDynamics)
#target_link_libraries(microstructureGenerator MoDELibDDD)





