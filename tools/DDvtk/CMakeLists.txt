cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Release)

PROJECT(DDvtk CXX)
set(CMAKE_CXX_STANDARD 17)
set(MoDELibRootDir "/home/ypachaur/Modelib_Development/MoDELib")

message("CMAKE_SYSTEM_NAME IS: ${CMAKE_SYSTEM_NAME}")
if(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
#	set(CMAKE_CXX_FLAGS "-fopenmp -O3 -Ofast")
elseif(${CMAKE_SYSTEM_NAME} MATCHES Linux)
	set(CMAKE_CXX_FLAGS "-fopenmp -march=native")
endif()

add_executable(DDvtk DDvtk.cxx)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

if(VTK_LIBRARIES)
  target_link_libraries(DDvtk ${VTK_LIBRARIES})
else()
  target_link_libraries(DDvtk vtkHybrid vtkWidgets)
endif()


include_directories(/home/ypachaur/ModelEssentials/Eigen)
#include_directories(../../header_symlinks)
include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/FEM/)
include_directories(${MoDELibRootDir}/include/Math/)
include_directories(${MoDELibRootDir}/include/PolycrystallineMaterials/)
include_directories(${MoDELibRootDir}/include/DislocationDynamics)
include_directories(${MoDELibRootDir}/include/DislocationDynamics/Visualization)
include_directories(${MoDELibRootDir}/include/IO)
include_directories(${MoDELibRootDir}/include/Geometry/)
include_directories(${MoDELibRootDir}/include/Mesh/)
include_directories(${MoDELibRootDir}/include/Lattices/)


add_subdirectory(${MoDELibRootDir}/src/Lattices Lattices)
target_link_libraries(DDvtk MoDELibLattices)

add_subdirectory(${MoDELibRootDir}/src/Mesh Mesh)
target_link_libraries(DDvtk MoDELibMesh)

add_subdirectory(${MoDELibRootDir}/src/Geometry Geometry)
target_link_libraries(DDvtk MoDELibGeometry)


