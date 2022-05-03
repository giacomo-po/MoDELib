cmake_minimum_required(VERSION 3.1.0)

project(DDqt VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_PREFIX_PATH "/Applications/Qt/6.1.3/macos/lib/cmake/")

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

if(CMAKE_VERSION VERSION_LESS "3.7.0")
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif()

set(EigenRootDir "/usr/local/include")
set(MoDELibRootDir "${CMAKE_CURRENT_SOURCE_DIR}/../../")


include_directories(${EigenRootDir})
include_directories(${MoDELibRootDir}/include/DislocationDynamics)
include_directories(${MoDELibRootDir}/include/Mesh/)
include_directories(${MoDELibRootDir}/include/IO)
include_directories(${MoDELibRootDir}/include/Utilities)
include_directories(${MoDELibRootDir}/include/Math/)
include_directories(${MoDELibRootDir}/include/Geometry/)
include_directories(${MoDELibRootDir}/include/vtk/)
include_directories(${MoDELibRootDir}/include/PolycrystallineMaterials/)
include_directories(${MoDELibRootDir}/include/Lattices/)


add_subdirectory(${MoDELibRootDir}/src/Geometry Geometry)
add_subdirectory(${MoDELibRootDir}/src/Mesh Mesh)
add_subdirectory(${MoDELibRootDir}/src/Lattices Lattices)
add_subdirectory(${MoDELibRootDir}/src/vtk vtk)


find_package(Qt6 COMPONENTS Widgets OpenGL REQUIRED)
message (STATUS "Qt6_VERSION: ${Qt6_VERSION}")
message (STATUS "Qt6_LIBRARY_DIRS: ${Qt6_LIBRARY_DIRS}")
if (NOT Qt6_FOUND)
  message("Skipping DDvtk: ${Qt6_NOT_FOUND_MESSAGE}")
  return()
endif()

find_package(VTK REQUIRED)
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
message (STATUS "VTK_LIBRARY_DIRS: ${VTK_LIBRARY_DIRS}")
if (NOT VTK_FOUND)
  message("Skipping DDvtk: ${VTK_NOT_FOUND_MESSAGE}")
  return()
endif()


add_executable(DDqt DDqt.cxx)

target_link_libraries(DDqt MoDELibMesh MoDELibGeometry MoDELibvtk ${VTK_LIBRARIES} Qt6::Widgets Qt6::OpenGL Qt6::Core)
vtk_module_autoinit(TARGETS DDqt MODULES ${VTK_LIBRARIES} )
