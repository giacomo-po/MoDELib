cmake_minimum_required(VERSION 3.0)

PROJECT(DDvtk CXX)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(DDvtk DDvtk.cxx)

include_directories(/usr/local/include)
include_directories(../../header_symlinks)

set(CMAKE_BUILD_TYPE Release)

set_property(TARGET DDvtk PROPERTY CXX_STANDARD 14)

if(VTK_LIBRARIES)
  target_link_libraries(DDvtk ${VTK_LIBRARIES})
else()
  target_link_libraries(DDvtk vtkHybrid vtkWidgets)
endif()
