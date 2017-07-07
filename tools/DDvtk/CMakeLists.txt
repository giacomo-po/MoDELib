
cmake_minimum_required(VERSION 2.8)

PROJECT(DDvtk)
#PROJECT(DDvtk CXX) # forces g++, but linker fails


find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(DDvtk DDvtk)

include_directories(/Users/Giacomo/Documents/eigen)
include_directories(../../)


set(CMAKE_BUILD_TYPE Release)

set_property(TARGET DDvtk PROPERTY CXX_STANDARD 11)


if(VTK_LIBRARIES)
  target_link_libraries(DDvtk ${VTK_LIBRARIES})
else()
  target_link_libraries(DDvtk vtkHybrid vtkWidgets)
endif()
