
cmake_minimum_required(VERSION 2.8)

PROJECT(DDviewer)
#PROJECT(DDviewer CXX) # forces g++, but linker fails


find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(DDviewer DDviewer)

include_directories(/usr/local/include)
include_directories(../..)


set(CMAKE_BUILD_TYPE Release)

set_property(TARGET DDviewer PROPERTY CXX_STANDARD 11)


if(VTK_LIBRARIES)
  target_link_libraries(DDviewer ${VTK_LIBRARIES})
else()
  target_link_libraries(DDviewer vtkHybrid vtkWidgets)
endif()
