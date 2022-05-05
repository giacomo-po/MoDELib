include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/Geometry/)
include_directories(${MoDELibRootDir}/include/Mesh/)
include_directories(${MoDELibRootDir}/include/vtk/)
include_directories(${MoDELibRootDir}/include/DislocationDynamics/)
include_directories(${MoDELibRootDir}/include/PolycrystallineMaterials/)
include_directories(${MoDELibRootDir}/include/Lattices/)
include_directories(${MoDELibRootDir}/include/IO/)
include_directories(${MoDELibRootDir}/include/Math/)


set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

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


add_library(MoDELibvtk STATIC 
	${MoDELibRootDir}/include/vtk/SimplicialMeshActor.h
	SimplicialMeshActor.cpp
	${MoDELibRootDir}/include/vtk/DDqtMainWindow.h
	DDqtMainWindow.cpp
	${MoDELibRootDir}/include/vtk/DDqtVTKwidget.h
	DDqtVTKwidget.cpp
	${MoDELibRootDir}/include/vtk/DDconfigVtk.h
	DDconfigVtk.cpp
	${MoDELibRootDir}/include/vtk/NetworkNodeActor.h
	NetworkNodeActor.cpp
	${MoDELibRootDir}/include/vtk/NetworkLinkActor.h
	NetworkLinkActor.cpp
	${MoDELibRootDir}/include/vtk/InclusionActor.h
	InclusionActor.cpp
			)

target_link_libraries(MoDELibvtk MoDELibMesh MoDELibGeometry ${VTK_LIBRARIES} Qt6::Widgets Qt6::OpenGL Qt6::Core)