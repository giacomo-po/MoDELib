include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/Geometry/)
include_directories(${MoDELibRootDir}/include/Mesh/)
include_directories(${MoDELibRootDir}/include/vtk/)
include_directories(${MoDELibRootDir}/include/DislocationDynamics/)
include_directories(${MoDELibRootDir}/include/PolycrystallineMaterials/)
include_directories(${MoDELibRootDir}/include/Lattices/)

find_package(VTK REQUIRED)
find_package(Qt6 COMPONENTS Widgets OpenGL Core REQUIRED)


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
			)

target_link_libraries(MoDELibvtk ${VTK_LIBRARIES} Qt6::Widgets Qt6::OpenGL Qt6::Core)