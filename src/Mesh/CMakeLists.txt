include_directories(${MoDELibRootDir}/include/IO)
include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/Math/)
include_directories(${MoDELibRootDir}/include/Geometry/)
include_directories(${MoDELibRootDir}/include/Mesh/)
add_library(MoDELibMesh 
	STATIC 
		GmshReader.cpp
		MeshPlane.cpp
		PlanarMeshFace.cpp
		MeshBoundarySegment.cpp
		Simplex.cpp
		MeshRegion.cpp
		triangle.cpp
		SimplicialMesh.cpp
		TriangularMesh.cpp
			)

