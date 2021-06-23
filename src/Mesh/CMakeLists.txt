include_directories(${MoDELibRootDir}/include/Utilities/)
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

