include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/Lattices)
include_directories(${MoDELibRootDir}/include/GlidePlanes)
include_directories(${MoDELibRootDir}/include/DislocationMobilities)
include_directories(${MoDELibRootDir}/include/PolycrystallineMaterials/)
include_directories(${MoDELibRootDir}/include/GlidePlanes/)
include_directories(${MoDELibRootDir}/include/FEM/)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsIO)
target_sources(${PROJECT_NAME} PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/DefectiveCrystalParameters.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/DDconfigIO.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/DDauxIO.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/DDtraitsIO.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/DDbaseIO.cpp
	)

