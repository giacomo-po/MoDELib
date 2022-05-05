include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/GlidePlanes)
include_directories(${MoDELibRootDir}/include/DislocationMobilities)
include_directories(${MoDELibRootDir}/include/PolycrystallineMaterials/)
include_directories(${MoDELibRootDir}/include/GlidePlanes/)
include_directories(${MoDELibRootDir}/include/FEM/)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsIO)
add_library(MoDELibDDDIO 
	STATIC 
	DefectiveCrystalParameters.cpp
	DDconfigIO.cpp
	
			)

