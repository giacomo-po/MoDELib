include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/DislocationMobilities)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsIO)
include_directories(${MoDELibRootDir}/include/DislocationMicrostructure/)

add_library(MoDELibDDDmicrostructure 
	STATIC 
	MicrostructureGenerator.cpp
			)

