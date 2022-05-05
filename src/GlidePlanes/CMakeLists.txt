include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/DislocationMobilities)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsIO)
add_library(MoDELibGlidePlanes 
	STATIC 
	GlidePlane.cpp
	PeriodicGlidePlaneFactory.cpp
	PeriodicGlidePlane.cpp
			)

