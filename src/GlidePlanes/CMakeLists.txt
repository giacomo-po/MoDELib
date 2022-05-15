include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/DislocationMobilities)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsIO)
include_directories(${MoDELibRootDir}/include/GlidePlanes)
add_library(MoDELibGlidePlanes 
	STATIC 
	GlidePlane.cpp
	PeriodicGlidePlaneFactory.cpp
	PeriodicGlidePlane.cpp
			)

