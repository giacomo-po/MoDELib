include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/DislocationMobilities)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsIO)
add_library(MoDELibDDDmobilities 
	STATIC 
		DislocationMobilityBase.cpp
		DislocationMobilityBCC.cpp
		DislocationMobilityFCC.cpp
		DislocationMobilityHEXbasal.cpp
		DislocationMobilityHEXprismatic.cpp
		DislocationMobilityHEXpyramidal.cpp
			)

