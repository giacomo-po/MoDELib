include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/DislocationMobilities)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsIO)
include_directories(${MoDELibRootDir}/include/GlidePlanes)

target_sources(${PROJECT_NAME} PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/GlidePlane.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GlidePlaneFactory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PeriodicGlidePlaneFactory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PeriodicGlidePlane.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GlidePlaneNoise.cpp
			)

