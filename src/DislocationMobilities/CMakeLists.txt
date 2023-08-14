include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/DislocationMobilities)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsIO)
target_sources(${PROJECT_NAME} PRIVATE  
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationMobilityBase.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationMobilityBCC.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationMobilityFCC.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationMobilityHEXbasal.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationMobilityHEXprismatic.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationMobilityHEXpyramidal.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationMobilitySelector.cpp
		)

