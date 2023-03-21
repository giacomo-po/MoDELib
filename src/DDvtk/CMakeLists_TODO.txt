include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/DislocationMobilities)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsIO)
include_directories(${MoDELibRootDir}/include/DislocationMicrostructure/)

target_sources(${PROJECT_NAME} PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/MicrostructureGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/MicrostructureGeneratorBase.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PeriodicDipoleGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PeriodicLoopGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PrismaticLoopGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/InclusionsGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/IrradiationDefectsGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/VTKGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/VTKsegments.cpp
			)

