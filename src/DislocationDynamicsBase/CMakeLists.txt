include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/LoopNetwork/)
include_directories(${MoDELibRootDir}/include/Quadrature/)
include_directories(${MoDELibRootDir}/include/ParticleInteraction/)
include_directories(${MoDELibRootDir}/include/SpaceDecomposition/)
include_directories(${MoDELibRootDir}/include/DiscreteCrackMechanics/)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsBase/)

target_sources(${PROJECT_NAME} PRIVATE 
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationLoopPatches.cpp
#		${CMAKE_CURRENT_SOURCE_DIR}/DislocationQuadraturePoint.cpp
		)
