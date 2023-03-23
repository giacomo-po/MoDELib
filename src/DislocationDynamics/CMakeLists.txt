include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/LoopNetwork/)
include_directories(${MoDELibRootDir}/include/Quadrature/)
include_directories(${MoDELibRootDir}/include/ParticleInteraction/)
include_directories(${MoDELibRootDir}/include/SpaceDecomposition/)
include_directories(${MoDELibRootDir}/include/DiscreteCrackMechanics/)
include_directories(${MoDELibRootDir}/include/DislocationDynamics/)

target_sources(${PROJECT_NAME} PRIVATE 
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationLocalReference.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/ConfinedDislocationObject.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DefectiveCrystal.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationNetwork.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationGlideSolver.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationNode.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DDtimeIntegrator.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationSegment.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationLoop.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationLoopLink.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationLoopNode.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationNetworkRemesh.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationCrossSlip.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/CrossSlipModels.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/ExternalLoadControllerBase.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/SimplexBndNormal.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationQuadraturePoint.cpp
		)
