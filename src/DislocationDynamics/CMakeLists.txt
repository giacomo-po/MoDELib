include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/LoopNetwork/)
include_directories(${MoDELibRootDir}/include/Quadrature/)
include_directories(${MoDELibRootDir}/include/ParticleInteraction/)
include_directories(${MoDELibRootDir}/include/SpaceDecomposition/)
include_directories(${MoDELibRootDir}/include/DiscreteCrackMechanics/)
include_directories(${MoDELibRootDir}/include/DislocationDynamics/)

add_library(MoDELibDDD 
	STATIC 
	DislocationLocalReference.cpp
		ConfinedDislocationObject.cpp
		DefectiveCrystal.cpp
		DislocationNetwork.cpp
		DislocationGlideSolver.cpp
		DislocationNode.cpp
		DDtimeIntegrator.cpp
		DislocationSegment.cpp
		DislocationLoop.cpp
		DislocationLoopLink.cpp
		DislocationLoopNode.cpp
		DislocationNetworkRemesh.cpp
		DislocationCrossSlip.cpp
		EshelbyInclusion.cpp
		ExternalLoadControllerBase.cpp
		SimplexBndNormal.cpp
			)

