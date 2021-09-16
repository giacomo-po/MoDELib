include_directories(${MoDELibRootDir}/include/Utilities/)
add_library(MoDELibDDD 
	STATIC 
	DislocationLocalReference.cpp
	DefectiveCrystalParameters.cpp
		ConfinedDislocationObject.cpp
		DislocationNetwork.cpp
		DislocationNode.cpp
		DDtimeIntegrator.cpp
		DislocationSegment.cpp
		DislocationLoop.cpp
		DislocationLoopLink.cpp
		DislocationLoopNode.cpp
		DislocationNetworkRemesh.cpp
		DislocationCrossSlip.cpp
			)

