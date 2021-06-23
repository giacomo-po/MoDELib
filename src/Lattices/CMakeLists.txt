include_directories(${MoDELibRootDir}/include/Utilities/)
add_library(MoDELibLattices 
	STATIC 
		CSL.cpp
		DSCL.cpp
		Lattice.cpp
		LatticeDirection.cpp
		LatticeGCD.cpp
		LatticeLine.cpp
		LatticePlane.cpp
		LatticePlaneBase.cpp
		LatticeTransitionMatrix.cpp
		LatticeVector.cpp
		LLL.cpp
		PeriodicLatticeInterpolant.cpp
		RationalLatticeDirection.cpp
			)

