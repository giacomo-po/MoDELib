include_directories(${MoDELibRootDir}/include/Math/)
include_directories(${MoDELibRootDir}/include/Utilities/)
add_library(MoDELibLattices 
	STATIC  BiCrystal.cpp
		Lattice.cpp
		LatticeDirection.cpp
		LatticeCore.cpp
		LatticeVector.cpp
		LLL.cpp
		RationalLatticeDirection.cpp
		RationalMatrix.cpp
		ReciprocalLatticeDirection.cpp
		ReciprocalLatticeVector.cpp
		RLLL.cpp
		PeriodicLatticeInterpolant.cpp
		LatticePlane.cpp
		LatticePlaneBase.cpp
			)

