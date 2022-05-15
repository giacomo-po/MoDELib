include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/Lattices/)
include_directories(${MoDELibRootDir}/include/Math/)
include_directories(${MoDELibRootDir}/include/Geometry/)

#add_library(MoDELibLattices STATIC 
target_sources(${PROJECT_NAME} PRIVATE 
		${CMAKE_CURRENT_SOURCE_DIR}/CSL.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DSCL.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Lattice.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/LatticeDirection.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/LatticeGCD.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/LatticeLine.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/LatticePlane.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/LatticePlaneBase.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/LatticeTransitionMatrix.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/LatticeVector.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/LLL.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/PeriodicLatticeInterpolant.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/RationalLatticeDirection.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/RationalMatrix.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/RationalOrthogonalApproximation.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/ReciprocalLatticeDirection.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/ReciprocalLatticeVector.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/RLLL.cpp
			)

