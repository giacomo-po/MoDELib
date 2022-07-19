include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/Geometry/)
include_directories(${MoDELibRootDir}/include/PolycrystallineMaterials/)
target_sources(${PROJECT_NAME} PRIVATE 
		${CMAKE_CURRENT_SOURCE_DIR}/BCClattice.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/FCClattice.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Grain.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/GrainBoundary.cpp
#		${CMAKE_CURRENT_SOURCE_DIR}/HEXlattice.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Polycrystal.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/PolycrystallineMaterialBase.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/SingleCrystalBase.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/SlipSystem.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/GammaSurface.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/SecondPhase.cpp
		)

