include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/LoopNetwork/)
include_directories(${MoDELibRootDir}/include/Quadrature/)
include_directories(${MoDELibRootDir}/include/ParticleInteraction/)
include_directories(${MoDELibRootDir}/include/SpaceDecomposition/)
include_directories(${MoDELibRootDir}/include/DiscreteCrackMechanics/)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsBase/)

target_sources(${PROJECT_NAME} PRIVATE 
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationLoopPatches.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationFieldBase.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/EshelbyInclusionBase.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/SphericalInclusion.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/PolyhedronInclusion.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/StraightDislocationSegment.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/StressStraight.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/VoigtTraits.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DDconfigFields.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationDynamicsBase.cpp
		)
