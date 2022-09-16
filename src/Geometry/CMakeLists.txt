include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/Geometry/)
target_sources(${PROJECT_NAME} PRIVATE 
		${CMAKE_CURRENT_SOURCE_DIR}/FiniteLineSegment.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/LineLineIntersection.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/LineSimplexIntersection.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/PlanarPolygon.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Plane.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/PlaneLineIntersection.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/PlanePlaneIntersection.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/PlaneSegmentIntersection.cpp
#		${CMAKE_CURRENT_SOURCE_DIR}/PolygonClipper.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/SegmentSegmentDistance.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/PlanesIntersection.cpp
		)

