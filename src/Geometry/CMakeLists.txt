include_directories(${MoDELibRootDir}/include/Utilities/)
add_library(MoDELibGeometry
	STATIC
		FiniteLineSegment.cpp
		LineLineIntersection.cpp
		LineSimplexIntersection.cpp
		PlanarPolygon.cpp
		Plane.cpp
		PlaneLineIntersection.cpp
		PlanePlaneIntersection.cpp
		PlaneSegmentIntersection.cpp
		PolygonClipper.cpp
		SegmentSegmentDistance.cpp
			)

