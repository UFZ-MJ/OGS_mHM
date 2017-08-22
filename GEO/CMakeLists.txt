# Source files
set( SOURCES
	Color.cpp
	geo_dom.cpp
	geo_lib.cpp
	geo_lin.cpp
	geo_mathlib.cpp
	geo_ply.cpp
	geo_pnt.cpp
	geo_sfc.cpp
	geo_vol.cpp
	Point.cpp
	Polyline.cpp
	PolylineSmoother.cpp
	Polygon.cpp
	ProjectData.cpp
	SimplePolygonTree.cpp
	Surface.cpp
	AxisAlignedBoundingBox.cpp
	Station.cpp
	GEOObjects.cpp
	PointVec.cpp
	GeoType.cpp
	BruteForceClosestPair.cpp
	DivideAndConquerClosestPair.cpp
)

# Header files
set( HEADERS
	Color.h
	geo_dom.h
	geo_lib.h
	geo_lin.h
	geo_mathlib.h
	geo_ply.h
	geo_pnt.h
	geo_sfc.h	
	geo_vol.h
	AxisAlignedBoundingBox.h
	TemplatePoint.h
	Point.h
	PointWithID.h
	Polyline.h
	PolylineSmoother.h
	Polygon.h
	ProjectData.h
	SimplePolygonTree.h
	Triangle.h
	Surface.h
	PointVec.h
	TemplateVec.h
	PolylineVec.h
	SurfaceVec.h
	PropertyBounds.h
	Station.h	
	GEOObjects.h
	GeoObject.h
	GeoType.h
	QuadTree.h
	ClosestPair.h
	BruteForceClosestPair.h
	DivideAndConquerClosestPair.h
)

# Create the library
add_library( GEO STATIC
	${SOURCES}
	${HEADERS}
)


include_directories(
	.
	../Base
	../FEM
	../FileIO
	../MathLib
	../MSH
)


target_link_libraries (
	GEO
	Base
	MathLib
)

