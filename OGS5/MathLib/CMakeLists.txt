# Source files
set( HEADERS
	Vector3.h
	Matrix.h
	MathTools.h
	AnalyticalGeometry.h
	max.h
	LinkedTriangle.h
	CubicSpline.h
	LinearInterpolation.h
	EarClippingTriangulation.h
	# directory LinAlg
	LinAlg/LinearSolver.h
	LinAlg/DirectLinearSolver.h
	LinAlg/DenseDirectLinearSolver.h
	LinAlg/GaussAlgorithm.h
	LinAlg/TriangularSolve.h
	LinAlg/IterativeLinearSolver.h
)

set( SOURCES
	MathTools.cpp
	AnalyticalGeometry.cpp
	LinkedTriangle.cpp
	CubicSpline.cpp
	LinearInterpolation.cpp
	EarClippingTriangulation.cpp
	# directory LinAlg
	LinAlg/GaussAlgorithm.cpp
	LinAlg/TriangularSolve.cpp
)

# Create the library
add_library( MathLib STATIC ${HEADERS} ${SOURCES} )

set_target_properties(MathLib PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	../Base
	../GEO
	.
	./LinAlg
)

target_link_libraries (
	MathLib
       Base
)

