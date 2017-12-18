
INCLUDE_DIRECTORIES(
	./
	../../Base
	../../GEO
	../../FileIO
	../../MathLib
	../../FEM
	../../MSH
	../../Qt/OpenSG
)

# Create executables
ADD_EXECUTABLE( testPointInTriangle testPointInTriangle.cpp )

TARGET_LINK_LIBRARIES( testPointInTriangle
	FEM
        MSH
       	GEO
        MathLib
       	FileIO
        Base
       	OGSProject
	${QT_LIBRARIES}
)
