# Source files
SET( SOURCES
	GMSHInterface.cpp
	GMSInterface.cpp
	GocadInterface.cpp
	OGSIOVer4.cpp
	PetrelInterface.cpp	
	StationIO.cpp
	#FEMIO.cpp
	FEMIO/GeoIO.cpp
)

# Header files
SET( HEADERS
	GMSHInterface.h
	GMSInterface.h
	GocadInterface.h
	OGSIOVer4.h
	PetrelInterface.h
	StationIO.h
	#FEMIO.h
	FEMIO/GeoIO.h
	MathIO/CRSIO.h
)

IF (OGS_USE_NETCDF)
	FIND_PACKAGE( VTK REQUIRED )
	INCLUDE( ${VTK_USE_FILE} )
	SET (HEADERS NetCDFInterface.h ${HEADERS})
	SET (SOURCES NetCDFInterface.cpp ${SOURCES})
ENDIF (OGS_USE_NETCDF)

IF (Shapelib_FOUND)
	SET( HEADERS ${HEADERS} SHPInterface.h )
	SET (SOURCES ${SOURCES} SHPInterface.cpp )
	INCLUDE_DIRECTORIES (${Shapelib_INCLUDE_DIR})
ENDIF (Shapelib_FOUND)

IF (QT4_FOUND)
	SET( HEADERS ${HEADERS} XMLInterface.h )
	SET (SOURCES ${SOURCES} XMLInterface.cpp )
ENDIF (QT4_FOUND)

# Create the library
ADD_LIBRARY( FileIO STATIC
	${SOURCES}
	${HEADERS}
)

INCLUDE_DIRECTORIES(
	.
	./FEMIO
	./MathIO
	../Base
	../FEM
	../GEO
	../MathLib
	../MSH
)

TARGET_LINK_LIBRARIES( FileIO
#	Base
#	MathLib
	GEO
	FEM
)

IF (Shapelib_FOUND)
	TARGET_LINK_LIBRARIES( FileIO ${Shapelib_LIBRARIES} )
ENDIF (Shapelib_FOUND)

IF(OGS_PACKAGING)
	# Install xml schema files
	FILE (GLOB XML_SCHEMA_FILES . "*.xs?")
	INSTALL (FILES ${XML_SCHEMA_FILES} DESTINATION bin)
ENDIF(OGS_PACKAGING)
