
INCLUDE_DIRECTORIES(
	./
	../../Base
	../../GEO
	../../FileIO
	../../MathLib
	../../FEM
	../../MSH
)

# Create executables
ADD_EXECUTABLE( ExtractMeshNodes mainExtractMeshNodes.cpp ExtractMeshNodes.h ExtractMeshNodes.cpp )
ADD_EXECUTABLE( ModifyMeshProperties mainModifyMeshProperties.cpp ModifyMeshProperties.h ModifyMeshProperties.cpp )
ADD_EXECUTABLE( filterMeshNodes filterMeshNodes.cpp )
ADD_EXECUTABLE( extractSurface mainExtractSurface.cpp ExtractSurface.h ExtractSurface.cpp )

TARGET_LINK_LIBRARIES( ExtractMeshNodes
	GEO
	MSH
	FileIO
	FEM
	OGSProject
)

TARGET_LINK_LIBRARIES( ModifyMeshProperties
	GEO
	MSH
	MathLib
	FileIO
	FEM
	OGSProject
)

TARGET_LINK_LIBRARIES( filterMeshNodes
	Base
)

TARGET_LINK_LIBRARIES( extractSurface
	GEO
	MSH
	FileIO
	FEM
	OGSProject
	${QT_LIBRARIES}
)


IF (OGS_PACKAGING)
	INSTALL ( TARGETS ExtractMeshNodes RUNTIME DESTINATION bin COMPONENT
ExtractMeshNodes )
ENDIF (OGS_PACKAGING)
