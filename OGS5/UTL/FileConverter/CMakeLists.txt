
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

## Shape to GLI file converter ##
IF (Shapelib_FOUND)
	
	INCLUDE_DIRECTORIES (${Shapelib_INCLUDE_DIR})
	
	# Create executables
	ADD_EXECUTABLE( ConvertSHPToGLI ConvertSHPToGLI.cpp )

	TARGET_LINK_LIBRARIES( ConvertSHPToGLI
		FEM
	        MSH
        	#GEO
	        MathLib
        	FileIO
	        Base
        	OGSProject
		${Shapelib_LIBRARIES}
		${QT_LIBRARIES}
	)

	IF (OGS_PACKAGING)
		INSTALL ( TARGETS ConvertSHPToGLI RUNTIME DESTINATION bin COMPONENT ConvertSHPToGLI )
	ENDIF (OGS_PACKAGING)
	
ENDIF (Shapelib_FOUND)

## VTK to OpenSG file converter ##
IF (OGS_USE_OPENSG)
	
	FIND_PACKAGE( VTK REQUIRED )
	INCLUDE( ${VTK_USE_FILE} )
	
	ADD_EXECUTABLE (ConvertVtkToOsg ConvertVtkToOsg.cpp)
	
	USE_OPENSG ( ConvertVtkToOsg )
	INCLUDE_DIRECTORIES ( ../OpenSG )
	TARGET_LINK_LIBRARIES ( ConvertVtkToOsg OgsOpenSG vtkRendering )
	
	IF (OGS_PACKAGING)
		INSTALL ( TARGETS ConvertVtkToOsg RUNTIME DESTINATION bin COMPONENT ConvertVtkToOsg )
	ENDIF (OGS_PACKAGING)
	
ENDIF (OGS_USE_OPENSG)
