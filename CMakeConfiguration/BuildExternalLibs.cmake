# Build all libraries in ExternalLibs/

if (NOT ExternalLibs_BUILT)

	# Built shapelib
	message(STATUS "Building shapelib...")
	
	# Windows
	if (WIN32)
		if (VS32)
		execute_process (
			COMMAND vcvars32.bat
			COMMAND nmake /f makefile.vc
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ExternalLibs/shapelib )
		
		else (VS32)
			if (VS64)
				execute_process (
					COMMAND vcvars32.bat
					COMMAND nmake /f makefile.vc
					WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ExternalLibs/shapelib )
			
			endif (VS64)
		endif (VS32)
	endif (WIN32)
	
	# Linux, Mac
	if ( UNIX )
 	    if ( APPLE )
 	        # OSX
			execute_process (
				COMMAND sh ./build_linux.sh
				WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ExternalLibs/shapelib )
 	    else ( APPLE )
 	        # Linux / Unix
			execute_process (
			COMMAND sh ./build_linux.sh
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ExternalLibs/shapelib )
 	    endif ( APPLE )
	endif ( UNIX )

	
	set (ExternalLibs_BUILT ON CACHE BOOL INTERNAL)
	
endif (NOT ExternalLibs_BUILT)
