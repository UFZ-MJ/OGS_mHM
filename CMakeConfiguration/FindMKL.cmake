
if (NOT MKL_FOUND)
	
	#Extract MKL
	FIND_PATH (MKL_INCLUDE_DIR_FOUND mkl.h ${PROJECT_SOURCE_DIR}/../Libs/MKL/include)
	IF (NOT MKL_INCLUDE_DIR_FOUND)
		FIND_PATH (MKL_DIR_FOUND mkl-include.tgz ${PROJECT_SOURCE_DIR}/../Libs/MKL)
		IF (MKL_DIR_FOUND)
			MESSAGE (STATUS "Uncompressing MKL...")
			EXECUTE_PROCESS (COMMAND tar xvzf mkl-include.tgz WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../Libs/MKL/ OUTPUT_QUIET)
			IF (HAVE_64_BIT)
				EXECUTE_PROCESS (COMMAND tar xvzf mkl-64.tgz WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../Libs/MKL/ OUTPUT_QUIET)
			ELSE (HAVE_64_BIT)
				EXECUTE_PROCESS (COMMAND tar xvzf mkl-32.tgz WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../Libs/MKL/ OUTPUT_QUIET)
			ENDIF (HAVE_64_BIT)
		ELSE (MKL_DIR_FOUND)
			MESSAGE (STATUS "MKL archives in ../Libs/ not found")
		ENDIF (MKL_DIR_FOUND)
	ENDIF (NOT MKL_INCLUDE_DIR_FOUND)

	include(LibFindMacros)
	
	find_path( MKL_INCLUDE_DIR
		NAMES mkl.h
		PATHS ${CMAKE_SOURCE_DIR}/../Libs/MKL/include)

	if ( UNIX )
		# Tell if the unix system is on 64-bit base
		if(CMAKE_SIZEOF_VOID_P MATCHES "8")
		    if(CMAKE_C_COMPILER MATCHES "icc")
		        find_library(MKL_LIBRARIES
				    NAMES mkl_solver_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread
				    PATHS ${CMAKE_SOURCE_DIR}/../Libs/MKL/64 )
		    else(CMAKE_C_COMPILER MATCHES "icc")
			    find_library(MKL_LIBRARIES
				    NAMES mkl_solver_lp64 mkl_intel_lp64 mkl_gnu_thread mkl_core
				    PATHS ${CMAKE_SOURCE_DIR}/../Libs/MKL/64 )	
		    endif(CMAKE_C_COMPILER MATCHES "icc") 
		else (CMAKE_SIZEOF_VOID_P MATCHES "8")
			find_library(MKL_LIBRARIES
				NAMES mkl_solver mkl_intel mkl_gnu_thread mkl_core
				PATHS ${CMAKE_SOURCE_DIR}/../Libs/MKL/32 )	
		endif (CMAKE_SIZEOF_VOID_P MATCHES "8")	
	endif ( UNIX )

	# Set the include dir variables and the libraries and let libfind_process do the rest.
	# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
	if (NOT MKL_LIBRARIES STREQUAL "MKL_LIBRARIES-NOTFOUND" AND NOT MKL_INCLUDE_DIR STREQUAL "MKL_INCLUDE_DIR-NOTFOUND")
		set(MKL_PROCESS_INCLUDES MKL_INCLUDE_DIR)
		set(MKL_PROCESS_LIBS MKL_LIBRARIES)
		libfind_process(MKL)
	else (NOT MKL_LIBRARIES STREQUAL "MKL_LIBRARIES-NOTFOUND" AND NOT MKL_INCLUDE_DIR STREQUAL "MKL_INCLUDE_DIR-NOTFOUND")
		message (STATUS "Warning: MKL not found!")
	endif (NOT MKL_LIBRARIES STREQUAL "MKL_LIBRARIES-NOTFOUND" AND NOT MKL_INCLUDE_DIR STREQUAL "MKL_INCLUDE_DIR-NOTFOUND")
	
endif (NOT MKL_FOUND)
