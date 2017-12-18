#### OGS-5 Project ####
## Programmed by Lars Bilke
## Modified by:
##            WW 20.04.2010. 

# Specify minimum CMake version
cmake_minimum_required(VERSION 2.6)	# or set to 2.4 for older systems (like JUGENE...)

# Project name
project( OGS-5 )

option ( CMAKE_CMD_ARGS "Build ogs5 storing the CMake command line args (hint: use cmake.ogs.sh" OFF )

if( CMAKE_CMD_ARGS )
  message( STATUS "CMake args: ${CMAKE_CMD_ARGS}")
endif( CMAKE_CMD_ARGS )

# OGS version information. Adjust these if you release a new version.
SET (OGS_VERSION_MAJOR 5)
SET (OGS_VERSION_MINOR 1)
SET (OGS_VERSION_PATCH 00)
SET (OGS_RELEASE_PERSONS "KR/TF/LB")
SET (OGS_VERSION "${OGS_VERSION_MAJOR}.${OGS_VERSION_MINOR}.${OGS_VERSION_PATCH}(${OGS_RELEASE_PERSONS})")
# should we use ISO dates: 2010-12-31 ?
SET (OGS_DATE "03.03.2011")

# Adds useful macros and variables
INCLUDE( CMakeConfiguration/Macros.cmake )

# Suppress waring on add_subdirectory(dir) where dir contains no CMakeLists.txt
IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.7)
	CMAKE_POLICY(SET CMP0014 OLD)
ENDIF ()

IF (WIN32)
	## For Visual Studio compiler
	IF (MSVC)
		ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
		ADD_DEFINITIONS(-D_CRT_NONSTDC_NO_WARNINGS)
		ADD_DEFINITIONS(-D_CRT_XNONSTDC_NO_WARNINGS)
		SET ( CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE TYPE INTERNAL FORCE )
		# Sets warning level 3 and ignores some warnings
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /wd4290")
	ELSE (MSVC)
#FOR CYGWIN.  25.02.2010. WW
		MESSAGE (STATUS "Might be GCC under cygwin.")
#		MESSAGE (FATAL_ERROR "Aborting: On Windows only the Visual Studio compiler is supported!")
	ENDIF (MSVC)
ENDIF (WIN32)

### For GNU C/CXX. WW
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
	IF( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
		MESSAGE(STATUS "Set GCC release flags")
		SET(CMAKE_CXX_FLAGS "-O3 -DNDEBUG")
	ENDIF()
	# -g
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wall -Wextra -fno-nonansi-builtins")
        ADD_DEFINITIONS(
		-DGCC -std=c++11
	)
	OPTION (OGS_PROFILE "Enables compiling with flags set for profiling with gprof." OFF)
	IF (OGS_PROFILE)
		SET(PROFILE_FLAGS "-pg -fno-omit-frame-pointer -O2 -DNDEBUG -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROFILE_FLAGS}")
		SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PROFILE_FLAGS}ctest")
	ENDIF (OGS_PROFILE)
ENDIF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)


# Set build directories
SET( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )
SET( LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )
IF (MSVC)
	SET(OGS_EXECUTABLE ${EXECUTABLE_OUTPUT_PATH}/release/ogs)
ELSE (MSVC)
	SET(OGS_EXECUTABLE ${EXECUTABLE_OUTPUT_PATH}/ogs)
ENDIF (MSVC)

# Find precompiled libraries (for BRNS GEMS LIS)
FIND_PATH (OGS_PRECOMPILED_LIBS_DIR_FOUND GEMS3_rl.lib ${PROJECT_SOURCE_DIR}/../Libs/precompiled)
IF (OGS_PRECOMPILED_LIBS_DIR_FOUND)
	INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/../Libs/precompiled)
	LINK_DIRECTORIES (${PROJECT_SOURCE_DIR}/../Libs/precompiled)
ELSE (OGS_PRECOMPILED_LIBS_DIR_FOUND)
	IF (WIN32)
		IF (OGS_FEM_BRNS OR OGS_FEM_GEMS OR OGS_FEM_CHEMAPP)
			MESSAGE (FATAL_ERROR "Precompiled libraries not found! Make sure to also checked out the trunk/Libs directory beneath your sources directory.")
		ENDIF (OGS_FEM_BRNS OR OGS_FEM_GEMS OR OGS_FEM_CHEMAPP)
	ELSE (WIN32)
		IF (OGS_FEM_LIS)
			MESSAGE (FATAL_ERROR "Precompiled libraries not found! Make sure to also checked out the trunk/Libs directory beneath your sources directory.")
		ENDIF (OGS_FEM_LIS)	
	ENDIF (WIN32)
ENDIF (OGS_PRECOMPILED_LIBS_DIR_FOUND)

FIND_PATH (LIB_DIR_FOUND GEMS3_rl.lib ${PROJECT_SOURCE_DIR}/LIB)
IF (LIB_DIR_FOUND)
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/LIB )
	LINK_DIRECTORIES( ${PROJECT_SOURCE_DIR}/LIB )
ENDIF (LIB_DIR_FOUND)

# Set additional CMake modules path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeConfiguration")


###-----WW------------------------------------------------------
#OPTION Parallel
OPTION(PARALLEL_USE_MPI OFF)
OPTION(PARALLEL_USE_OPENMP OFF)
IF(_OPENMP)
	SET(PARALLEL_USE_OPENMP ON)
ENDIF(_OPENMP)
IF(PARALLEL_USE_OPENMP)
	SET(_OPENMP ON)
ENDIF(PARALLEL_USE_OPENMP)

# Configurations
OPTION(OGS_FEM "Default configuration" OFF)
OPTION(OGS_FEM_SP "Sparse matrix configuration" OFF)
OPTION(OGS_FEM_GEMS "GEMS configuration" OFF)
OPTION(OGS_FEM_PQC "Phreeqc configuration" OFF)
OPTION(OGS_FEM_BRNS "BRNS configuration" OFF)
OPTION(OGS_FEM_CHEMAPP "CHEMAPP configuration" OFF)
OPTION(OGS_USE_QT "User interface configuration" OFF)
OPTION(OGS_FEM_LIS "Library of Iterative Solvers for Linear Systems configuration" OFF)

IF (UNIX) # Only supported on Linux
	OPTION(OGS_FEM_MKL "Math kernel library configuration" OFF)
	OPTION(OGS_FEM_MPI "Message passing interface configuration" OFF)
ENDIF (UNIX)

# Set default configuration when no other config is given
IF (NOT OGS_FEM AND NOT OGS_FEM_SP AND NOT OGS_FEM_GEMS AND NOT OGS_FEM_PQC AND NOT OGS_FEM_BRNS AND NOT OGS_FEM_CHEMAPP AND NOT OGS_FEM_LIS AND NOT OGS_FEM_MKL AND NOT OGS_FEM_MPI AND NOT OGS_USE_QT)
	MESSAGE (STATUS "No configuration specified. Default confuguration is used.")
	SET (OGS_FEM ON)
ENDIF (NOT OGS_FEM AND NOT OGS_FEM_SP AND NOT OGS_FEM_GEMS AND NOT OGS_FEM_PQC AND NOT OGS_FEM_BRNS AND NOT OGS_FEM_CHEMAPP AND NOT OGS_FEM_LIS AND NOT OGS_FEM_MKL AND NOT OGS_FEM_MPI AND NOT OGS_USE_QT)

# Additional options
OPTION(OGS_PACKAGING "Creating installers / packages" OFF)
IF (OGS_PACKAGING)
	OPTION (OGS_PACKAGING_ZIP "Do you want to package as zip?" OFF)
ENDIF (OGS_PACKAGING)
OPTION(OGS_PYTHON "Enable python scripting interface" OFF)
OPTION(OGS_BUILD_UTILITIES "Build additional utility programs?" OFF)
OPTION(OGS_CMAKE_DEBUG "Show additional cmake debug information" OFF)
OPTION(OGS_ONE_BIG_GTEST "Should all gtests be collected to one ctest" ON)
OPTION(OGS_USE_NETCDF "NetCDF interface" OFF)
OPTION(OGS_USE_OPENSG "VTK to OpenSG converter" OFF)
OPTION(OGS_BUILD_INFO "Should build information be generated" ON)

# Check if a valid OGS configuration is given (see Macros.cmake)
CHECK_CONFIG()

# Check for number of processors
CHECK_PROCESSOR_COUNT()

###LIBRARY
FIND_PACKAGE (PythonInterp)
FIND_PACKAGE( Shapelib )

IF (OGS_PYTHON)
	FIND_PACKAGE (PythonLibs 2.5 REQUIRED)
ENDIF (OGS_PYTHON)

# Find Qt4 library
IF(NOT OGS_DONT_USE_QT)
	FIND_PACKAGE( Qt4 4.5)
ENDIF(NOT OGS_DONT_USE_QT)

IF ( QT4_FOUND )
	# OPTION(OGS_GUI OFF )
	# this is needed to correctly link the qt libraries through target_link_libraries
	# By default only QtCore and QtGui modules are enabled
	# other modules must be enabled like this:
	SET( QT_USE_QTOPENGL TRUE )
	SET( QT_USE_QTSQL TRUE )
	SET( QT_USE_QTTEST TRUE )
	SET( QT_USE_QTXML TRUE )
	IF (QT_QTXMLPATTERNS_FOUND)
		set( QT_USE_QTXMLPATTERNS TRUE )
	ENDIF (QT_QTXMLPATTERNS_FOUND)
	INCLUDE( ${QT_USE_FILE} )
	IF (OGS_USE_QT)
		SET(OGS_USE_NETCDF ON)
	ENDIF (OGS_USE_QT)
ENDIF (QT4_FOUND )

IF (OGS_BUILD_INFO)
	find_package( Subversion )
	if(Subversion_FOUND)
	  #Subversion_WC_INFO(${CMAKE_SOURCE_DIR} "svn") , not a working copy
	  execute_process(
	    COMMAND "svnversion" "--no-newline" ${CMAKE_SOURCE_DIR}
	    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	    OUTPUT_VARIABLE svn_revision
	    ERROR_VARIABLE svn_error
	    )
	  set( SVN_REVISION "${svn_revision}" )
	endif( Subversion_FOUND )

	if( NOT MSVC )
	  # see if there is a file that will work with find_package
	  find_file( GIT FindGit.cmake CMakeConfiguration/ )
	  if( GIT STREQUAL GIT-NOTFOUND )
	    message( WARNING "there is no FindGit.cmake file, so I am downloading one that is licensed under the GPL, copyright Szilard Pall (pszilard@cbr.su.se).  Do not commit it to the repository." )
	    file( DOWNLOAD http://repo.or.cz/w/gromacs.git/blob_plain/f81efa276ee7e1c41197d9d8b56b91c7764e0641:/cmake/FindGit.cmake ${PROJECT_SOURCE_DIR}/CMakeConfiguration/FindGit.cmake TIMEOUT 10 )
	  endif( GIT STREQUAL GIT-NOTFOUND )
	    execute_process(
	      COMMAND "date" "+%Y-%m-%d %H:%M:%S"
	      OUTPUT_VARIABLE build_timestamp
	      OUTPUT_STRIP_TRAILING_WHITESPACE
	      )
	    if( build_timestamp )
	       set( BUILD_TIMESTAMP "${build_timestamp}" )
	       message( STATUS "Build timestamp: ${BUILD_TIMESTAMP}" )
	    endif()
	  find_package( Git )
	  if(Git_FOUND)
	    execute_process(
	      COMMAND "git" "describe" "--always"
	      COMMAND "xargs" "-Ixx" "git" "name-rev" "xx"
	      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	      OUTPUT_VARIABLE git_commit_info
	      OUTPUT_STRIP_TRAILING_WHITESPACE
	      ERROR_VARIABLE not_git_repo
	      )
	    if( git_commit_info )
	      #message( STATUS "Git executable version: " ${Git_VERSION} )
	      message( STATUS "Git commit: " ${git_commit_info} )
	      set( GIT_COMMIT_INFO "${git_commit_info}" )
	    else( git_commit_info )
	      if( not_git_repo )
		message( STATUS "Git not versioning the source" )
	      endif( not_git_repo )
	    endif( git_commit_info )
	  endif( Git_FOUND )
	endif( NOT MSVC )
ENDIF (OGS_BUILD_INFO)

# This is for Configure.h which is generated later
INCLUDE_DIRECTORIES( ${PROJECT_BINARY_DIR}/Base )

IF(OGS_FEM)
   MESSAGE (STATUS  "Configuring for FEM command line" )	
ENDIF(OGS_FEM)

IF(OGS_FEM_SP)
	MESSAGE (STATUS  "Configuring for FEM command line with new sparse matrix" )	
	SET( PRJ_EXT "-SP" )
	INCLUDE(configs/sparse.cmake)
ENDIF(OGS_FEM_SP)

IF(OGS_FEM_MPI)
	MESSAGE (STATUS  "Configuring for FEM command line with MPI" )
	SET( PRJ_EXT "-MPI" )
	SET(PARALLEL_USE_MPI ON)	
	INCLUDE(configs/sparse.cmake)
ENDIF(OGS_FEM_MPI)

IF(OGS_FEM_GEMS)
	MESSAGE (STATUS  "Configuring for FEM command line with GEMS" )
	SET( PRJ_EXT "-GEMS" )
	INCLUDE(configs/gems.cmake)
ENDIF(OGS_FEM_GEMS)

IF(OGS_FEM_PQC)
	MESSAGE (STATUS  "Configuring for FEM command line with PQC" )	
	SET( PRJ_EXT "-PQC" )
	INCLUDE(configs/PQC.cmake)
    ADD_SUBDIRECTORY ( PQC )
ENDIF(OGS_FEM_PQC)

IF(OGS_FEM_BRNS)
	MESSAGE (STATUS  "Configuring for FEM command line with BRNS" )
	SET( PRJ_EXT "-BRNS" )
	INCLUDE(configs/BRNS.cmake)
ENDIF(OGS_FEM_BRNS)

IF(OGS_FEM_LIS)
	MESSAGE (STATUS  "Configuring for FEM command line with LIS" )
	SET( PRJ_EXT "-LIS" )
	SET(PARALLEL_USE_OPENMP ON)
	INCLUDE(configs/vc_lis.cmake)
ENDIF(OGS_FEM_LIS)

IF(OGS_FEM_MKL)
	MESSAGE (STATUS  "Configuring for FEM command line with MKL" )
	SET( PRJ_EXT "-MKL" )
	SET(PARALLEL_USE_OPENMP ON)
	INCLUDE(configs/pardiso.cmake)
ENDIF(OGS_FEM_MKL)

IF(OGS_FEM_CHEMAPP)
	MESSAGE (STATUS  "Configuring for FEM command line with CHEMAPP" )
	SET( PRJ_EXT "-CHEMAPP" )
	INCLUDE(configs/chemapp.cmake)
ENDIF(OGS_FEM_CHEMAPP)

INCLUDE(FEM.cmake )

IF(OGS_USE_QT)
   MESSAGE (STATUS  "Configuring for OGS with GUI" )	
   SET(OGS_FEM ON )
   #INCLUDE(FEM.cmake )
   INCLUDE(GUI.cmake )
ENDIF(OGS_USE_QT)

# Create Configure.h
CONFIGURE_FILE (Base/Configure.h.in ${PROJECT_BINARY_DIR}/Base/Configure.h)
CONFIGURE_FILE (Base/BuildInfo.h.in ${PROJECT_BINARY_DIR}/Base/BuildInfo.h)

IF (MSVC)
	ADD_SUBDIRECTORY (CMakeConfiguration)
ENDIF (MSVC)

IF(CMAKE_MINOR_VERSION GREATER 6)
	IF (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" AND NOT OGS_FEM_GEMS )
		#SET (GTEST_MSVC_SEARCH MT)
		SET (GTEST_ROOT ${PROJECT_SOURCE_DIR}/../Libs/gtest)
	ENDIF (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" AND NOT OGS_FEM_GEMS )
	FIND_PACKAGE (GTest)
ENDIF(CMAKE_MINOR_VERSION GREATER 6)

IF(DEFINED BENCHMARK_DIR)
	FIND_PATH (BENCHMARK_DIR_FOUND copy.py ${BENCHMARK_DIR})
ELSE()
	FIND_PATH (BENCHMARK_DIR_FOUND copy.py ${PROJECT_SOURCE_DIR}/../benchmarks)
ENDIF()

IF (BENCHMARK_DIR_FOUND OR GTEST_FOUND)
	ENABLE_TESTING()
ENDIF (BENCHMARK_DIR_FOUND OR GTEST_FOUND)


IF(DEFINED EXAMPLEDATA_DIR)
	FIND_PATH (EXAMPLEDATA_DIR_FOUND points.gli ${EXAMPLEDATA_DIR})
ELSE()
	FIND_PATH (EXAMPLEDATA_DIR_FOUND points.gli ${PROJECT_SOURCE_DIR}/../ExampleData)
ENDIF()


IF (EXAMPLEDATA_DIR_FOUND)
	INSTALL (DIRECTORY ${EXAMPLEDATA_DIR_FOUND} DESTINATION . PATTERN .svn EXCLUDE)
ENDIF (EXAMPLEDATA_DIR_FOUND)

IF (GTEST_FOUND)
	include_directories(${GTEST_INCLUDE_DIRS})
	add_subdirectory(tests)
ENDIF (GTEST_FOUND)

IF (OGS_PACKAGING)
	# Add packaging configuration
	MESSAGE (STATUS "Enabling packaging")
	INCLUDE( CMakeConfiguration/Pack.cmake)
ENDIF (OGS_PACKAGING)

IF (OGS_CMAKE_DEBUG)
	INCLUDE (CMakeConfiguration/PrintVariables.cmake)
ENDIF (OGS_CMAKE_DEBUG)	

IF (BENCHMARK_DIR_FOUND)
	STRING (COMPARE EQUAL "$ENV{HOSTNAME}" "dev2.intern.ufz.de" HOST_IS_DEV2)
	MESSAGE (STATUS "Benchmarking through CTest enabled")
	IF (${HOST_IS_DEV2})
		MESSAGE (STATUS "On dev2: file comparing enabled")
	ENDIF (${HOST_IS_DEV2})
	
	FILE (GLOB BENCHMARK_CONFIGS "${PROJECT_SOURCE_DIR}/../benchmarks/*.cmake")
	FOREACH (BENCHMARK_CONFIG ${BENCHMARK_CONFIGS})
		INCLUDE ("${BENCHMARK_CONFIG}")
	ENDFOREACH(BENCHMARK_CONFIG)
	
	UNSET (COPY_BENCHMARKS_TO_REF CACHE)
ENDIF (BENCHMARK_DIR_FOUND)

IF (OGS_BUILD_UTILITIES)
	ADD_SUBDIRECTORY (UTL/MSHGEOTOOLS/)
	ADD_SUBDIRECTORY (UTL/FileConverter/)
	ADD_SUBDIRECTORY (tests/GEO/)
ENDIF (OGS_BUILD_UTILITIES)
