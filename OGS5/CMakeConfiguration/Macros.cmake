# check 64 bit
if( CMAKE_SIZEOF_VOID_P EQUAL 4 )
 set( HAVE_64_BIT 0 )
 set( BITS 32 )
else( CMAKE_SIZEOF_VOID_P EQUAL 4 )
 set( HAVE_64_BIT 1 )
 set( BITS 64)
endif( CMAKE_SIZEOF_VOID_P EQUAL 4 )


# Visual Studio detection
if (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" 
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008"
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10"
  OR ${CMAKE_GENERATOR} STREQUAL "NMake Makefiles")
  
  set (VS32 TRUE)
  set (VS64 FALSE)
  message (STATUS "Generator: Visual Studio 32 Bit")

endif (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" 
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008"
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10"
  OR ${CMAKE_GENERATOR} STREQUAL "NMake Makefiles")
  
if (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005 Win64" 
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008 Win64"
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10 Win64")
  
  set (VS32 FALSE)
  set (VS64 TRUE)
  message (STATUS "Generator: Visual Studio 64 Bit")

endif (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005 Win64" 
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008 Win64"
  OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10 Win64")

# Convert environment variables
if (NOT $ENV{LIBRARIES_DIR} STREQUAL "")
	STRING(REGEX REPLACE "\\\\" "/" LIBRARIES_DIR $ENV{LIBRARIES_DIR}) 
endif (NOT $ENV{LIBRARIES_DIR} STREQUAL "")

MACRO(COPY_FILE_IF_CHANGED in_file out_file target)
    IF(${in_file} IS_NEWER_THAN ${out_file})    
        #message("Copying file: ${in_file} to: ${out_file}")
        ADD_CUSTOM_COMMAND (
                TARGET     ${target}
                POST_BUILD
                COMMAND    ${CMAKE_COMMAND}
                ARGS       -E copy ${in_file} ${out_file}
        )
        ENDIF(${in_file} IS_NEWER_THAN ${out_file})
ENDMACRO(COPY_FILE_IF_CHANGED)

MACRO(COPY_FILE_INTO_DIRECTORY_IF_CHANGED in_file out_dir target)
        GET_FILENAME_COMPONENT(file_name ${in_file} NAME)
        COPY_FILE_IF_CHANGED(${in_file} ${out_dir}/${file_name}
${target})      
ENDMACRO(COPY_FILE_INTO_DIRECTORY_IF_CHANGED)

#Copies all the files from in_file_list into the out_dir. 
# sub-trees are ignored (files are stored in same out_dir)
MACRO(COPY_FILES_INTO_DIRECTORY_IF_CHANGED in_file_list out_dir target)
    FOREACH(in_file ${in_file_list})
                COPY_FILE_INTO_DIRECTORY_IF_CHANGED(${in_file}
${out_dir} ${target})
        ENDFOREACH(in_file)     
ENDMACRO(COPY_FILES_INTO_DIRECTORY_IF_CHANGED)

MACRO(COPY_FILE_INTO_EXECUTABLE_DIRECTORY in_file target)
	IF (WIN32)
		COPY_FILE_INTO_DIRECTORY_IF_CHANGED(
			${CMAKE_CURRENT_SOURCE_DIR}/${in_file}
			${EXECUTABLE_OUTPUT_PATH}/Debug
			${target}
		)
		COPY_FILE_INTO_DIRECTORY_IF_CHANGED(
			"${CMAKE_CURRENT_SOURCE_DIR}/${in_file}"
			${EXECUTABLE_OUTPUT_PATH}/Release
			${target}
		)
	ELSE (WIN32)
		COPY_FILE_INTO_DIRECTORY_IF_CHANGED(
			${CMAKE_CURRENT_SOURCE_DIR}/${in_file}
			${EXECUTABLE_OUTPUT_PATH}
			${target}
		)
	ENDIF (WIN32)
ENDMACRO(COPY_FILE_INTO_EXECUTABLE_DIRECTORY)

# Adds a benchmark run.
# authorName Your short name
# benchmarkName Relative path in benchmarks directory
# ogsConfiguration E.g. "OGS_FEM"
# Additional arguments add output files to compare
FUNCTION (ADD_BENCHMARK authorName benchmarkName ogsConfiguration)
  
  SET (CONFIG_MATCH FALSE)
  IF (${ogsConfiguration} STREQUAL "OGS_FEM" AND OGS_FEM)
    SET (CONFIG_MATCH TRUE)
  ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM" AND OGS_FEM)
  IF (${ogsConfiguration} STREQUAL "OGS_FEM_SP" AND OGS_FEM_SP)
    SET (CONFIG_MATCH TRUE)
  ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_SP" AND OGS_FEM_SP)
  IF (${ogsConfiguration} STREQUAL "OGS_FEM_GEMS" AND OGS_FEM_GEMS)
    SET (CONFIG_MATCH TRUE)
  ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_GEMS" AND OGS_FEM_GEMS)
  IF (${ogsConfiguration} STREQUAL "OGS_FEM_BRNS" AND OGS_FEM_BRNS)
    SET (CONFIG_MATCH TRUE)
  ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_BRNS" AND OGS_FEM_BRNS)
  IF (${ogsConfiguration} STREQUAL "OGS_FEM_PQC" AND OGS_FEM_PQC)
    SET (CONFIG_MATCH TRUE)
  ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_PQC" AND OGS_FEM_PQC)
  IF (UNIX) # Only supported on Linux
    IF (${ogsConfiguration} STREQUAL "OGS_FEM_LIS" AND OGS_FEM_LIS)
      SET (CONFIG_MATCH TRUE)
    ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_LIS" AND OGS_FEM_LIS)
    IF (${ogsConfiguration} STREQUAL "OGS_FEM_MKL" AND OGS_FEM_MKL)
      SET (CONFIG_MATCH TRUE)
    ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_MKL" AND OGS_FEM_MKL)
    IF (${ogsConfiguration} STREQUAL "OGS_FEM_MPI" AND OGS_FEM_MPI)
      SET (CONFIG_MATCH TRUE)
    ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_MPI" AND OGS_FEM_MPI)
  ENDIF (UNIX)

  IF (CONFIG_MATCH)
    IF (WIN32)
      SET (ogsExe ${EXECUTABLE_OUTPUT_PATH}/Release/ogs)
    ELSE (WIN32)
      SET (ogsExe ${EXECUTABLE_OUTPUT_PATH}/ogs)
    ENDIF(WIN32)
    
    STRING (REGEX MATCH "[^/]+$" benchmarkStrippedName ${benchmarkName})
    STRING (LENGTH ${benchmarkName} benchmarkNameLength)
    STRING (LENGTH ${benchmarkStrippedName} benchmarkStrippedNameLength)
    MATH (EXPR substringLength ${benchmarkNameLength}-${benchmarkStrippedNameLength}) 
    STRING (SUBSTRING ${benchmarkName} 0 ${substringLength} benchmarkDir)
    STRING (REPLACE "/" "_" benchmarkNameUnderscore ${benchmarkName})
    STRING (REPLACE "_LONG_" "_" benchmarkNameUnderscore ${benchmarkNameUnderscore})
    STRING (REPLACE "_EXCEEDING_" "_" benchmarkNameUnderscore ${benchmarkNameUnderscore})

    # Delete output files on testing
    FOREACH (entry ${ARGN})
	SET (FILES_TO_DELETE "${FILES_TO_DELETE} ${entry}")
    ENDFOREACH (entry ${ARGN})

    # Adds a benchmark run. This calls AddTest.cmake to execute several steps.
    ADD_TEST (
      ${authorName}_BENCHMARK_${benchmarkName}
      ${CMAKE_COMMAND}
	-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
	-DEXECUTABLE_OUTPUT_PATH=${EXECUTABLE_OUTPUT_PATH}
	-DbenchmarkStrippedName=${benchmarkStrippedName}
	-DbenchmarkDir=${benchmarkDir}
	-DFILES_TO_DELETE=${FILES_TO_DELETE}
	-P ${PROJECT_SOURCE_DIR}/CMakeConfiguration/AddBenchmark.cmake
    )
    
    # compare file differences with python script only on dev2.intern.ufz.de
	IF(HOST_IS_DEV2)
    	IF (PYTHONINTERP_FOUND)
    	  FILE (REMOVE ${PROJECT_SOURCE_DIR}/../benchmarks/results/temp/temp_${benchmarkNameUnderscore}.txt)
    	  FOREACH (entry ${ARGN})
    	    FILE (APPEND ${PROJECT_SOURCE_DIR}/../benchmarks/results/temp/temp_${benchmarkNameUnderscore}.txt "${entry}\n")
    	  ENDFOREACH (entry ${ARGN})
    	  ADD_TEST (
    	    ${authorName}_FILECOMPARE_${benchmarkName}  
    	    ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/../benchmarks/results
    	    ${PYTHON_EXECUTABLE}
    	    ${PROJECT_SOURCE_DIR}/scripts/compare.py
    	    temp/temp_${benchmarkNameUnderscore}.txt
    	    ../../benchmarks_ref/
    	    ${authorName}_${benchmarkNameUnderscore}.html
    	    ../
    	  )
    	# compare files with builtin cmake command
    	ELSE (PYTHONINTERP_FOUND)
    	  FOREACH (OUTPUTFILE ${ARGN})
    	    ADD_TEST (
    	      ${authorName}_FILECOMPARE_${OUTPUTFILE}
    	      ${CMAKE_COMMAND} -E compare_files ${PROJECT_SOURCE_DIR}/../benchmarks/${OUTPUTFILE} ${PROJECT_SOURCE_DIR}/../benchmarks_ref/${OUTPUTFILE}
    	    )
    	  ENDFOREACH (OUTPUTFILE ${ARGN})
    	ENDIF (PYTHONINTERP_FOUND)
	ENDIF(HOST_IS_DEV2)
  
  # copy benchmark output files to reference directory
  IF (COPY_BENCHMARKS_TO_REF)
    FOREACH (entry ${ARGN})
      CONFIGURE_FILE( ${PROJECT_SOURCE_DIR}/../benchmarks/${entry} ${PROJECT_SOURCE_DIR}/../benchmarks_ref/${entry} COPYONLY)
    ENDFOREACH (entry ${ARGN})
  ENDIF (COPY_BENCHMARKS_TO_REF)
  
  ENDIF (CONFIG_MATCH)   
  
ENDFUNCTION (ADD_BENCHMARK authorName benchmarkName ogsConfiguration filesToCompare)

# Checks if a valid ogs configuration is given
FUNCTION(CHECK_CONFIG)
	
	SET(configs
		"${OGS_USE_QT}"
		"${OGS_FEM}"
		"${OGS_FEM_SP}"
		"${OGS_FEM_MPI}"
		"${OGS_FEM_GEMS}"
		"${OGS_FEM_BRNS}"
		"${OGS_FEM_MKL}"
		"${OGS_FEM_PQC}"
		"${OGS_FEM_LIS}"
		"${OGS_FEM_CHEMAPP}")
	
	SET(counter 0)
	
	FOREACH(config ${configs})
		IF (config)
			MATH(EXPR counter "${counter} + 1")
		ENDIF (config)
	ENDFOREACH()
	IF (counter EQUAL 0)
		MESSAGE(STATUS "No configuration specified. Assuming default configuration...")
		SET(OGS_FEM ON)
	ENDIF (counter EQUAL 0)
	
	IF (counter GREATER 1)
		MESSAGE(FATAL_ERROR "Error: More than one OGS configuration given. Please use only one of the following configurations:
			OGS_USE_QT (GUI configuration)
			OGS_FEM (Default FEM configuration)
			OGS_FEM_SP
			OGS_FEM_MPI
			OGS_FEM_GEMS
			OGS_FEM_BRNS
			OGS_FEM_MKL
			OGS_FEM_PQC
			OGS_FEM_LIS
			OGS_FEM_CHEMAPP")
	ENDIF (counter GREATER 1)
	
ENDFUNCTION()

# Sets variable PROCESSOR_COUNT
# Can be used with: cmake --build .
FUNCTION(CHECK_PROCESSOR_COUNT)
	
	IF(NOT DEFINED PROCESSOR_COUNT)

	  # Unknown:
	  SET(PROCESSOR_COUNT 0)

	  # Linux:
	  SET(cpuinfo_file "/proc/cpuinfo")
	  IF(EXISTS "${cpuinfo_file}")
	    FILE(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
	    LIST(LENGTH procs PROCESSOR_COUNT)
	  ENDIF()

	  # Mac:
	  IF(APPLE)
	    FIND_PROGRAM(cmd_sys_pro "system_profiler")
	    IF(cmd_sys_pro)
	      EXECUTE_PROCESS(COMMAND ${cmd_sys_pro} OUTPUT_VARIABLE info)
	      STRING(REGEX REPLACE "^.*Total Number Of Cores: ([0-9]+).*$" "\\1"
	        PROCESSOR_COUNT "${info}")
	    ENDIF()
	  ENDIF()

	  # Windows:
	  IF(WIN32)
	    SET(PROCESSOR_COUNT "$ENV{NUMBER_OF_PROCESSORS}")
	  ENDIF()

	ENDIF()
	
	IF(PROCESSOR_COUNT)
	  SET(CTEST_BUILD_FLAGS "-j${PROCESSOR_COUNT}")
	ENDIF()
	
	# MESSAGE(STATUS "No. of processors: ${PROCESSOR_COUNT}")
	
ENDFUNCTION()
