# Source files
set( HEADERS
	binarySearch.h
	Configure.h.in
	DateTools.h
	FileFinder.h
	FileTools.h
	quicksort.h
	printList.h
	uniqueListInsert.h
	StringTools.h
	swap.h
	wait.h	
)

set( SOURCES
	DateTools.cpp
	StringTools.cpp
	binarySearch.cpp
)

# Create the library
add_library( Base STATIC ${HEADERS} ${SOURCES} )

set_target_properties(Base PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	../GEO
	../MathLib
	.
)

#target_link_libraries (
#        Base
#        GEO
#)


IF(OGS_BUILD_CORRECT_SCIENTIFIC)
	# Create utility executable
	add_executable( correctScientificNotation CorrectScientificNotationMain.cpp )
ENDIF(OGS_BUILD_CORRECT_SCIENTIFIC)

