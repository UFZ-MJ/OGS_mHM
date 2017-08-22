# Source files
set( SOURCES
	gs_project.cpp
)


# Header files
set( HEADERS
	gs_project.h
)

# and include additional header
include_directories(
	../Base
	../MathLib
	../GEO
	../FileIO
	../MSH
	../FEM
)


# Create the library
add_library( OGSProject
	${SOURCES}
	${HEADERS}
)

#target_link_libraries( OGSProject
#
#)
