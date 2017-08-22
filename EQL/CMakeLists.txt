# Source files
SET ( SOURCES 
	cacint.h
	cacint.cpp
	eqlink.h
	eqlink.cpp
)

# Create the library
ADD_LIBRARY ( EQL STATIC ${SOURCES} )

INCLUDE_DIRECTORIES (
	../Base
	../Mathlib
	../GEO
	../MSH
	../FEM
)
