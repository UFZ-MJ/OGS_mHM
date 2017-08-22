# Source files

set( SOURCES
	admin.h
	msh_core.h
	msh_edge.h
	msh_elem.h
	msh_lib.h
	msh_mesh.h
	msh_node.h
	MSHEnums.h
	MeshNodesAlongPolyline.h
	MeshQualityChecker.h
	
	msh_core.cpp
	msh_edge.cpp
	msh_elem.cpp
	msh_gen.cpp
	msh_lib.cpp
	msh_mat.cpp
	msh_mesh.cpp
	msh_node.cpp
	MSHEnums.cpp
	MeshNodesAlongPolyline.cpp
	MeshQualityChecker.cpp
	#msh_pnt.cpp		# references removed elements.h nodes.h
	
	Delaunay/DT_Matrix.h
	Delaunay/tricommon.h
	Delaunay/tri_algebra.cpp
	Delaunay/tri_geometry.cpp
	Delaunay/tri_main.cpp
	Delaunay/tri_triangulation.cpp
	Delaunay/tri_utility.cpp
	
	#NW/delaunay2D.h
	#NW/delaunay3D.h
	#NW/delaunay2D.cpp
	#NW/delaunay3D.cpp
	
	#TM/dtm_calculate.h
	#TM/dtm_crowd.h
	#TM/dtm_element.h
	#TM/dtm_error.h
	#TM/dtm_figure.h
	#TM/dtm_fixed.h
	#TM/dtm_gpnode.h
	#TM/dtm_gptetra.h
	#TM/dtm_gptriangle.h
	#TM/dtm_group.h
	#TM/dtm_laplas.h
	#TM/dtm_neighbor.h
	#TM/dtm_node.h
	#TM/dtm_plane.h
	#TM/dtm_point.h
	#TM/dtm_stdio.h
	#TM/dtm_surface.h
	#TM/dtm_tetgen.h
	#TM/dtm_tetra.h
	#TM/dtm_timer.h
	#TM/dtm_triangle.h
	#TM/dtmesh.h
	#TM/dtm_calculate.cpp
	#TM/dtm_crowd.cpp
	#TM/dtm_error.cpp
	#TM/dtm_figure.cpp
	#TM/dtm_gpnode.cpp
	#TM/dtm_gptetra.cpp
	#TM/dtm_gptriangle.cpp
	#TM/dtm_laplas.cpp
	#TM/dtm_neighbor.cpp
	#TM/dtm_node.cpp
	#TM/dtm_plane.cpp
	#TM/dtm_point.cpp
	#TM/dtm_stdio.cpp
	#TM/dtm_surface.cpp
	#TM/dtm_tetgen.cpp
	#TM/dtm_tetra.cpp
	#TM/dtm_triangle.cpp
	#TM/dtmesh.cpp
	#TM/main.cpp
)

include_directories(
	.
	../Base
	../GEO
	../MathLib
	../FileIO
	../FEM
	../LIB
	../OGSProject
)

# Create the library
add_library( MSH STATIC ${SOURCES} )

target_link_libraries (
	MSH
	FileIO
	GEO
)
