# Source files

set( SOURCES
	Stiff_Bulirsch-Stoer.h
	break.h
	display.h
	eos.h
	fem_ele.h
	fem_ele_std.h
	fem_ele_vec.h
	FEMCondition.h
	FEMEnums.h
	femlib.h
	files0.h
	#intrface.h
	makros.h
	mathlib.h
	matrix_class.h
	memory.h
	par_ddc.h
	pcs_dm.h
	problem.h
	prototyp.h
	rf_bc_new.h
	rf_fct.h
	rf_fluid_momentum.h
	rf_ic_new.h
	rf_kinreact.h
	rf_mfp_new.h
	rf_mmp_new.h
	rf_msp_new.h
	rf_node.h
	rf_num_new.h
	rf_out_new.h
	rf_pcs.h
	rf_random_walk.h
	rf_react.h
	rf_st_new.h
	rf_tim_new.h
	rfmat_cp.h
	solver.h
	timer.h
	tools.h
	vtk.h
	GeoInfo.h
	ProcessInfo.h
	DistributionInfo.h
	GeoInfo.cpp
	ProcessInfo.cpp
	DistributionInfo.cpp
	Stiff_Bulirsch-Stoer.cpp
	break.cpp
	display.cpp
	eos.cpp
	fem_ele.cpp
	fem_ele_std.cpp
	fem_ele_vec.cpp
	femlib.cpp
	FEMCondition.cpp
	FEMEnums.cpp
	files0.cpp
	#intrface.cpp
	mathlib.cpp
	matrix_class.cpp
	memory.cpp
	par_ddc.cpp
	pcs_dm.cpp
	problem.cpp
	rf_bc_new.cpp
	rf_fct.cpp
	rf_fluid_momentum.cpp
	rf_ic_new.cpp
	rf_kinreact.cpp
	rf_mfp_new.cpp
	rf_mmp_new.cpp
	rf_msp_new.cpp
	rf_node.cpp
	rf_num_new.cpp
	rf_out_new.cpp
	rf_pcs.cpp
	rf_random_walk.cpp
	rf_react.cpp
	rf_st_new.cpp
	rf_tim_new.cpp
	rfmat_cp.cpp	
	timer.cpp
	tools.cpp
	vtk.cpp
)

## New sparse class. WW
IF(OGS_FEM_SP OR OGS_FEM_MPI OR OGS_FEM_GEMS OR OGS_FEM_BRNS OR OGS_FEM_PQC OR OGS_FEM_MKL OR OGS_FEM_LIS)
	SET( SOURCES ${SOURCES} equation_class.h equation_class.cpp )
ELSE(OGS_FEM_SP OR OGS_FEM_MPI OR OGS_FEM_GEMS OR OGS_FEM_BRNS OR OGS_FEM_PQC OR OGS_FEM_MKL OR OGS_FEM_LIS)
	SET( SOURCES ${SOURCES} solver.h solver.cpp matrix.h matrix.cpp)
ENDIF(OGS_FEM_SP OR OGS_FEM_MPI OR OGS_FEM_GEMS OR OGS_FEM_BRNS OR OGS_FEM_PQC OR OGS_FEM_MKL OR OGS_FEM_LIS)


IF(OGS_FEM_GEMS)
	SET( SOURCES ${SOURCES} rf_REACT_GEM.h rf_REACT_GEM.cpp ) 
ENDIF( OGS_FEM_GEMS )

IF(OGS_FEM_BRNS)
	SET( SOURCES ${SOURCES} rf_REACT_BRNS.h rf_REACT_BRNS.cpp ) 
ENDIF(OGS_FEM_BRNS)

include_directories(
	.
	../Base
	../MathLib
	../GEO
	../FileIO
	../MSH
	../GEM
	../OGSProject
)

IF (LIBPQC)
  INCLUDE_DIRECTORIES(
    ../PQC
    )
  ADD_DEFINITIONS ( -DLIBPHREEQC )
  IF (LIBPQC_DEBUG)
    ADD_DEFINITIONS ( -DMDL_DEBUG )
  ENDIF (LIBPQC_DEBUG)
ENDIF (LIBPQC)

# Create the library
add_library( FEM STATIC ${SOURCES} )

target_link_libraries( 
	FEM
	Base
	FileIO
 	GEO
 	MSH
	OGSProject
)

IF (LIBPQC)
   target_link_libraries( FEM
	GEO
	MSH
	PQC
)
ENDIF (LIBPQC)

IF (OGS_FEM_CHEMAPP)
	INCLUDE_DIRECTORIES ( ../EQL )
	TARGET_LINK_LIBRARIES ( FEM EQL )
ENDIF (OGS_FEM_CHEMAPP)
