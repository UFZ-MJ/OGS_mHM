# Source files

set( SOURCES
	array.h
	databr.h
	datach.h
	gdatastream.h
	gstring.h
	io_arrays.h
	jama_cholesky.h
	jama_lu.h
	m_const.h
	m_param.h
	ms_multi.h
	node.h
	num_methods.h
	s_fgl.h
	tnt.h
	tnt_array1d.h
	tnt_array2d.h
	tnt_i_refvec.h
	v_user.h
	verror.h
)


# Create the library
add_library( GEM STATIC ${SOURCES} )

SET_TARGET_PROPERTIES(GEM PROPERTIES LINKER_LANGUAGE CXX)

# Link Qt library
#target_link_libraries( GEM
#	${QT_LIBRARIES}
#)
