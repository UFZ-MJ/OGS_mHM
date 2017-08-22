cmake_minimum_required( VERSION 2.8 )

find_package( GTest REQUIRED )

set( TFILES
  a.bc
  a.gli
  a.ic
  a.mfp
  a.mmp
  a.msh
  a.msp
  a.num
  a.out
  a.pcs
  a.tim
  )

UPDATE_MODEL_FILES ( 
  ${PROJECT_BINARY_DIR}/tests/data/bmskel
  TFILES )
