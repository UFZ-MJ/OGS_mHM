## OGS-5 ##
see https://svn.ufz.de/ogs for build instructions

## For those who are interested in ogs-pardiso on linux, make sure you have not only checked out the sources/ directory but also the Libs/ directory and then compile with the option

cd sources
mkdir Build
cd Build
cmake .. -DOGS_FEM_MKL=ON

NOTE: Since the executable, ogs, is created by linking the libraries dynamically, users need to add the Libs/MKL/[Bits]/ directory to their shell environment. Otherwise, the executable complains missing libraries.

## For OGS-GEM user, use the following cmake option,
cd sources
mkdir Build
cd Build
cmake .. -DOGS_FEM_GEMS=ON 

## For using the GUI use the following cmake option ##

cmake .. -DOGS_USE_QT=ON

