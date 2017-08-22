#!/usr/bin/env bash

# Return code. Will be set to 1 (error) when a build failed
returncode=0

# Goto sources directory
cd ..

# Executables will be copied to Release/
mkdir -p Release

# Iterate over configurations
for config in "" "SP" "MPI" "GEMS" "PQC" "BRNS" "MKL" "LIS"
do
	if [ "$config" = "" ] ; then
		config_cmake="OGS_FEM"
		exe_name="ogs"
		build_dir="build_fem"
	else
		config_low=$( echo "$config"|tr -s '[:upper:]' '[:lower:]' )
		config_cmake="OGS_FEM_$config"
		exe_name="ogs_$config_low"
		build_dir="build_$config_low"
	fi

	# Cleanup
	rm -rf $build_dir

	# Create build directory
	mkdir -p $build_dir && cd $build_dir

	# Run CMake
	../scripts/cmake.ogs.sh -D$config_cmake=ON -DOGS_DONT_USE_QT=ON -DMPI_INCLUDE_PATH=/opt/openmpi-1.4.1/include ..

	# Build
	make -j

	# Remember exit code
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi

	# Copy executable
	cp bin/ogs ../Release/$exe_name

	cd ..

done

echo "exit code is ${returncode}"
exit ${returncode}