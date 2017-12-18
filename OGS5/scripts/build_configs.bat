:: Setup Visual Studio environment
call "%VS80COMNTOOLS%\..\..\VC\bin\vcvars32.bat"
call "%VS90COMNTOOLS%\..\..\VC\bin\vcvars32.bat"

:: Goto sources directory
cd ..

:: Cleanup
rd /S /Q Release build_fem build_gems build_pqc build_brns

:: Executables will copied to Release directory
mkdir Release

:: Build FEM
mkdir build_fem
cd build_fem
cmake -DOGS_FEM=ON -DOGS_DONT_USE_QT=ON ..
devenv OGS-FEM-5.sln /Build Release
cd bin\Release
copy /Y ogs.exe ..\..\..\Release\
cd ..\..\..\

:: Build FEM_GEMS
mkdir build_gems
cd build_gems
cmake -DOGS_FEM_GEMS=ON -DOGS_DONT_USE_QT=ON ..
devenv OGS-FEM-5-GEMS.sln /Build Release
cd bin\Release
ren ogs.exe ogs_gems.exe
copy /Y ogs_gems.exe ..\..\..\Release\
cd ..\..\..\

:: Build FEM_PQC
mkdir build_pqc
cd build_pqc
cmake -DOGS_FEM_PQC=ON -DOGS_DONT_USE_QT=ON ..
devenv OGS-FEM-5-PQC.sln /Build Release
cd bin\Release
ren ogs.exe ogs_pqc.exe
copy /Y ogs_pqc.exe ..\..\..\Release\
cd ..\..\..\

:: Build FEM_BRNS
mkdir build_brns
cd build_brns
cmake -DOGS_FEM_BRNS=ON -DOGS_DONT_USE_QT=ON ..
devenv OGS-FEM-5-BRNS.sln /Build Release
cd bin\Release
ren ogs.exe ogs_brns.exe
copy /Y ogs_brns.exe ..\..\..\Release\
cd ..\..\..\
