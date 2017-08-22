@ECHO OFF
ECHO Please enter a path where to store the OGS libraries (e.g. "C:\Libs"):
SET /P LIBDIR=[Enter path]
mkdir %LIBDIR%
cd %LIBDIR%
ECHO Downloading QtPatcher ...
wget http://dl.dropbox.com/u/9346484/QtPatcher.exe

:: Qt
mkdir Qt
cd Qt
ECHO Downloading Qt ...
wget http://dl.dropbox.com/u/9346484/qt-4.6.3_VS2005.7z
7za x -y qt-4.6.3_VS2005.7z
del /S /Q qt-4.6.3_VS2005.7z
:: wget http://dl.dropbox.com/u/9346484/qt-4.6.3-min_VS2005.7z
:: 7za x -y qt-4.6.3-min_VS2005.7z
:: del /S /Q qt-4.6.3-min_VS2005.7z
cd ..

:: VTK
mkdir VTK
cd VTK
ECHO Downloading VTK ...
wget http://dl.dropbox.com/u/9346484/vtk-5.6.0_VS2005.7z
7za x -y vtk-5.6.0_VS2005.7z
del /S /Q vtk-5.6.0_VS2005.7z
cd ..

:: GeoTiff
ECHO Downloading GeoTIFF
wget http://dl.dropbox.com/u/9346484/libgeotiff_VS2005.7z
7za x -y libgeotiff_VS2005.7z
del /S /Q libgeotiff_VS2005.7z

::Tiff
ECHO Downloading TIFF
wget http://dl.dropbox.com/u/9346484/libtiff_VS2005.7z
7za x -y libtiff_VS2005.7z
del /S /Q libtiff_VS2005.7z

:: Shapelib
ECHO Downloading Shapelib
wget http://dl.dropbox.com/u/9346484/shapelib_VS2005.7z
7za x -y shapelib_VS2005.7z
del /S /Q shapelib_VS2005.7z

:: Set environment variables
wget http://dl.dropbox.com/u/9346484/SetEnv.exe
setenv -ua OGS_LIBS %LIBDIR%
setenv -ua QTDIR ~OGS_LIBS~\Qt
setenv -ua Path %"~OGS_LIBS~\Qt\bin;~OGS_LIBS~\VTK\bin"
setenv -ua QMAKESPEC ~OGS_LIBS~\Qt\mkspecs\win32-msvc2005
del /S /Q SetEnv.exe

:: Run QtPatcher
QtPatcher.exe