# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 20:37:59 2017

@author: miao
"""
import os
import sys
import numpy as np

if (os.name == 'posix'):
    raster_lib_path = os.environ['HOME'] + '/Python_Scripts/rasterpy'
elif (os.name == 'nt'):
    ogs_lib_path = 'C:/Users/group/Source/lib/rasterpy'
sys.path.append(raster_lib_path)

from rasterpy import *

#mHM = Netcdf(ncfile='/home/miao/Naegelstedt/Raster/netcdf/baseflow_monthly_mm_mon-1.nc',var='routed_baseflow')
mHM = Netcdf(ncfile='/home/miao/Documents/Naegelstedt/Raster/Precipitation/pre.nc',var='routed_baseflow')

mHM.rasterin('/home/miao/Documents/Naegelstedt/Raster/mHMbfrc') # a necessary step to define self.nrows, self.ncols
for index in range(576,576+37):
    raster = mHM.nc2ascii(timeindex=index)
    #mHM.writeraster(raster,'/home/miao/Naegelstedt/Raster/River_Baseflow/timeseries/original/mHMbfrc{0}'.format(index))
    print raster
#mHM.writeraster(mergedras,'/home/miao/Documents/Naegelstedt/Raster/mHMbfrc2')
