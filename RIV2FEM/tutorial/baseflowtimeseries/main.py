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

mHM = Netcdf(ncfile='/home/miao/Naegelstedt/Raster/netcdf/baseflow_monthly_mm_mon-1.nc',var='routed_baseflow')

mHM.rasterin('/home/miao/Naegelstedt/Raster/mHMbfrc') #a necessary step to define self.nrows, self.ncols

bfseries = mHM.timeseries(50,40)
#mHM.writeraster(raster,'/home/miao/Naegelstedt/Raster/timeseries/mHMbfrc{0}'.format(index))
print bfseries

with open('/home/miao/Naegelstedt/Raster/row50_col40_bfts','w') as outfile:
    outfile.write('\n'.join(str(data) for data in bfseries))