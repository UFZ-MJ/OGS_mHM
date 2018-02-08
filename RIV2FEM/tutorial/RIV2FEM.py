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

timeindex = (576, 600)

for ind in range(timeindex[0],timeindex[1]+1):
    mHM = Raster()
    #gwrech = mHM.rasterin('/home/miao/Documents/Naegelstedt/Raster/recharge-monthly.asc') #read recharge
    gwdis = mHM.rasterin('/home/miao/Naegelstedt/Raster/timeseries/original/mHMbfrc{0}'.format(ind)) #read discharge
    gwdis = mHM.extract(gwdis,threshold = 50)
    
    #agwdis = mHM.clipdata(gwdis, domain=[0,31])
    # remove areas close to boundaries
    gwdis = mHM.removedata(gwdis,rowdomain=[7,8],coldomain=[12,19])
    gwdis = mHM.removedata(gwdis,rowdomain=[20,21],coldomain=[18,23])
    #gwdis = mHM.removedata(gwdis,rowdomain=[20,22],coldomain=[30,32])
    gwdis = mHM.removedata(gwdis,rowdomain=[28,30],coldomain=[35,35])
    
    gwdis = mHM.scaleraster(gwdis,-1)
    
    #gwdis = mHM.removedata(gwdis,rowdomain=[56,61],coldomain=[55,58])
    
    #newgwdis = mHM.modifyraster(gwdis, 0.2) #modify discharge by multiplying a factor
    #newgwdis = mHM.logfit(gwdis, domain = [5,1210], valuerange = [5, 135])
    
    #mergedras = mHM.mergeraster(discharge=gwdis,recharge=0)
    #maskedarr = mHM.maskraster(gwdis)
    
    #mHM.plot(maskedarr)
    
    mHM.writeraster(gwdis,'/home/miao/Naegelstedt/Raster/timeseries/cliped/mHMbfrc{0}'.format(ind))
