# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 20:37:59 2017

@author: miao
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

if (os.name == 'posix'):
    raster_lib_path = os.environ['HOME'] + '/Python_Scripts/rasterpy'
elif (os.name == 'nt'):
    ogs_lib_path = 'C:/Users/group/Source/lib/rasterpy'
sys.path.append(raster_lib_path)
#matplotlib.style.use('ggplot')

from rasterpy import *

prodir = '/home/miao/Documents/Naegelstedt/Simulation/SC1_GW_ex_old'
outfile = prodir + '/waterbalance_timeseries_1.png'

#model = pd.DataFrame({'Model':wattabs})
config = {'start': '1975',
            'end': '2005'}
rng = pd.date_range(config['start'],config['end'],freq='M')
#model.index = pd.DatetimeIndex(pd.Series(rng))

inwaterlist = []
outwaterlist = []
#
#for i in range(600-12*30,600):
#    mHM = Raster()
#    gwdis = mHM.rasterin('/home/miao/Documents/Naegelstedt/Raster/River_Baseflow/timeseries/original/mHMbfrc{0}'.format(i)) #read discharge
#    gwdis = mHM.clipdata(gwdis,domain=[30,10000])
#    gwdis = mHM.calrivermean(gwdis)
#    # remove areas close to boundaries
#    #gwdis = mHM.removedata(gwdis,rowdomain=[7,8],coldomain=[12,19])
#    gwdis = mHM.removedata(gwdis,rowdomain=[19,21],coldomain=[16,20])
#    #gwdis = mHM.removedata(gwdis,rowdomain=[20,22],coldomain=[30,32])
#    #gwdis = mHM.removedata(gwdis,rowdomain=[28,30],coldomain=[35,35])
#    gwdis = mHM.scaleraster(gwdis,-.04)
#    	 
#    #print gwdis
#    gwrec = mHM.rasterin('/home/miao/Documents/Naegelstedt/Raster/Recharge/transience/temp/mHM_recahrge_{0}.ASC'.format(i))
#    gwrec = mHM.scaleraster(gwrec,.001)
#    	
#    mergedras = mHM.mergeraster(gwdis,gwrec)
#    invol, outvol = mHM.waterbalance(mergedras)
#    inwaterlist.append(invol)
#    outwaterlist.append(outvol)
    	#mHM.writeraster(mergedras,'/home/miao/Naegelstedt/Raster/MergedBaseflowandGWRecharge/original/mHMbfgr{0}'.format(i))

#for i in range(600-12*30,340):
#    mHM = Raster()
#    raster = mHM.rasterin('/home/miao/Documents/Naegelstedt/Raster/MergedBaseflowandGWRecharge/original/mHMbfgr{0}'.format(i)) #read discharge
#    invol, outvol = mHM.waterbalance(raster)
#    inwaterlist.append(invol)
#    outwaterlist.append(outvol)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#tableau20 = colorscheme("tableau20")
    
for i in range(612-3-12*30,612-3):
    mHM = Raster()
    raster = mHM.rasterin(prodir+'/mHMhcbf{0}.asc'.format(i)) #read discharge
    invol, outvol = mHM.waterbalance(raster,method='mean')
    inwaterlist.append(invol+2.5)
    outwaterlist.append(outvol)


inwater = np.array(inwaterlist)
outwater = np.array(outwaterlist)  
to3out = outwater[:-2]
last3out = outwater[-2:]
outwater = np.concatenate([last3out,to3out])
outwater = (outwater-outwater.mean())*1.2 + outwater.mean()
 
data = pd.DataFrame({'recharge':inwater,'discharge':outwater})
data.index = pd.DatetimeIndex(pd.Series(rng))
#data = data.resample("a", how='mean')

fig=plt.figure(figsize=(15,5))  
ax = fig.add_subplot(111)  
# Ensure that the axis ticks only show up on the bottom and left of the plot.  
# Ticks on the right and top of the plot are generally unnecessary chartjunk.  
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Make sure your axis ticks are large enough to be easily read.  
# You don't want your viewers squinting to read your plot. 
ax.tick_params(labelsize=16)

# plot reference lines
ax.set_ylabel(r"Value [mm]",fontsize=20)
#ax.set_xlabel(r"Monitoring well index",fontsize=24)
ax.plot(data.index,data['discharge'], ls='-',lw=2, label='baseflow')

ax.plot(data.index,data['recharge'], ls='--',lw=2, label='recharge')
#ax.plot(new.index, new['Model'], c=tableau20[4],ls='--',lw=3, label ='simulated')
ax.set(ylim=(-1,40))  
ax.set(xlim=('1975','2005'))  
ax.grid(ls='--')
legend = ax.legend(loc='upper right', shadow=False, fontsize = 22)

plt.savefig(outfile,dpi=300)


