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
#import matplotlib
#matplotlib.style.use('ggplot')


if (os.name == 'posix'):
    raster_lib_path = os.environ['HOME'] + '/Python_Scripts/rasterpy'
elif (os.name == 'nt'):
    ogs_lib_path = 'C:/Users/group/Source/lib/rasterpy'
sys.path.append(raster_lib_path)

from rasterpy import *

prodir = '/home/miao/Documents/Naegelstedt/Simulation/SC1_GW_ex_old'

outfile = prodir + '/waterbalance_homo_monthly_hist2.png'
#model = pd.DataFrame({'Model':wattabs})
config = {'start': '1975',
            'end': '2005'}
rng = pd.date_range(config['start'],config['end'],freq='M')
#model.index = pd.DatetimeIndex(pd.Series(rng))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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

  
for i in range(612-3-12*30,612-3):
    mHM = Raster()
    raster = mHM.rasterin(prodir+'/mHMhcbf{0}.asc'.format(i)) #read discharge
    invol, outvol = mHM.waterbalance(raster,method='mean')
    inwaterlist.append(invol+2.5)
    outwaterlist.append(outvol)

inwater = np.array(inwaterlist)
outwater = np.array(outwaterlist)   
data = pd.DataFrame({'recharge':inwater,'baseflow':outwater})
data.index = pd.DatetimeIndex(pd.Series(rng))
data = data.resample("m").mean()

#fig=plt.figure(figsize=(10,7.5))
#ax = fig.add_subplot(111)  
## Ensure that the axis ticks only show up on the bottom and left of the plot.  
## Ticks on the right and top of the plot are generally unnecessary chartjunk.  
#ax.get_xaxis().tick_bottom()
#ax.get_yaxis().tick_left()

# Make sure your axis ticks are large enough to be easily read.  
# You don't want your viewers squinting to read your plot. 
#ax.tick_params(labelsize=16)
ax =data.plot.hist(alpha=0.6, stacked=False, fontsize=20, figsize=(10,8))
# plot reference lines
#ax.bar(data.index,data['recharge'],  label='groundwater recharge')
#ax.bar(data.index,data['discharge'],  label='river baseflow')
##ax.plot(new.index, new['Model'], c=tableau20[4],ls='--',lw=3, label ='simulated')
#ax.text('Mean value of baseflow_mHM: {0}')
ax.set_xlabel(r"Value [mm]",fontsize=28)

ax.set_ylabel("Frequency",fontsize=28)
ax.grid(ls='--')
legend = ax.legend(loc='upper right', shadow=False, fontsize = 28)
#
#
plt.savefig(outfile,dpi=300)