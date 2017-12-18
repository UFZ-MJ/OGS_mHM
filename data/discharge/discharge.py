# -*- coding: utf-8 -*-
"""
Created on Thu Aug 03 22:59:15 2017

@author: miao
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 11:00:43 2016

@author: miao
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
#from sklearn.linear_model import LinearRegression as LinReg

#plt.style.use('ggplot')

font = {'family': 'sans-serif',
        'weight': 'normal',
        'size': 14}

matplotlib.rc('font', **font)

path = '/home/miao/Documents/Naegelstedt/Raster/Discharge/daily_discharge.out'
outfile = '/home/miao/Documents/Naegelstedt/Raster/Discharge/discharge1.png'
#lottype = 'L'  # Scatter or Line

# ==============================================================================
#  read the csv file as dataframe
# ==============================================================================
df = pd.read_table(path,delim_whitespace=True)
#df = pd.to_datetime(df)

df_date = df['Day'].apply(str)+'-' + df['Mon'].apply(str)+'-'+df['Year'].apply(str)
df_date = pd.to_datetime(df_date,format="%d-%m-%Y")
df.index = pd.DatetimeIndex(df_date)
df=df.drop(['No', 'Day', 'Mon','Year'], axis=1)

df_month = df.resample('M').mean()
df_month1 = df_month['1975':]
Rcor = df_month1.corr()

obsmean = df_month1['Qobs_0009304070'].mean()
#obsmedi = df_month1['Qobs_0009304070'].median()
modmean = df_month1['Qsim_0009304070'].mean()
#modmedi = df_month1['Qsim_0009304070'].median()

#Nash–Sutcliffe model efficiency coefficient
df_month1['NS1'] = (df_month1['Qsim_0009304070']-df_month1['Qobs_0009304070'])**2
df_month1['NS2'] = (df_month1['Qobs_0009304070']-obsmean)**2

NF = 1 - df_month1['NS1'].sum()/df_month1['NS2'].sum()



plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig=plt.figure(figsize=(20,4.5))  
ax = fig.add_subplot(111)  

# Ensure that the axis ticks only show up on the bottom and left of the plot.  
# Ticks on the right and top of the plot are generally unnecessary chartjunk.  
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Make sure your axis ticks are large enough to be easily read.  
# You don't want your viewers squinting to read your plot. 
ax.tick_params(labelsize=16)

# plot reference lines
ax.plot(df_month1.index, df_month1['Qsim_0009304070'], ls='-',lw=3, label='simulated')
ax.plot(df_month1.index, df_month1['Qobs_0009304070'], ls='--',lw=3, label ='observed')

legend = ax.legend(loc='upper right', shadow=True, fontsize = 20)
ax.set_ylabel(r"Discharge Q [m^{3} s^{-1}]",fontsize=20)

ax.set(xlim=('1975','2005'))
ax.set(ylim=(0,22))    
ax.text('1976',17, r'NSE=%.2f'%NF + '\n'+ r'R_{cor}=%.2f'%Rcor['Qsim_0009304070'][0], fontsize=22)

plt.grid(ls='--')
plt.savefig(outfile,dpi=300)
plt.show()