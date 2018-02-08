# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 17:34:15 2017

@author: miao
"""
import numpy as np
from scipy.optimize import curve_fit
#from pylab import *

class Raster(object):
    def __init__(self):
        self.header = ['' for ii in range(6)] #define the outfile header
    
    def rasterin(self, rasterfile):        
        '''
        Usage
        --------
        read raster file into a 2D numpy array.
        
        Output
        --------
        rasterarr -- a masked numpy array.
        '''
        f_id = open(rasterfile)
        #f_lines = f_id.readlines()
        for line_i in range(6):
            line = f_id.next().strip()
            self.header[line_i] = line + '\n'
            line = line.split()
            if   line[0] == 'ncols':
                self.ncols = int(line[1])         
            elif line[0] == 'nrows':
                self.nrows = int(line[1])
            elif line[0] == 'xllcorner':
                self.xllcorner = float(line[1])
            elif line[0] == 'yllcorner':
                self.yllcorner = float(line[1])
            elif line[0] == 'NODATA_value':
                self.NODATA_value = int(line[1])
        rasterarr = np.zeros((self.nrows, self.ncols))
        self.mask = np.zeros((self.nrows, self.ncols))
        
        for line_i in range(6, 6 + self.nrows):
            line = f_id.next().strip()
            line = line.split()
            for col_i in range(0, self.ncols):
                if abs(float(line[col_i]) - float(self.NODATA_value))<1e-4:
                    self.mask[line_i - 6][col_i] = 1
                #if float(line[col_i]) > threshold:
                rasterarr[line_i - 6][col_i] = float(line[col_i])
        #        rasterarr[line_i - 6][col_i] = line[col_i]
        f_id.close()
        rasterarr = np.array(rasterarr)
        return rasterarr
        
    def maskraster(self, inraster):
        maskedarr = np.ma.array(inraster, mask = self.mask)
        return maskedarr
        
    def scaleraster(self, inraster, factor):
        '''
        Usage
        --------
        Scale raster file by multiplying a factor.
        
        Input
        --------
        inraster - a np array of raster values.
        
        Output
        --------
        outsterarr -- a masked numpy array.
        '''
        outraster = inraster
        for val in np.nditer(outraster, op_flags=['readwrite']):
            if abs(val-self.NODATA_value) > 1e-4:
                val[...] = val * factor
        return outraster
        
    def calrivermean(self, inraster):
        '''
        Usage
        --------
        Calculate mean baseflow value on rivers (not on whole catchment).

        Input
        --------
        domain(list)  - define the domain of input value.
        valuerange(list)  -define the domain of output value.
        '''
        rvnodenum = np.count_nonzero(inraster>0)

        meanvalue = np.amax(inraster) / rvnodenum
        
        for val in np.nditer(inraster, op_flags=['readwrite']):
            if val > 0:
                val[...] = meanvalue

        return inraster

    def logfit(self, inraster, domain = [], valuerange= []):
        '''
        Usage
        --------
        Adjust extremely large values to the overall reasonable range using logarithmic fit.
        
        Input
        --------
        domain(list)  - define the domain of input value.
        valuerange(list)  -define the domain of output value.
        '''
        def func(x,a,b):
            return a+b*np.log(x)
            
        popt, pcov = curve_fit(func, domain, valuerange)
        for val in np.nditer(inraster, op_flags=['readwrite']):
            if val > domain[0] and val < domain[1]:
                val[...] = func(val, *popt)
        return inraster

    def mergeraster(self, raster1, raster2):
        '''
        Usage
        --------
        Merge two rasters.
        
        Output
        --------
        rasterarr --Merged masked array.
        '''
        outraster = raster1 + raster2
        for val in np.nditer(outraster,op_flags=['readwrite']):
            if abs(val-2*self.NODATA_value) < 1e-4:
                val[...]= self.NODATA_value # change NODATA_value to 0
        return outraster
        
    def removedata(self, inraster, rowdomain=[5,14], coldomain=[14,26]):
        '''
        Usage
        --------
        Set 0 at points in defined domain.
        
        Input
        --------
        rowdomain(list)   - the domain of rows
        coldomain(list)   - the domain of columns
        
        Output
        --------
        rasterarr --Merged masked array.
        '''
        it = np.nditer(inraster, flags=['multi_index'], op_flags=['readwrite']) # see: https://docs.scipy.org/doc/numpy/reference/arrays.nditer.html
        for val in it:
            if it.multi_index[0]>=rowdomain[0] and it.multi_index[0]<=rowdomain[1]\
            and it.multi_index[1]>=coldomain[0] and it.multi_index[1]<=coldomain[1]:
                val[...] = 0
        return inraster
        
    def clipdata(self, inraster, domain=[]):
        '''
        Usage
        --------
        Clip data which are beyond threshold.
        
        Input
        --------
        domain(list)  -- domain of the raster value to be cliped.
        
        '''
        outraster = np.copy(inraster)
        for val in np.nditer(outraster,op_flags=['readwrite']):
            if (val<=domain[0] or val>=domain[1]) and abs(val-self.NODATA_value)>1e-2:
                val[...]= 0 # change NODATA_value to 0
        
        return outraster
    
    def waterbalance(self, inraster):
        effeinraster = inraster[abs(inraster-self.NODATA_value)>1e-4]
        inwatervol = np.sum(effeinraster[effeinraster>0])
        outwatervol = -np.sum(effeinraster[effeinraster<0])
        print ("Total amount of groundwater recharge: {0}m\nTotal amount of groundwater discharge: {1}m".format(float(inwatervol),float(outwatervol)))
        #print (targetarr)
        with open('waterbalancelog','a') as outfile:
		    outfile.write("Total amount of groundwater recharge: {0}\nTotal amount of groundwater discharge: {1}\n".format(float(inwatervol),float(outwatervol)))
        return inwatervol, outwatervol

    def writeraster(self, inraster, outfile):
        
        outarr = ''.join(x for x in self.header)
        for rowi in xrange(inraster.shape[0]):
            row = inraster[rowi,]
            outrow = '  '.join([str(col) for col in row])
            outarr += (outrow + '\n' )
        with open(outfile,'w') as out:
            out.write(outarr)
        print (' Raster file successfully written as: ' + outfile)
    
    def plot(self, inraster):

        imshow(inraster)
        colorbar()
        show()
        
    def extract(self, accubf, thres):
        '''
        Usage
        --------
        Extract value by subtracting the adjacent value. This is the counter way of accumulate. The
        inraster should be the original mHM baseflow raster.
        target - nodes around rivers (exclude river itsself)
        
        Input
        --------
        domain(list)  -- domain of the raster value to be cliped.
        '''
#        outraster = inraster
        '''
#        def find_nearest(array,value):
#            
#            idx = (np.abs(array-value)).argmin() # find the nearest value
#              
#            return array.flat[idx]

        '''            
#        it = np.nditer(inraster, flags=['multi_index'])
#        #rimask = np.zeros((self.nrows, self.ncols))
#        for val in it:
#            if abs(val-self.NODATA_value) > 1e-4 and abs(val) > threshold: #check if this point is a valid river
#                adjarr = inraster[it.multi_index[0]-1:it.multi_index[0]+2,it.multi_index[1]-1:it.multi_index[1]+2]# extract the nearest blocks
#                for candi in np.nditer(adjarr,op_flags=['readwrite']):
#                    if candi > threshold/5:
#                        candi[...] = 0
#                outraster[it.multi_index] = np.sum(adjarr)
#            elif abs(val-self.NODATA_value) > 1e-4:
#                outraster[it.multi_index] = 0
#        processedras = np.copy(clipedraster)
#        it = np.nditer(clipedraster, flags=['multi_index'])
#        for val in it:
#            if val > 0: #check if this point is a valid river
#                adjarr = clipedraster[it.multi_index[0]-1:it.multi_index[0]+2,it.multi_index[1]-1:it.multi_index[1]+2]# extract the nearest blocks
#                effecadjarr = adjarr[adjarr<val]
#                effecadjarr = effecadjarr[effecadjarr>self.NODATA_value]# focus on nodes whose values are smaller than val
#                sortedeffcadjarr = np.sort(effecadjarr)
#                if len(sortedeffcadjarr)==0 or abs(sortedeffcadjarr[-1]-0)<1e-3:
#                    processedras[it.multi_index] = 0
##                elif abs(sortedeffcadjarr[1]-0)<1e-3: # only one adjacent node has valid value
##                    processedras[it.multi_index] = val-sortedeffcadjarr[0]
#                elif np.sum(sortedeffcadjarr[-3:])<val: # adjacent value too small
#                    processedras[it.multi_index] = val-np.sum(sortedeffcadjarr[-3:])
#                elif np.sum(sortedeffcadjarr[-2:])<val: # adjacent value too small
#                    processedras[it.multi_index] = val-np.sum(sortedeffcadjarr[-2:])
##                elif sortedeffcadjarr[-1]<val*0.9:
##                    processedras[it.multi_index] = 0 # prevent abnormal values
#                else:
#                    processedras[it.multi_index] = val-sortedeffcadjarr[-1]

        bf1 = np.copy(accubf)                
        mask1 = np.zeros((self.nrows, self.ncols))+1 #mask all the nodes
        it = np.nditer(bf1, flags=['multi_index'])
        for val in it:
            if val > thres:
                #tararr = bf[it.multi_index[0]-1:it.multi_index[0]+2,it.multi_index[1]-1:it.multi_index[1]+2]# extract the nearest blocks
                mask1[it.multi_index[0]-1:it.multi_index[0]+2,it.multi_index[1]-1:it.multi_index[1]+2] = 0
        target= np.ma.array(bf1, mask = mask1)
        target= np.ma.masked_where(abs(target-self.NODATA_value)<1e-3,target)
        target= np.ma.masked_where(target>thres,target)
        
        # set up a masked array of river
        bf2 = np.copy(accubf)
        mask2 = np.zeros((self.nrows, self.ncols))+1
        mask2[bf2>thres]=0
        river = np.ma.array(bf2,mask = mask2)
        
        deaccriv = river.copy()
        it = np.nditer(target,flags=['multi_index'])
        for val in it:
            if val > thres:
                tararr = target[it.multi_index[0]-1:it.multi_index[0]+2,it.multi_index[1]-1:it.multi_index[1]+2]# extract the nearest blocks
                val1 = np.sum(tararr)
                deaccriv[it.multi_index] = val1
                target[it.multi_index[0]-1:it.multi_index[0]+2,it.multi_index[1]-1:it.multi_index[1]+2] = np.ma.masked
        deaccriv = np.ma.filled(deaccriv,0)       
        #return np.ma.array(outraster,mask=rimask)
        return deaccriv
