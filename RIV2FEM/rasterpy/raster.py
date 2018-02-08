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
                if abs(float(line[col_i]) - float(self.NODATA_value))<1e-2:
                    self.mask[line_i - 6][col_i] = 1
                #if float(line[col_i]) > threshold:
                rasterarr[line_i - 6][col_i] = float(line[col_i])
        #        rasterarr[line_i - 6][col_i] = line[col_i]
        f_id.close()
        rasterarr = np.array(rasterarr)
        return rasterarr
        
    def maskraster(self, inraster):
        maskedarr = np.ma.masked_values(inraster, self.NODATA_value)
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
        #meanvalue = np.sum(inraster[inraster>0]) / rvnodenum

        for val in np.nditer(inraster, op_flags=['readwrite']):
            if val > 0:
                val[...] = meanvalue

        return inraster
    
    def calrecmean(self, inraster):
        '''
        Usage
        --------
        Calculate mean recharge value over the whole catchment).

        Input
        --------
        domain(list)  - define the domain of input value.
        valuerange(list)  -define the domain of output value.
        '''
        rvnodenum = np.count_nonzero(inraster>0)

        #meanvalue = np.amax(inraster) / rvnodenum
        meanvalue = np.sum(inraster[inraster>0]) / rvnodenum

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
                val[...]= self.NODATA_value 
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
    
    def waterbalance(self, inraster, method = 'sum'):
        effeinraster = self.maskraster(inraster)#extract effective raster
        cellnumber = effeinraster.count() #size of masked array
        if method == 'sum':
            inwatervol = np.sum(effeinraster[effeinraster>0])*500*500
            outwatervol = -np.sum(effeinraster[effeinraster<0])*500*500
            print ("Total amount of groundwater recharge: {0}m3\nTotal amount of groundwater discharge: {1}m3".format(float(inwatervol),float(outwatervol)))
            with open('waterbalancelog','a') as outfile:
                outfile.write ("Total amount of groundwater recharge: {0}m3\nTotal amount of groundwater discharge: {1}m3".format(float(inwatervol),float(outwatervol)))
        #print (targetarr)
        elif method == 'mean':
            inwatervol = np.sum(effeinraster[effeinraster>0])*1000/cellnumber
            outwatervol = -np.sum(effeinraster[effeinraster<0])*1000/cellnumber
            print ("Mean groundwater recharge: {0}mm/month\nMean groundwater discharge: {1}mm/month".format(float(inwatervol),float(outwatervol)))
            with open('waterbalancelog','a') as outfile:
                outfile.write("Mean groundwater recharge: {0}mm/month\nMean groundwater discharge: {1}mm/month".format(float(inwatervol),float(outwatervol)))
        return inwatervol, outwatervol
    
    def mask_river (self, inraster):
        '''
        Input
        ------
        inraster  -- the input raster file (not the output flx.asc file)
        '''
        #effeinraster = inraster[abs(inraster-self.NODATA_value)>1e-4] #extract effective raster
        #cellnumber = effeinraster.size
        
        valid_ras = self.maskraster(inraster)
        #valid_ras_siz = valid_ras.count()
                
        baseflow_ras = np.ma.masked_greater(valid_ras, 0) # remove all cells that are with positive values
        baseflow_mask = baseflow_ras.mask
        #baseflow_ras = np.ma.masked_values(baseflow_ras, self.NODATA_value)
        
        #ogs_flx_raster = self.rasterin([:-3]+'flx.asc')
        #ogs_baseflow = np.ma.array(ogs_flx_raster, mask = baseflow_mask)
        
        return baseflow_mask
    
    def baseflow_ogs(self, inraster, riv_mask, method = 'mean'):
        valid_ras = self.maskraster(inraster)
        valid_ras_siz = valid_ras.count()
        
        baseflow_ras = np.ma.array(valid_ras, mask=riv_mask)
        if method == 'mean':
            norm_bf = np.sum(baseflow_ras)*1000/valid_ras_siz
            print ("Mean normalized baseflow: {0}mm".format(norm_bf))
        return norm_bf    
                             
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
        
    def extract(self, accubf, oribf, thres):
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
        mask2[bf2>thres]=0 #only river is not masked
        river = np.ma.array(bf2,mask = mask2)
        
        deaccriv = river.copy()
        it = np.nditer(target,flags=['multi_index'])
        for val in it:
            if val > thres:
                tararr = target[it.multi_index[0]-1:it.multi_index[0]+2,it.multi_index[1]-1:it.multi_index[1]+2]# extract the nearest blocks
                val1 = np.sum(tararr)
                deaccriv[it.multi_index] = val1
                target[it.multi_index[0]-1:it.multi_index[0]+2,it.multi_index[1]-1:it.multi_index[1]+2] = np.ma.masked
        # original river baseflow
        oririv = np.ma.array(oribf,mask=mask2)
        deaccriv = deaccriv + oririv
        
        deaccriv = np.ma.filled(deaccriv,0)       
        #return np.ma.array(outraster,mask=rimask)
        return deaccriv
    
    def river_nodes (self, geofile, rivname):
        '''
        Usage
        ------
        Calculate number of nodes within a geometry in .gli file.
        
        '''
        with open (geofile,'r') as infile:
            content = infile.readlines()
            riv_len = 0
            
            for (ind,line) in enumerate(content):
                words = line.split()
                if words[0].startswith(rivname):
                    for line in content[ind+2:]:
                        words = line.split()
                        
                        if words[0].startswith('#POLYLINE'):
                            break
                        elif words[0].startswith('#STOP'):
                            break
                        riv_len += 1
        return riv_len
        