# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 16:27:17 2017

@author: miao
"""

#import ufz
from .raster import *

class Netcdf(Raster):
    '''
    process with netcdf data.
    '''
    def __init__(self, ncfile, var):
        Raster.__init__(self)
        self.ncfile = ncfile
        self.var = var
        
    def nc2ascii(self,timeindex):
        
        Q = ufz.readnc(self.ncfile,var=self.var)[timeindex]

        return Q.data

    def timeseries(self,rown,coln):
        '''
        time series of a special point

        :param rown: the row index
        :param coln: the colomn index
        :return: a list
        '''
        Q = ufz.readnc(self.ncfile,var=self.var)[:,rown,coln]
        return Q
