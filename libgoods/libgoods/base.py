#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from netCDF4 import Dataset, MFDataset, num2date
from matplotlib import path
from libgoods import nctools

class nc():
    
    def __init__(self,FileName=None,GridFileName=None):
        
        if FileName is not None:
            self.FileName = FileName
            if isinstance(FileName,list):
                self.Dataset = MFDataset(FileName)
            else:
                self.Dataset = Dataset(FileName)
        else:
            self.Dataset = None
            
        if GridFileName is not None:
            self.GridFileName = GridFileName
            self.GridDataset = Dataset(GridFileName)
        else:
            self.GridDataset = None
        
        self.data_vars = ['u_velocity',
                          'v_velocity',
                          'ice_u',
                          'ice_v',
                          'ice_thickness',
                          'ice_fraction']
                          
        self.dlx = 0 #default does not cross dateline
        
    def update(self,FileName):
        '''
        Change nc Dataset to point to a new nc file or url without reinitializing everything (retain grid info)
        '''
        if isinstance(FileName,list):
            self.Dataset = MFDataset(FileName)
        else:
            self.Dataset = Dataset(FileName)
    
    def when(self):

        print('Start date: ', num2date(self.time[0], self.time_units))
        try:
            print('End date: ', num2date(self.time[-1], self.time_units))
        except IndexError:
            printnum2date(self.time[0], self.time_units)
     