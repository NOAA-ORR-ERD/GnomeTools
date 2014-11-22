#!/usr/bin/env python
import numpy as np
from netCDF4 import Dataset

def latlon_psi_2_rho(self,latr,lonr):
        lonp = (lonr[0:-1,0:-1]+lonr[1:,1:])*0.5
        latp = (latr[0:-1,0:-1]+latr[1:,1:])*0.5
        return lonp, latp
        
class cgrid:
    """
    A class for dealing with curvilinear grid model output and converting to GNOME format
    Requires passing in a var_map dict so the 
    variable names can be customized for different models or datasets
            
    """
    
    def __init__(self,FileName=None):
            
        if FileName is not None:
            self.FileName = FileName
            self.Dataset = Dataset(FileName)
            self.data = dict()
            self.atts = dict()
            
    def get_dimensions(self,var_map):
        
        self.time = self.Dataset.variables[var_map['time']]
   
        self.atts['time'] = self.time.__dict__ 
        self.data['time'] = self.time[:]
        
        lat = self.Dataset.variables[var_map['latitude']]
        self.atts['lat'] = lat.__dict__
        self.data['lat'] = lat[:]
        
        lon = self.Dataset.variables[var_map['longitude']]
        self.atts['lon'] = lon.__dict__
        lon = lon[:]
        self.data['lon'] = (lon > 180).choose(lon,lon-360)
        

    
    def get_arakawa_grid_info(self,yindex=None,xindex=None):
        
        if xindex is None and yindex is None:
            x1 = 0; x2 = len(self.data['lon'])
            y1 = 0; y2 = len(self.data['lat'])
        else:
            y1 = yindex[0]; y2 = yindex[1]
            x1 = xindex[0][0]; x2 = xindex[0][1]
            
        self.grid['angle'] = self.Dataset.variables['angle'][y1+1:y2,x1+1:x2] 
        self.grid['mask'] = self.Dataset.variables['mask_rho'][y1+1:y2,x1+1:x2]
        self.grid['cosa'] = np.cos(self.grid['angle']) * self.grid['mask']
        self.grid['sina'] = np.sin(self.grid['angle']) * self.grid['mask']
        self.grid['hc'] = self.Dataset.variables['hc'][:]
        self.grid['depth'] = self.Dataset.variables['h'][y1+1:y2,x1+1:x2]
        self.grid['Cs_r'] = self.Dataset.variables['Cs_r'][:]
        self.grid['sc_r'] = self.Dataset.variables['s_rho'][:]
    
    def get_data_arakawa(self,var_map,tindex=None,yindex=None,xindex=None,is3d=0):
        
        '''
        In this case, lon/lat on psi (P) grid, u on u-grid, v on v-grid
        
        '''
        if tindex is None:
            self.data['time_ss'] = self.data['time']
            t1 = 0; t2 = len(self.data['time']); ts = 1
        else:
            t1 = tindex[0]; t2 = tindex[1]; ts = tindex[2]
            self.data['time_ss'] = self.data['time'][t1:t2:ts]
            
        if xindex is None and yindex is None:
            x1 = 0; x2 = len(self.data['lon']); step = 1
            y1 = 0; y2 = len(self.data['lat'])
        else:
            y1 = yindex[0]; y2 = yindex[1]; step = yindex[2]
            x1 = xindex[0][0]; x2 = xindex[0][1]
        
        u = self.Dataset.variables[var_map['u_velocity']]
        self.atts['u'] = u.__dict__
        v = self.Dataset.variables[var_map['v_velocity']]
        self.atts['v'] = v.__dict__

        if is3d == 1: 
            u_on_upts = u[t1:t2+1:ts,:,y1:y2+1,x1:x2]
            v_on_vpts = v[t1:t2+1:ts,:,y1:y2,x1:x2+1]
        else:
            u_on_upts = u[t1:t2+1:ts,-1,y1:y2+1,x1:x2]
            v_on_vpts = v[t1:t2+1:ts,-1,y1:y2,x1:x2+1]

        #replace nans or fill values with 0 for interpolating to rho grid
        u_on_upts = (np.isnan(u_on_upts)).choose(u_on_upts,0)
        v_on_vpts = (np.isnan(v_on_vpts)).choose(v_on_vpts,0)
        u_on_upts = (u_on_upts > 1e10).choose(u_on_upts,0)
        v_on_vpts = (v_on_vpts > 1e10).choose(v_on_vpts,0)
        
        #calculate u/v on rho points -- we lose exterior most u/v values
        u = (u_on_upts[:,:,1:-1,:-1] + u_on_upts[:,:,1:-1,1:])/2. 
        v = (v_on_vpts[:,:,:-1,1:-1] + v_on_vpts[:,:,1:,1:-1])/2.
        
        #rotate vectors to north/east
        u_rot = np.zeros(np.shape(u)); v_rot = np.zeros(np.shape(v));
        for ii in range(u.shape[1]):
            u_rot[:,ii,:,:] = u[:,ii,:,:] * self.grid['cosa'] - v[:,ii,:,:] *self.grid['sina']
            v_rot[:,ii,:,:] = u[:,ii,:,:] * self.grid['sina'] + v[:,ii,:,:] * self.grid['cosa']