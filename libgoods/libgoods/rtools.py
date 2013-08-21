#!/usr/bin/env python
import numpy as np
from netCDF4 import Dataset, num2date, date2num, date2index
from libgoods import nctools

class rgrid:
    """
    A class for dealing with unstructured grid model output and converting to GNOME format
    Although I use variable names consistent with FVCOM, by passing in a var_map dict the 
    variable names can be customized for SELFE or ADCIRC
    
    Right now the attribute specifying whether the elements are orderd clockwise or counter
    clockwise needs to be manually added before writing to GNOME format (GNOME requres this, 
    but its not often specified in the model output)
        
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
        
    def get_data(self,var_map,tindex=None,yindex=None,xindex=None,is3d=0):
        
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

        if len(xindex) == 1:
            if len(u.shape) == 4:
                if is3d == 1:
                    self.data['u'] = u[t1:t2:ts,:,y1:y2:step,x1:x2:step]
                    self.data['v'] = v[t1:t2:ts,:,y1:y2:step,x1:x2:step]
                else:    
                    self.data['u'] = u[t1:t2:ts,0,y1:y2:step,x1:x2:step]
                    self.data['v'] = v[t1:t2:ts,0,y1:y2:step,x1:x2:step]
            elif len(u.shape) == 3:
                self.data['u'] = u[t1:t2:ts,y1:y2:step,x1:x2:step]
                self.data['v'] = v[t1:t2:ts,y1:y2:step,x1:x2:step]
            else: #u/v variables neither 3d or 4d -- add some error catching here
                pass   
        else: #this subset "wraps around" in longitude
            x3 = xindex[1][0]; x4 = xindex[1][1]
            ax = len(u.shape)-1
            if len(u.shape) == 4:
                #this does not include option to extract 3d (surface layer only)
                #if we do get a 3d region that crosses prime meridian -- add option to extract 3d here
                self.data['u'] = np.concatenate((u[t1:t2:ts,0,y1:y2:step,x1:x2:step],u[t1:t2:ts,0,y1:y2:step,x3:x4:step]),axis=ax)
                self.data['v'] = np.concatenate((v[t1:t2:ts,0,y1:y2:step,x1:x2:step],v[t1:t2:ts,0,y1:y2:step,x3:x4:step]),axis=ax)
            elif len(u.shape) == 3:
                self.data['u'] = np.concatenate((u[t1:t2:ts,y1:y2:step,x1:x2:step],u[t1:t2:ts,y1:y2:step,x3:x4:step]),axis=ax)
                self.data['v'] = np.concatenate((v[t1:t2:ts,y1:y2:step,x1:x2:step],v[t1:t2:ts,y1:y2:step,x3:x4:step]),axis=ax)
            else: #u/v variables neither 3d or 4d -- add some error catching here
                pass 
            
        #check for nans
        try:
            ufill = self.atts['u']['_FillValue']
            vfill = self.atts['v']['_FillValue']
        except KeyError:
            try:
                ufill = self.atts['u']['missing_value']
                vfill = self.atts['v']['missing_value']
            except KeyError:
                ufill = 999.
                vfill = 999.
        self.u = (self.u == np.NaN).choose(self.u,ufill)
        self.v = (self.v == np.NaN).choose(self.v,vfill)
            
        if has_attr(self,date_line):
            self.data['lon_ss'] = (self.data['lon_ss'] < 0).choose(self.data['lon_ss'],self.data['lon_ss']+360)
            
        # Check that longitude is increasing
        if self.data['lon_ss'][0] > self.data['lon_ss'][-1]:
            sort_id = self.data['lon_ss'].argsort()
            self.data['lon_ss'] = self.data['lon_ss'][sort_id]
            if self.data['u'].ndim == 3:
                self.data['u'] = self.data['u'][:,:,sort_id]
                self.data['v'] = self.data['v'][:,:,sort_id]
            elif u.ndim == 2:
                self.data['u'] = self.data['u'][:,sort_id]
                self.data['v'] = self.data['v'][:,sort_id]

        # Check that latitude vector is from S->N so that
        if self.data['lat_ss'][0] > self.data['lat_ss'][-1]:
            sort_id = self.data['lat_ss'].argsort()
            self.data['lat_ss'] = self.data['lat_ss'][sort_id]
            if u.ndim == 3:
                self.data['u'] = self.data['u'][:,sort_id,:]
                self.data['v'] = self.data['v'][:,sort_id,:]
            elif u.ndim == 2:
                self.data['u'] = self.data['u'][sort_id,:]
                self.data['v'] = self.data['v'][sort_id,:]

        
    