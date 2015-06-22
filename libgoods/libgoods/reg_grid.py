#!/usr/bin/env python
import numpy as np
from netCDF4 import Dataset

class rgrid:
    """
    A class for dealing with regular grid model output and converting to GNOME format
    By passing in a var_map dict the 
    variable names can be customized for different models or datasets
            
    """
    
    def __init__(self,FileName=None):
            
        if FileName is not None:
            self.FileName = FileName
            self.Dataset = Dataset(FileName)
            self.data = dict()
            self.atts = dict()
            
    def update(self,FileName):
        #point to a new nc file or url without reinitializing everything
        self.Dataset = Dataset(FileName)
            
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
        
    def subset(self,bbox,stride=1,dl=0,lat='lat',lon='lon'):
        '''
        bbox = [slat,wlon,nlat,elon]
        can pass in lat/lon names to specify which grid the subset is done on (for c-grids)
        '''
        
        subset_lat = np.nonzero(np.logical_and(self.data['lat']>=bbox[0],self.data['lat']<=bbox[2]))[0]
        subset_lon = np.nonzero(np.logical_and(self.data['lon']>=bbox[1],self.data['lon']<=bbox[3]))[0]
        
        if stride >= len(subset_lat):
            stride = 1
        self.y = [subset_lat[0],subset_lat[-1]+1,stride]
        self.x = [subset_lon[0],subset_lon[-1]+1,stride]
        
        
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
            x1 = xindex[0]; x2 = xindex[1] #not dealing with dateline case
            self.data['lon_ss'] = self.data['lon'][x1:x2:step]
            self.data['lat_ss'] = self.data['lat'][y1:y2:step]
        
        u = self.Dataset.variables[var_map['u_velocity']]
        self.atts['u'] = u.__dict__
        v = self.Dataset.variables[var_map['v_velocity']]
        self.atts['v'] = v.__dict__


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
        self.data['u'] = (self.data['u'] == np.NaN).choose(self.data['u'],ufill)
        self.data['v'] = (self.data['v'] == np.NaN).choose(self.data['v'],vfill)
            
#        # Check that longitude is increasing
#        if self.data['lon_ss'][0] > self.data['lon_ss'][-1]:
#            sort_id = self.data['lon_ss'].argsort()
#            self.data['lon_ss'] = self.data['lon_ss'][sort_id]
#            if self.data['u'].ndim == 3:
#                self.data['u'] = self.data['u'][:,:,sort_id]
#                self.data['v'] = self.data['v'][:,:,sort_id]
#            elif u.ndim == 2:
#                self.data['u'] = self.data['u'][:,sort_id]
#                self.data['v'] = self.data['v'][:,sort_id]
#
#        # Check that latitude vector is from S->N so that
#        if self.data['lat_ss'][0] > self.data['lat_ss'][-1]:
#            sort_id = self.data['lat_ss'].argsort()
#            self.data['lat_ss'] = self.data['lat_ss'][sort_id]
#            if u.ndim == 3:
#                self.data['u'] = self.data['u'][:,sort_id,:]
#                self.data['v'] = self.data['v'][:,sort_id,:]
#            elif u.ndim == 2:
#                self.data['u'] = self.data['u'][sort_id,:]
#                self.data['v'] = self.data['v'][sort_id,:]

        
 
    def write_nc(self,ofn,is3d=False):
        """
      
        Write GNOME compatible netCDF file (netCDF3) from regular grid data
        
        
    
        """
        
        nc = Dataset(ofn,'w',format='NETCDF3_CLASSIC')
        
        # Global Attributes
        setattr(nc,'grid_type','regular')
    
        # determine if its a subset in time
        t_key = 'time'
        try:
            if self.data['u'].shape[0] == len(self.data['time_ss']):
                t_key = 'time_ss'
        except KeyError: #TODO -- if it has key time_ss but it doesn't match that case is not caught
            if self.data['u'].shape[0] != len(self.data['time']):
                raise Exception('Dimensions of u/v do not match time variable')
                
        lon_key = 'lon'; lat_key = 'lat'
        # determine if its a subset of the grid
        try:
            lon_ss_len = len(self.data['lon_ss'])
            if self.data['u'].shape[-1] == lon_ss_len:
                lon_key = 'lon_ss'; lat_key = 'lat_ss'
        except KeyError:
            lon_len = len(self.data['lon'])
            if self.data['u'].shape[-2:] != lon_len:
                print 'Dimensions dont match'
                raise Exception('Dimensions of u/v do not match grid variables')
                
        # add Dimensions
        nc.createDimension('lon',len(self.data[lon_key]))
        nc.createDimension('lat',len(self.data[lat_key]))
        nc.createDimension('time',None)
    
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
    
        # create variables
        nc_time = nc.createVariable('time','f8',('time',))
        nc_lon = nc.createVariable('lon','f4',('lon',))
        nc_lat = nc.createVariable('lat','f4',('lat',))
        if self.atts.has_key('wind'):
            nc_u = nc.createVariable('air_u','f4',('time','lat','lon'), \
                fill_value=ufill)
            nc_v = nc.createVariable('air_v','f4',('time','lat','lon'), \
                fill_value=vfill)
        else:
            nc_u = nc.createVariable('water_u','f4',('time','lat','lon'), \
                fill_value=ufill)
            nc_v = nc.createVariable('water_v','f4',('time','lat','lon'), \
                fill_value=vfill)
        
        # add data
        nc_lon[:] = self.data[lon_key]
        nc_lat[:] = self.data[lat_key]
        if len(self.data['u'].shape) == 2:
            nc_time[0] = self.data[t_key]
            nc_u[0,:] = self.data['u']
            nc_v[0,:] = self.data['v']
        else:
            nc_time[:] = self.data[t_key]
            nc_u[:] = self.data['u']
            nc_v[:] = self.data['v']
    
        # add variable attributes from 'atts' (nested dict object)
        for an_att in self.atts['time'].iteritems():
            setattr(nc_time,an_att[0],an_att[1])
            
        for an_att in self.atts['u'].iteritems():
            if an_att[0] != '_FillValue':
                setattr(nc_u,an_att[0],an_att[1])
    
        for an_att in self.atts['v'].iteritems():
            if an_att[0] != '_FillValue':
                setattr(nc_v,an_att[0],an_att[1])
        
                
        nc.close()   