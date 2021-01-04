from __future__ import print_function
import numpy as np
from netCDF4 import Dataset, MFDataset
from libgoods import base

class rect(base.nc):
    """
    A class for dealing with regular grid model output and converting to GNOME format
    By passing in a var_map dict the 
    variable names can be customized for different models or datasets
            
    """
            
    def get_dimensions(self,var_map,get_time=True,get_xy=True,get_z=False):
        
        if get_time:
            self.time_varname = var_map['time']
            self.time = self.Dataset.variables[var_map['time']][:]
            self.time_units = self.Dataset[var_map['time']].units
            self.time_dimension = self.Dataset.variables[var_map['time']].dimensions
        
        if get_xy:
            if self.GridDataset is not None:
                self.lon = self.GridDataset.variables[var_map['lon']][:]
                self.lat = self.GridDataset.variables[var_map['lat']][:]
            else:
                self.lon = self.Dataset.variables[var_map['lon']][:]
                self.lat = self.Dataset.variables[var_map['lat']][:]
            self.x = [0,len(self.lon),1]
            self.y = [0,len(self.lat),1]

            self.lon = (self.lon > 180).choose(self.lon,self.lon-360)
            
        if get_z:
            self.depth = self.Dataset.variables[var_map['z']]
            self.z = [0,self.depth.shape[0]]
            
        
    def subset(self,bbox,stride=1,dl=0):
        '''
        bbox = [slat,wlon,nlat,elon]

        '''       
        subset_lat = np.nonzero(np.logical_and(self.lat>=bbox[0],self.lat<=bbox[2]))[0]
            
        if dl==0:
            subset_lon = np.nonzero(np.logical_and(self.lon>=bbox[1],self.lon<=bbox[3]))[0]
        else:
            subset_lon = np.nonzero(np.logical_or(self.lon>=bbox[1],self.lon<=bbox[3]))[0]
            
        
        if stride >= len(subset_lat):
            stride = 1
            
        self.y = [subset_lat[0],subset_lat[-1]+1,stride]
        self.x = [subset_lon[0],subset_lon[-1]+1,stride]
        
    
    def write_nc(self,var_map,ofn='test.nc',t_index=None,d_index=0,grid_only=False,dl=0):       
        '''
        Write a curvilinear grid with u/v on center points (may or may not have grid
        cell lat/lons)
        '''       
        nc_in = self.Dataset
        nc_in.set_auto_maskandscale(False)
        
        if self.GridDataset is not None:
            nc_grid = self.GridDataset
            nc_grid.set_auto_maskandscale(False)
        else:
            nc_grid = nc_in
            
        nc_out = Dataset(ofn,'w',format='NETCDF3_CLASSIC')
                
        if not grid_only:
            if t_index is not None:
                self.time = self.time[t_index[0]:t_index[1]:t_index[2]]
            else:
                t_index = [0,len(self.time),1]               
            nc_out.createDimension(self.time_dimension[0],None)
            tvar = nc_out.createVariable(self.time_varname, 'f4', self.time_dimension)
            tvar[:] = self.time[:]
            tvar.setncattr('units',self.time_units)
        
        nc_out.createDimension('lat',self.y[1]-self.y[0])
        nc_out.createDimension('lon',self.x[1]-self.x[0])
        
        
        var_in  = nc_grid.variables[var_map['lon']]           
        var_out = nc_out.createVariable('lon', var_in.dtype,('lon'))
        if dl == 0:            
            var_out[:] = var_in[self.x[0]:self.x[1]:self.x[2]]      
        else:
            lon = var_in[self.x[0]:self.x[1]:self.x[2]]
            var_out[:] = (lon < 0 ).choose(lon,lon+360)
        # Copy NetCDF attributes:
        for attr in var_in.ncattrs():
            var_out.setncattr(attr, var_in.getncattr(attr))
                
        var_in  = nc_grid.variables[var_map['lat']]           
        var_out = nc_out.createVariable('lat', var_in.dtype,('lat'))
        var_out[:] = var_in[self.y[0]:self.y[1]:self.y[2]]          
        # Copy NetCDF attributes:
        for attr in var_in.ncattrs():
            var_out.setncattr(attr, var_in.getncattr(attr))
            
        if 'mask' in var_map:
            var_in  = nc_in.variables[var_map['mask']]
            coords = var_in.coordinates
            dims = (var_in.dimensions[-2],var_in.dimensions[-1])
            var_out = nc_out.createVariable('mask', var_in.dtype, dims)
            var_out[:] = var_in[self.y[0]:self.y[1]:self.y[2],self.x[0]:self.x[1]:self.x[2]]
            # Copy NetCDF attributes:
            for attr in var_in.ncattrs():
                var_out.setncattr(attr, var_in.getncattr(attr))
                
        if not grid_only:
            
            #determine which data is desired from var_map
            data_vars = []
            for var in self.data_vars:
                if var in var_map:        
                    data_vars.append(var_map[var])
            
            for var in data_vars:
                # Create variable in new file:
                var_in  = nc_in.variables[var]
                dims = (self.time_dimension[0],var_in.dimensions[-2],var_in.dimensions[-1])
                var_out = nc_out.createVariable(var, var_in.dtype, dims)
                # Copy data:
                if len(var_in.shape)==4:
                    var_out[:] = var_in[t_index[0]:t_index[1]:t_index[2],d_index,self.y[0]:self.y[1],self.x[0]:self.x[1]]
                elif len(var_in.shape)==3:
                    var_out[:] = var_in[t_index[0]:t_index[1]:t_index[2],self.y[0]:self.y[1],self.x[0]:self.x[1]]

                # Copy NetCDF attributes:
                if type(self.Dataset) is MFDataset:
                    var_in = Dataset(self.FileName[0]).variables[var]
                for attr in var_in.ncattrs():
                    if attr != '_FillValue':
                        var_out.setncattr(attr, var_in.getncattr(attr))
                                 

        nc_out.close()
        
