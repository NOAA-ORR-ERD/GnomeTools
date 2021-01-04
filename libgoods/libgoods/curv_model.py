#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from netCDF4 import Dataset, MFDataset
from matplotlib import path
from libgoods import nctools, base

class curv(base.nc):
         
    def get_dimensions(self,var_map,get_time=True,get_xy=True,get_z=False):
        '''
        Get the model dimensions (time,x,y)
        Can get just time dimension, just xy, or everything
        '''
        
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
            
            self.x = [0,self.lon.shape[1]]
            self.y = [0,self.lat.shape[0]]
            self.lon = (self.lon > 180).choose(self.lon,self.lon-360)
            
        if get_z:
            self.depth = self.Dataset.variables[var_map['z']]
            self.z = [0,self.depth.shape[0]]
            
           
    def subset(self,bbox,stride=1,dl=True):
        '''
        bbox = [slat,wlon,nlat,elon]
   
        '''
        glat = self.lat
        glon = self.lon
        
        sl = bbox[0]
        nl = bbox[2]
        wl = bbox[1]
        el = bbox[3]
        
        if (abs(np.nanmax(glat)-nl) < 1e-3) and (abs(np.nanmin(glon)-wl) < 1e-3): #original values
            self.y = [0,np.size(glat,0),1]
            self.x = [0,np.size(glat,1),1]
        else: #do subset
                     
            if not dl:
                [yvec,xvec] = np.where(np.logical_and(np.logical_and(glat>=sl,glat<=nl),np.logical_and(glon>=wl,glon<=el)))
            else:
                [yvec,xvec] = np.where(np.logical_and(np.logical_and(glat>=sl,glat<=nl),np.logical_or(glon>=wl,glon<=el)))
                #self.dlx = 1 #Could do this to make the writing not need the dl flag passed in explicity
                
            if len(yvec) > 2 and len(xvec) > 2:
                y1 = min(yvec)
                y2 = max(yvec)+1
                x1 = min(xvec)
                x2 = max(xvec)+1
                self.y = [y1,y2,stride]
                self.x = [x1,x2,stride]           
            else:
                self.y = [0,np.size(glat,0),1]
                self.x = [0,np.size(glat,1),1]
               
    def write_nc_gnome1(self,var_map,ofn='test.nc',t_index=None,is3d=False,grid_only=False,zind=-1,dl=False,wind=False):       
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
        setattr(nc_out,'grid_type','curv')
        
        if not grid_only:
            if t_index is not None:
                self.time = self.time[t_index[0]:t_index[1]:t_index[2]]
            else:
                t_index = [0,len(self.time),1]               
            nc_out.createDimension(self.time_dimension[0],None)
            tvar = nc_out.createVariable(self.time_varname, 'f8', self.time_dimension)
            tvar[:] = self.time[:].data
            tvar.setncattr('units',self.time_units)
        
        nc_out.createDimension('y',self.y[1]-self.y[0])
        nc_out.createDimension('x',self.x[1]-self.x[0])

        #nc_out.createDimension('s_rho',1)
        

        # List of variables to potentially copy
        grid_vars = ['lon','lat','mask','lonc','latc']

        for var in grid_vars:
            if var in var_map:
                # Create variable in new file:
                var_in  = nc_grid.variables[var_map[var]]           
                var_out = nc_out.createVariable(var, var_in.dtype, ('y','x'))
                # Copy data
                if var in ['lon','lonc']: #Need to only do this when not crossing -180
                    lon = var_in[self.y[0]:self.y[1],self.x[0]:self.x[1]]
                    if dl is True:
                        lon = (lon < 0).choose(lon,lon+360)                  
                    var_out[:] =  lon         
                # Copy NetCDF attributes:
                else:  
                    var_out[:] = var_in[self.y[0]:self.y[1],self.x[0]:self.x[1]]          
                # Copy NetCDF attributes:
                for attr in var_in.ncattrs():
                    var_out.setncattr(attr, var_in.getncattr(attr))

        if not grid_only:
                   
            for var in self.data_vars:
                if var in var_map:
                    # Create variable in new file:
                    var_in  = nc_in.variables[var_map[var]]
                    coords = var_in.coordinates
                    dims = (self.time_dimension[0],'y','x')
                    
                    if var.find('velocity') > 0:
                        uv = var.split('_')[0]
                        if wind == True:
                            var_out = nc_out.createVariable('air_' + uv, var_in.dtype, dims)
                        else:
                            var_out = nc_out.createVariable('water_' + uv, var_in.dtype, dims)
                        
                    # Copy data:
                    if len(var_in.shape) == 4:    
                        var_out[:] = var_in[t_index[0]:t_index[1]:t_index[2],zind,self.y[0]:self.y[1],self.x[0]:self.x[1]]
                    else:
                        var_out[:] = var_in[t_index[0]:t_index[1]:t_index[2],self.y[0]:self.y[1],self.x[0]:self.x[1]]
                    # Copy NetCDF attributes:
                    for attr in var_in.ncattrs():
                        if attr == 'coordinates':
                            var_out.setncattr(attr, coords.replace('z ',''))
                        else:
                            if attr != '_FillValue':
                                var_out.setncattr(attr, var_in.getncattr(attr))
                                      
        nc_out.close()
        
        
    def write_nc(self,var_map,ofn='test.nc',t_index=None,is3d=False,grid_only=False,zind=-1,dl=True):       
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
            tvar = nc_out.createVariable(self.time_varname, 'f8', self.time_dimension)
            tvar[:] = self.time[:].data
            tvar.setncattr('units',self.time_units)
        
        nc_out.createDimension('y',self.y[1]-self.y[0])
        nc_out.createDimension('x',self.x[1]-self.x[0])

        #nc_out.createDimension('s_rho',1)
        

        # List of variables to potentially copy
        grid_vars = ['lon','lat','mask','lonc','latc']

        for var in grid_vars:
            if var in var_map:
                # Create variable in new file:
                var_in  = nc_grid.variables[var_map[var]]           
                var_out = nc_out.createVariable(var_map[var], var_in.dtype, ('y','x'))
                # Copy data
                if var in ['lon','lonc']: #Need to only do this when not crossing -180
                    lon = var_in[self.y[0]:self.y[1],self.x[0]:self.x[1]]
                    if dl is True:
                        lon = (lon < 0).choose(lon,lon+360)
                    if flipud:
                        print('flipping')
                        var_out[:] = lon.transpose()
                    else:
                        var_out[:] =  lon         
                # Copy NetCDF attributes:
                else:
                    if flipud:
                        var_out[:] = var_in[self.y[0]:self.y[1],self.x[0]:self.x[1]].transpose()
                    else:
                        var_out[:] = var_in[self.y[0]:self.y[1],self.x[0]:self.x[1]]          
                # Copy NetCDF attributes:
                for attr in var_in.ncattrs():
                    var_out.setncattr(attr, var_in.getncattr(attr))

        if not grid_only:
                   
            for var in self.data_vars:
                if var in var_map:
                    # Create variable in new file:
                    var_in  = nc_in.variables[var_map[var]]
                    coords = var_in.coordinates
                    dims = (self.time_dimension[0],'y','x')

                    var_out = nc_out.createVariable(var_map[var], var_in.dtype, dims)
                    # Copy data:
                    if len(var_in.shape) == 4:    
                        if flipud:
                            for i,tt in enumerate(range(t_index[0],t_index[1],t_index[2])):
                                var_out[i,:] = np.flipud(var_in[tt,zind,self.y[0]:self.y[1],self.x[0]:self.x[1]])
                        else:
                            var_out[:] = var_in[t_index[0]:t_index[1]:t_index[2],zind,self.y[0]:self.y[1],self.x[0]:self.x[1]]
                    else:
                        if flipud:
                            for i,tt in enumerate(range(t_index[0],t_index[1],t_index[2])):
                                var_out[i,:] = np.flipud(var_in[tt,self.y[0]:self.y[1],self.x[0]:self.x[1]])
                        else:
                            var_out[:] = var_in[t_index[0]:t_index[1]:t_index[2],self.y[0]:self.y[1],self.x[0]:self.x[1]]
                    # Copy NetCDF attributes:
                    for attr in var_in.ncattrs():
                        if attr == 'coordinates':
                            var_out.setncattr(attr, coords.replace('z ',''))
                        else:
                            if attr != '_FillValue':
                                var_out.setncattr(attr, var_in.getncattr(attr))
                    
                    # Add location attribute
                    # c = coords.split(' ')[-1]
                    # if c == var_map['lon']:
                        # var_out.setncattr("location", "node")
                        
                    # elif c == var_map['lonc']:
                        # var_out.setncattr("location", "center")


        # grid_topo = nc_out.createVariable('grid',np.dtype(int),())
        # grid_topo.cf_role = "grid_topology"
        # grid_topo.topology_dimension = 2 ;
        # grid_topo.node_coordinates = "ulon ulat"
        # grid_topo.face_coordinates = "lon lat"
        # grid_topo.node_dimensions = "y x"
        # grid_topo.face_dimensions = "y: y (padding: low) x: x (padding: low)"                            
            
        nc_out.close()
        
        


    
