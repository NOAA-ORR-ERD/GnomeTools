#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from netCDF4 import Dataset, MFDataset
from matplotlib import path
from libgoods import nctools, curv_model


class roms(curv_model.curv):

    """
    A class for loadiing variables from ROMS model output on specified domain subsets that are
    needed for GNOME
  
    """
         
    def get_dimensions(self,tvar='ocean_time',get_time=True,get_xy=True):
         
        if get_time:
            self.time_varname = tvar
            self.time = self.Dataset.variables[tvar][:]
            self.time_units = self.Dataset[tvar].units
            self.time_dimension = self.Dataset.variables[tvar].dimensions
        
        if get_xy:
            if self.GridDataset is not None:
                self.lon = self.GridDataset.variables['lon_rho'][:]
                self.lat = self.GridDataset.variables['lat_rho'][:]
            else:
                self.lon = self.Dataset.variables['lon_rho'][:]
                self.lat = self.Dataset.variables['lat_rho'][:]
        
            self.x = [0,self.lon.shape[1]]
            self.y = [0,self.lat.shape[0]]
    
    def write_nc(self,ofn='test.nc',t_index=None,is3d=False,grid_only=False):       
        '''

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
        
        nc_out.createDimension('eta_rho',self.y[1]-self.y[0])
        nc_out.createDimension('xi_rho',self.x[1]-self.x[0])
        nc_out.createDimension('eta_psi',self.y[1]-self.y[0]-1)
        nc_out.createDimension('xi_psi',self.x[1]-self.x[0]-1)
        nc_out.createDimension('eta_u',self.y[1]-self.y[0])
        nc_out.createDimension('xi_u',self.x[1]-self.x[0]-1)
        nc_out.createDimension('eta_v',self.y[1]-self.y[0]-1)
        nc_out.createDimension('xi_v',self.x[1]-self.x[0])
        
        if is3d:
            nc_out.createDimension('s_rho',nc_grid.dimensions['s_rho'].size)
        
        # List of variables to copy (they have to be in nc_in...):
        grid_vars = ['lon_rho','lat_rho','mask_rho','lon_psi','lat_psi',\
        'mask_psi','lon_u','lat_u','mask_u','lon_v','lat_v','mask_v','angle']
        
        if is3d:
            grid_vars.extend(['hc','Cs_r','s_rho','h'])

        for var in grid_vars:
            # Create variable in new file:            
            var_in  = nc_grid.variables[var]
            #print(var_in)
            try:
                dim1 = var_in.dimensions[0]     
            except IndexError:
                print('exception with variable:', var_in)
                var_out = nc_out.createVariable(var, var_in.dtype, ())
                var_out[:] = var_in[:]
            else:  
                if dim1.find('s_rho')>=0:           
                    var_out = nc_out.createVariable(var, var_in.dtype, ('s_rho'))
                    # Copy data:
                    var_out[:] = var_in[:]
                elif dim1.find('_rho')>=0:          
                    var_out = nc_out.createVariable(var, var_in.dtype, ('eta_rho','xi_rho'))
                    # Copy data:
                    var_out[:] = var_in[self.y[0]:self.y[1],self.x[0]:self.x[1]]
                elif dim1.find('_psi')>=0:          
                    var_out = nc_out.createVariable(var, var_in.dtype, ('eta_psi','xi_psi'))
                    # Copy data:
                    var_out[:] = var_in[self.y[0]:self.y[1]-1,self.x[0]:self.x[1]-1]
                elif dim1.find('_u')>=0:          
                    var_out = nc_out.createVariable(var, var_in.dtype, ('eta_u','xi_u'))
                    # Copy data:
                    var_out[:] = var_in[self.y[0]:self.y[1],self.x[0]:self.x[1]-1]
                elif dim1.find('_v')>=0:          
                    var_out = nc_out.createVariable(var, var_in.dtype, ('eta_v','xi_v'))
                    # Copy data:
                    var_out[:] = var_in[self.y[0]:self.y[1]-1,self.x[0]:self.x[1]]               
                else:
                    print('Variable dimensions could not be determined - skipping',var)
                    continue

            # Copy NetCDF attributes:
            for attr in var_in.ncattrs():
                var_out.setncattr(attr, var_in.getncattr(attr))
        
        if not grid_only:
            data_vars = ['u','v','temp','salt']
            
            for var in data_vars:
                # Create variable in new file:
                print(var)
                var_in  = nc_in.variables[var]
                coords = var_in.coordinates
                
                if is3d:
                    dims = (self.time_dimension[0],var_in.dimensions[1],var_in.dimensions[2],var_in.dimensions[3])
                    if coords.find('lon_u')>-1:  
                        var_out = nc_out.createVariable(var, var_in.dtype, dims)
                    # Copy data:
                        var_out[:] = var_in[t_index[0]:t_index[1]:t_index[2],:,self.y[0]:self.y[1],self.x[0]:self.x[1]-1]
                    elif coords.find('lon_v')>-1:          
                        var_out = nc_out.createVariable(var, var_in.dtype, dims)
                        # Copy data:
                        var_out[:] = var_in[t_index[0]:t_index[1]:t_index[2],:,self.y[0]:self.y[1]-1,self.x[0]:self.x[1]]
                    elif coords.find('lon_rho')>-1:          
                        var_out = nc_out.createVariable(var, var_in.dtype, dims)
                        # Copy data:
                        var_out[:] = var_in[t_index[0]:t_index[1]:t_index[2],:,self.y[0]:self.y[1],self.x[0]:self.x[1]]
                    else:
                        print('Variable dimensions could not be determined - skipping',var)
                        continue
                
                else:
                    dims = (self.time_dimension[0],var_in.dimensions[2],var_in.dimensions[3])
                    if coords.find('lon_u')>-1:  
                        var_out = nc_out.createVariable(var, var_in.dtype, dims)
                        # Copy data:
                        var_out[:] = var_in[t_index[0]:t_index[1]:t_index[2],-1,self.y[0]:self.y[1],self.x[0]:self.x[1]-1]
                    elif coords.find('lon_v')>-1:          
                        var_out = nc_out.createVariable(var, var_in.dtype, dims)
                        # Copy data:
                        var_out[:] = var_in[t_index[0]:t_index[1]:t_index[2],-1,self.y[0]:self.y[1]-1,self.x[0]:self.x[1]]
                    else:
                        print('Variable dimensions could not be determined - skipping',var)
                        continue
                        
                # Copy NetCDF attributes:
                for attr in var_in.ncattrs():
                    if attr == 'coordinates':
                        var_out.setncattr(attr, coords.replace('s_rho ',''))
                    else:
                        if attr != '_FillValue':
                            var_out.setncattr(attr, var_in.getncattr(attr))

        nc_out.close()
    

            
    
    def interp_and_rotate(self,u,v,is3d=False):
        '''Calculate u/v on rho points -- we lose exterior most u/v values
        Then rotate to north/east
        '''
        #self.grid['angle'] = self.grid['angle'][1:-1,1:-1]
        #self.grid['mask'] = self.grid['mask'][1:-1,1:-1]
        cosa = (np.cos(self.grid['angle'][1:-1,1:-1]) * self.grid['mask_rho'][1:-1,1:-1])
        sina = (np.sin(self.grid['angle'][1:-1,1:-1]) * self.grid['mask_rho'][1:-1,1:-1])
        if is3d:
            u_rot = np.zeros([u.shape[0],u.shape[1],v.shape[2]-1,v.shape[2]-1])
            v_rot = np.zeros_like(u_rot)
            for z in u.shape[1]:
                u_rho = (u[:,z,1:-1,:-1] +u[:,z,:-1,1:])/2. 
                v_rho = (v[:,z,:-1,1:-1] + v[:,z,1:,1:-1])/2.
                u_rot[:,z,:,:] = u_rho * cosa - v_rho * sina
                v_rot[:,z,:,:] = u_rho * sina + v_rho * cosa
        else:
            u_rho = (u[...,1:-1,:-1] + u[...,1:-1,1:])/2. 
            v_rho = (v[...,:-1,1:-1] + v[...,1:,1:-1])/2.
            u_rot = u_rho * cosa - v_rho * sina
            v_rot = u_rho * sina + v_rho * cosa
            
        return u_rot, v_rot
        
    def reduce_latlon_mesh_for_GNOME(self):
           
        if 'lon_ss' in self.data: #subset
            self.data['lon_ss'] = self.data['lon_ss'][:-1,:-1]
            self.data['lat_ss'] = self.data['lat_ss'][:-1,:-1]   
        else:
            self.data['lon'] = self.data['lon'][:-1,:-1]
            self.data['lat'] = self.data['lat'][:-1,:-1]
            
        #self.grid['mask'] = self.grid['mask_rho'][1:-1,1:-1]

    def write_nc_gnome1(self,ofn='test.nc',t_index=None,is3d=False,grid_only=False):       
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
        #nc_out.createDimension('yc',self.y[1]-self.y[0]-1)
        #nc_out.createDimension('xc',self.x[1]-self.x[0]-1)
        #nc_out.createDimension('s_rho',1)
        
        # List of variables to potentially copy
        grid_vars = ['lon','lat','mask','lonc','latc']
        roms_vars = ['lon_psi','lat_psi','mask_psi','lon_rho','lat_rho']
        for ii,var in enumerate(grid_vars):
            # Create variable in new file:r
            var_in  = nc_grid.variables[roms_vars[ii]]           
            var_out = nc_out.createVariable(var, var_in.dtype, ('y','x'))
            # Copy data
            print(var)
            print(var_in[self.y[0]:self.y[1],self.x[0]:self.x[1]].shape)
            print(var_out.shape)
            var_out[:] = var_in[self.y[0]:self.y[1],self.x[0]:self.x[1]]          
            # Copy NetCDF attributes:
            for attr in var_in.ncattrs():
                var_out.setncattr(attr, var_in.getncattr(attr))

        if not grid_only:
            
            t1 = tindex[0]; t2 = tindex[1]; ts = tindex[2]
         
            y1 = self.y[0]; y2 = self.y[1]; step = self.y[2]
            x1 = self.x[0]; x2 = self.x[1]

            u = self.Dataset.variables['u']
            u.set_auto_maskandscale(False)
            
            v = self.Dataset.variables['v']
            v.set_auto_maskandscale(False)
                                  
            if len(u.shape) == 3: #no z dimension
                u_on_upts = u[t1:t2+1:ts,y1:y2+1,x1:x2]
                v_on_vpts = v[t1:t2+1:ts,y1:y2,x1:x2+1]          
            elif is3d: 
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
            u,v = self.interp_and_rotate(u_on_upts,v_on_vpts)
            
            dims = (self.time_dimension[0],u.dimensions[1],u.dimensions[2],u.dimensions[3])
            
            var_out = nc_out.createVariable('water_u', u.dtype, dims)
            var_out[:] = u

            # Copy NetCDF attributes:
            for attr in u.ncattrs():
                if attr != '_FillValue':
                    var_out.setncattr(attr, var_in.getncattr(attr))
                              
        nc_out.close()
        
    def get_grid_info(self,yindex=None,xindex=None,is3d=False):
        
        if xindex is None and yindex is None:
            x1 = 0; x2 = self.data['lon_psi'].shape[1]; step = 1
            y1 = 0; y2 = self.data['lon_psi'].shape[0]
        else:
            y1 = yindex[0]; y2 = yindex[1]; step = yindex[2]
            x1 = xindex[0]; x2 = xindex[1]
        
        #temp fix, need to look at this (mask)
        self.grid['mask'] = self.Dataset.variables['mask_rho'][y1+1:y2:step,x1+1:x2:step] 
        self.grid['mask_rho'] = self.Dataset.variables['mask_rho'][y1:y2+1:step,x1:x2+1:step] 
        self.grid['angle'] = self.Dataset.variables['angle'][y1:y2+1:step,x1:x2+1:step] 
        
        if is3d:
            self.grid['h'] = self.Dataset.variables['h'][y1:y2+1:step,x1:x2+1:step] 
            self.grid['hc'] = self.Dataset.variables['hc'][:]
            self.grid['Cs_r'] = self.Dataset.variables['Cs_r'][:]
            self.grid['sc_r'] = self.Dataset.variables['s_rho'][:]
            
        #load lat/lon for rho, u, and v grids
        #TODO: I don't think is right indexing for u/v grids (not used yet)
        for var in ['lat_rho','lat_u','lat_v','lon_rho','lon_u','lon_v']:
        #for var in ['lat_rho','lon_rho']:
            ds_var = self.Dataset.variables[var]
            self.atts[var] = {}
            for an_att in ds_var.ncattrs():
                self.atts[var][an_att] = getattr(ds_var,an_att)
            self.data[var] = ds_var[y1:y2+1:step,x1:x2+1:step] 
            
        self.data['lonc'] = self.data['lon_rho'][1:-1,1:-1]
        self.data['latc'] = self.data['lat_rho'][1:-1,1:-1]

    
    def get_data(self,var_map,tindex=None,yindex=None,xindex=None,is3d=False,interp=True):
        
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
            x1 = 0; x2 = self.data['lon_psi'].shape[1]; step = 1
            y1 = 0; y2 = self.data['lon_psi'].shape[0]
        else:
            y1 = yindex[0]; y2 = yindex[1]; step = yindex[2]
            x1 = xindex[0]; x2 = xindex[1]
            self.data['lon_ss'] = self.data['lon_psi'][y1:y2:step,x1:x2:step]
            self.data['lat_ss'] = self.data['lat_psi'][y1:y2:step,x1:x2:step]
        
        u = self.Dataset.variables['u']
        u.set_auto_maskandscale(False)
        self.atts['u'] = {}
        for an_att in u.ncattrs():
            self.atts['u'][an_att] = getattr(u,an_att) 
        v = self.Dataset.variables['v']
        v.set_auto_maskandscale(False)
        self.atts['v'] = {}
        for an_att in v.ncattrs():
            self.atts['v'][an_att] = getattr(v,an_att) 

        if len(u.shape) == 3: #no z dimension
            u_on_upts = u[t1:t2+1:ts,y1:y2+1,x1:x2]
            v_on_vpts = v[t1:t2+1:ts,y1:y2,x1:x2+1]          
        elif is3d: 
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
               
        if interp:
            self.data['u'],self.data['v'] = self.interp_and_rotate(u_on_upts,v_on_vpts)
            # self.data['lon_rho'] = self.data['lon_rho'][1:-1,1:-1]
            # self.data['lat_rho'] = self.data['lat_rho'][1:-1,1:-1]
            # self.grid['mask_rho'] = self.grid['mask_rho'][1:-1,1:-1]
        else:
            self.data['u'] = u_on_upts
            self.data['v'] = v_on_vpts