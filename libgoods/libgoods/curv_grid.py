#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from netCDF4 import Dataset, MFDataset
from matplotlib import path
from libgoods import nctools

class cgrid():
    
    def __init__(self,FileName=None):
        
        if FileName is not None:
            self.FileName = FileName
            if isinstance(FileName,list):
                self.Dataset = MFDataset(FileName)
            else:
                self.Dataset = Dataset(FileName)
        self.data = dict()
        self.atts = dict()
        self.grid = dict()
        self.dlx = 0 #default does not cross dateline
            
    def update(self,FileName):
        #point to a new nc file or url without reinitializing everything
            if isinstance(FileName,list):
                self.Dataset = MFDataset(FileName)
            else:
                self.Dataset = Dataset(FileName)
        
    def get_dimensions(self,var_map,get_time=True,get_xy=True,get_z=False):
        '''
        Get the model dimensions (time,x,y)
        Can get just time dimension, just xy, or everything
        '''
        if get_time:
            self.time = self.Dataset.variables[var_map['time']]  
            self.atts['time'] = {}
            for an_att in self.time.ncattrs():
                self.atts['time'][an_att] = getattr(self.time,an_att) 
            self.data['time'] = self.time[:]
    
        if get_xy:
            lon = self.Dataset.variables[var_map['longitude']]
            self.atts['lon'] = {}
            for an_att in lon.ncattrs():
                self.atts['lon'][an_att] = getattr(lon,an_att)
            self.data['lon'] = lon[:]
            
            lat = self.Dataset.variables[var_map['latitude']]
            self.atts['lat'] = {}
            for an_att in lat.ncattrs():
                self.atts['lat'][an_att] = getattr(lat,an_att)
            self.data['lat'] = lat[:]      
            
        if get_z:
            z = self.Dataset.variables[var_map['z']]
            self.atts['z'] = {}
            for an_att in z.ncattrs():
                self.atts['z'][an_att] = getattr(z,an_att)
            self.data['z'] = z[:]    
            
    
    def subset_pt_in_poly(self,bbox,stride=1,lat='lat',lon='lon'):
        '''
        bbox = [slat,wlon,nlat,elon]
        can pass in lat/lon names to specify which grid the subset is done on (for c-grids)
        '''
        glat = self.data[lat]
        glon = self.data[lon]
        
        bbox=np.array(bbox)
        mypath=np.array([bbox[[1,3,3,1]],bbox[[0,0,2,2]]]).T
        p = path.Path(mypath)
        points = np.vstack((glon.flatten(),glat.flatten())).T   
        n,m = np.shape(glon)
        inside = p.contains_points(points).reshape((n,m))
        ii,jj = np.meshgrid(range(m),range(n))
        self.x = [min(ii[inside]),max(ii[inside])+1,stride]
        self.y = [min(jj[inside]),max(jj[inside])+1,stride]
        
    def subset(self,bbox,stride=1,dl=0,lat='lat',lon='lon'):
        '''
        bbox = [slat,wlon,nlat,elon]
        can pass in lat/lon names to specify which grid the subset is done on (for c-grids)
        '''
        glat = self.data[lat]
        glon = self.data[lon]
        
        sl = bbox[0]
        nl = bbox[2]
        wl = bbox[1]
        el = bbox[3]
        
        if (abs(np.nanmax(glat)-nl) < 1e-3) and (abs(np.nanmin(glon)-wl) < 1e-3): #original values
            self.y = [0,np.size(glat,0),1]
            self.x = [0,np.size(glat,1),1]
        else: #do subset
            
            if dl == 0:
                [yvec,xvec] = np.where(np.logical_and(np.logical_and(glat>=sl,glat<=nl),np.logical_and(glon>=wl,glon<=el)))
            else:
                [yvec,xvec] = np.where(np.logical_and(np.logical_and(glat>=sl,glat<=nl),np.logical_or(glon>=wl,glon<=el)))
                self.dlx = 1
                
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
               
    def get_grid_info(self,var_map,grid_vars=['mask','depth'],yindex=None,xindex=None):
    
        if xindex is None and yindex is None:
            x1 = 0; x2 = self.data['lon'].shape[1]; step = 1
            y1 = 0; y2 = self.data['lon'].shape[0]
        else:
            y1 = yindex[0]; y2 = yindex[1]; step = yindex[2]
            x1 = xindex[0]; x2 = xindex[1]
        
        for var in grid_vars: 
        #these are 2D spatial vars (need to add vertical info) for 3d
            try:
                self.grid[var] = self.Dataset.variables[var_map[var]][y1:y2:step,x1:x2:step]
            except KeyError:
                print('KeyError', var)
                pass

    
    def get_data(self,var_map,tindex=None,yindex=None,xindex=None,zindex=0,is3d=False,extra_2dvars=[]):
    
        ''' 
        var_map is a dict mapping model variable names to common names
        tindex can be used to subset in time --> tindex = [start,stop,step]
        xindex and yindex are for subsetting grid 
        zindex is for z layer (surface at zindex = 0 or -1)
        
        '''
        if tindex is None:
            self.data['time_ss'] = self.data['time']
            t1 = 0; t2 = len(self.data['time']); ts = 1
        else:
            t1 = tindex[0]; t2 = tindex[1]; ts = tindex[2]
            self.data['time_ss'] = self.data['time'][t1:t2:ts]
            
        if xindex is None and yindex is None:
            x1 = 0; x2 = self.data['lon'].shape[1]; step = 1
            y1 = 0; y2 = self.data['lon'].shape[0]
        else:
            y1 = yindex[0]; y2 = yindex[1]; step = yindex[2]
            x1 = xindex[0]; x2 = xindex[1]
            self.data['lon_ss'] = self.data['lon'][y1:y2:step,x1:x2:step]
            self.data['lat_ss'] = self.data['lat'][y1:y2:step,x1:x2:step]
        

        u = self.Dataset.variables[var_map['u_velocity']]
        u.set_auto_maskandscale(False)
        self.atts['u'] = {}
        for an_att in u.ncattrs():
            self.atts['u'][an_att] = getattr(u,an_att) 

        v = self.Dataset.variables[var_map['v_velocity']]
        v.set_auto_maskandscale(False)
        self.atts['v'] = {}
        for an_att in v.ncattrs():
            self.atts['v'][an_att] = getattr(v,an_att) 
        
        if is3d: 
            zindex=list(range(v.shape[1]))
            
        try:
            self.data['u'] = u[t1:t2:ts,zindex,y1:y2:step,x1:x2:step]
            self.data['v'] = v[t1:t2:ts,zindex,y1:y2:step,x1:x2:step]
        except ValueError: #data does not have vertical dimension
            self.data['u'] = u[t1:t2:ts,y1:y2:step,x1:x2:step]
            self.data['v'] = v[t1:t2:ts,y1:y2:step,x1:x2:step]
        
        if 'lonc' in var_map:
            lonc = self.Dataset.variables[var_map['lonc']]
            latc = self.Dataset.variables[var_map['latc']]
            self.atts['lonc'] = {}
            for an_att in lonc.ncattrs():
                self.atts['lonc'][an_att] = getattr(lonc,an_att)
            self.data['lonc'] = lonc[y1:y2:step,x1:x2:step]
            self.atts['latc'] = {}
            for an_att in latc.ncattrs():
                self.atts['latc'][an_att] = getattr(latc,an_att)
            self.data['latc'] = lonc[y1:y2:step,x1:x2:step]
            
        for var in extra_2dvars:
            nc_var = self.Dataset.variables[var]
            nc_var.set_auto_maskandscale(False)
            self.atts[var] = {}
            for an_att in nc_var.ncattrs():
                self.atts[var][an_att] = getattr(nc_var,an_att)
            self.data[var] = nc_var[t1:t2:ts,y1:y2:step,x1:x2:step]
        
    
    def make_vel_mask(self):
        try:
            fill_val = self.atts['u']['_FillValue']
        except KeyError:
            fill_val = self.atts['u']['missing_value']
        if len(self.data['u'].shape) == 4:
            u0 = self.data['u'][0,0,:,:]
        else:
            u0 = self.data['u'][0,:,:]
        self.grid['mask'] = (u0==fill_val).choose(1,0)
        
    def write_nc(self,ofn,is3d=False,gui_gnome=False,extra_2dvars=[],grid_only=False):
        """
        Write GNOME compatible netCDF file (netCDF3)
        * velocities are on center points and rotated to north/east
        * lat/lon can EITHER have the same dimensions as u/v or be one larger
          in both x/y dimensions (i.e. it defines the stencil) - GNOME will
          treat these two cases differently
        
        """
        
        nc = Dataset(ofn,'w',format='NETCDF3_CLASSIC')
        
        # Global Attributes
        setattr(nc,'grid_type','curvilinear')
        
        lon_key = 'lon'; lat_key = 'lat'
        if not grid_only:
            # test u/v dimensions
            if self.data['u'].shape != self.data['v'].shape:
                raise Exception('u/v dimensions differ')
            
            # determine if its a subset in time
            t_key = 'time'
            try:
                if self.data['u'].shape[0] == len(self.data['time_ss']):
                    t_key = 'time_ss'
            except KeyError: #TODO -- if it has key time_ss but it doesn't match that case is not caught
                if self.data['u'].shape[0] != len(self.data['time']):
                    raise Exception('Dimensions of u/v do not match time variable')
                
            # determine if its a subset of the grid
            try:
                lon_ss_shape = self.data['lon_ss'].shape
                lon_ss_shape_red = (lon_ss_shape[0]-1,lon_ss_shape[1]-1)
                if self.data['u'].shape[-2:] == lon_ss_shape or \
                    self.data['u'].shape[-2:] == lon_ss_shape_red:
                    lon_key = 'lon_ss'; lat_key = 'lat_ss'
            except KeyError:
                lon_shape = self.data['lon'].shape
                lon_shape_red = (lon_shape[0]-1,lon_shape[1]-1)
                if self.data['u'].shape[-2:] != lon_shape and \
                    self.data['u'].shape[-2:] != lon_shape_red:
                    raise Exception('Dimensions of u/v do not match grid variables')
    
        x = self.data[lon_key].shape[1]
        y = self.data[lat_key].shape[0]
        
        if grid_only:
            xc = x
            yc = y
            center_only = True
        else:
            if self.data[lon_key].shape == self.data['u'].shape[-2:]:
                #lat/lon are centerpoints
                center_only = True
                xc = x
                yc = y
            else:
                center_only = False
                xc = x-1
                yc = y-1
        
        # add Dimensions
        if not center_only:
            nc.createDimension('x',x)
            nc.createDimension('y',y)
        nc.createDimension('xc',xc)
        nc.createDimension('yc',yc)
        if is3d:
            zdim = 'sigma'
            zvar = 'sigma'
            try:            
                if self.atts['z']['standard_name'] == 'depth':
                    zdim= 'levels'
                    zvar = 'depth_levels'
            except KeyError:
                pass       
            nc.createDimension(zdim,len(self.data['z']))        
        nc.createDimension('time',None)

        if gui_gnome:
            nc_lonc = nc.createVariable('lon','f4',('yc','xc'))
            nc_latc = nc.createVariable('lat','f4',('yc','xc'))
        else:
            nc_lonc = nc.createVariable('lonc','f4',('yc','xc'))
            nc_latc = nc.createVariable('latc','f4',('yc','xc'))
        if not center_only:
            nc_lon = nc.createVariable('lon','f4',('y','x'))
            nc_lat = nc.createVariable('lat','f4',('y','x'))
            setattr(nc_lon,'standard_name','longitude')
            setattr(nc_lat,'standard_name','latitude')
            setattr(nc_lon,'units','degrees_east')
            setattr(nc_lat,'units','degrees_north')
        setattr(nc_lonc,'standard_name','longitude of center points')
        setattr(nc_latc,'standard_name','latitude of center points')
        setattr(nc_lonc,'units','degrees_east')
        setattr(nc_latc,'units','degrees_north')
        
        if not grid_only:

            nc_time = nc.createVariable('time','f4',('time',))
        
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
                    
            if 'wind' in self.atts:
                nc_u = nc.createVariable('air_u','f4',('time','yc','xc'), \
                    fill_value=ufill)
                nc_v = nc.createVariable('air_v','f4',('time','yc','xc'), \
                    fill_value=vfill)
            elif is3d:
                nc_u = nc.createVariable('water_u','f4',('time',zdim,'yc','xc'), \
                    fill_value=ufill)
                nc_v = nc.createVariable('water_v','f4',('time',zdim,'yc','xc'), \
                    fill_value=vfill)
            else:
                nc_u = nc.createVariable('water_u','f4',('time','yc','xc'), \
                    fill_value=ufill)
                nc_v = nc.createVariable('water_v','f4',('time','yc','xc'), \
                    fill_value=vfill)
        
        lon = self.data[lon_key]
        if self.dlx:
            lon = (lon < 0).choose(lon,lon+360)
        try:
            lonc = self.data['lonc']
            if self.dlx:
                lonc = (lonc < 0).choose(lonc,lonc+360)
        except KeyError:
            pass
        
         # add data
        if center_only:
            nc_lonc[:] = lon
            nc_latc[:] = self.data[lat_key]
        else:
            nc_lon[:] = lon
            nc_lat[:] = self.data[lat_key]
            nc_lonc[:] = lonc
            nc_latc[:] = self.data['latc']
        
        if not grid_only:
        #!!!!!!!!!!!!Add 3d

            if len(self.data['u'].shape) == 2:
                nc_time[0] = self.data[t_key]
                nc_u[0,:] = self.data['u']
                nc_v[0,:] = self.data['v']
            else:
                nc_time[:] = self.data[t_key]
                nc_u[:] = self.data['u']
                nc_v[:] = self.data['v']
        
            if is3d:
                nc_z = nc.createVariable(zvar,'f4',(zdim))
                nc_z[:] = self.data['z']
                nc_depth = nc.createVariable('depth','f4',('yc','xc'))
                nc_depth[:] = self.grid['depth']
                
            if 'mask' in self.grid:
                nc_mask = nc.createVariable('mask','f4',('yc','xc'))
                nc_mask[:] = self.grid['mask']
                setattr(nc_mask,'standard_name','mask on center points')
                setattr(nc_mask,'coordinates',u'latc lonc')
    
            # add variable attributes from 'atts' (nested dict object)
            for key,val in self.atts['time'].items():
                if not key.startswith('_'):
                    setattr(nc_time,key,val)
                
            self.atts['u']['coordinates'] = u'latc lonc'
            for key,val in self.atts['u'].items():
                if not key.startswith('_'):
                    setattr(nc_u,key,val)
            setattr(nc_u,'time','time')
                    
            self.atts['v']['coordinates'] = u'latc lonc'
            for key,val in self.atts['v'].items():

                if not key.startswith('_'):
                    setattr(nc_v,key,val)
            setattr(nc_v,'time','time')
                    
            for var in extra_2dvars:
                nc_var = nc.createVariable(var,'f4',('time','yc','xc'))
                nc_var[:] = self.data[var]
                setattr(nc_var,'coordinates',u'time latc lonc')
                for key,val in self.atts[var].items():
                    if not key.startswith('_'):
                        setattr(nc_var,key,val)
            
    
        nc.close()
    
class roms(cgrid):
    """
    A class for dealing with curvilinear grid model output and converting to GNOME format
    Requires passing in a var_map dict so the 
    variable names can be customized for different models or datasets
            
   """
         
    def get_dimensions(self,var_map,get_time=True,get_xy=True):
        
        if get_time:
            self.time = self.Dataset.variables[var_map['time']]
            self.atts['time'] = {}
            for an_att in self.time.ncattrs():
                self.atts['time'][an_att] = getattr(self.time,an_att) 
            self.data['time'] = self.time[:]
        
        if get_xy:
            #Load or create P grid lat/lon (sometimes not included in ROMS output)
            try:
                lon_psi = self.Dataset.variables['lon_psi']
                self.atts['lon_psi'] ={}            
                for an_att in lon_psi.ncattrs():
                    self.atts['lon_psi'][an_att] = getattr(lon_psi,an_att) 
                self.data['lon_psi'] = lon_psi[:]
                lat_psi = self.Dataset.variables['lat_psi']
                self.atts['lat_psi'] = {}
                for an_att in lat_psi.ncattrs():
                    self.atts['lat_psi'][an_att] = getattr(lat_psi,an_att) 
                self.data['lat_psi'] = lat_psi[:]
            except KeyError:
                print('Using rho grid to create P grid')
                lon_rho = self.Dataset.variables['lon_rho'][:]
                lat_rho = self.Dataset.variables['lat_rho'][:]
                self.data['lon_psi'] = (lon_rho[0:-1,0:-1]+lon_rho[1:,1:])*0.5
                self.atts['lon_psi'] = {'long_name':'longitude of PSI-points'}
                self.data['lat_psi'] = (lat_rho[0:-1,0:-1]+lat_rho[1:,1:])*0.5
                self.atts['lat_psi'] = {'long_name':'latitude of PSI-points'}
                
            self.data['lon'] = self.data['lon_psi']
            self.data['lat'] = self.data['lat_psi']

            
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

    def write_nc_native(self,ofn,is3d=False):
        """
      
        Write GNOME compatible netCDF file (netCDF3)
        Maintain u and v on sepearate grids
        
        %TODO: once this is implemented in GNOME, rename to "write_nc"
        to replace method in base class
        
        
        """
        nc = Dataset(ofn,'w',format='NETCDF3_CLASSIC')
        
        # Global Attributes
        setattr(nc,'grid_type','curvilinear')
        setattr(nc,'grid_case','arakawa-c')
    
        x = self.data['lon_psi'].shape[1]
        y = self.data['lat_psi'].shape[0]

        xc = x+1 #rho on center points plus outside boundary
        yc = y+1
        
        # determine if its a subset in time
        t_key = 'time'
        try:
            if self.data['u'].shape[0] == len(self.data['time_ss']):
                t_key = 'time_ss'
        except KeyError:
            if self.data['u'].shape[0] != len(self.data['time']):
                raise Exception('Dimensions of u/v do not match time variable')
        

        # add Dimensions
        nc.createDimension('x',x)
        nc.createDimension('y',y)
        nc.createDimension('xc',xc)
        nc.createDimension('yc',yc)
        nc.createDimension('ocean_time',None)
    
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
        nc_time = nc.createVariable('ocean_time','f4',('ocean_time',))
        nc_lonp = nc.createVariable('lon_psi','f4',('y','x'))
        nc_latp = nc.createVariable('lat_psi','f4',('y','x'))
        nc_lonr = nc.createVariable('lon_rho','f4',('yc','xc'))
        nc_latr = nc.createVariable('lat_rho','f4',('yc','xc'))
        nc_lonu = nc.createVariable('lon_u','f4',('yc','x'))
        nc_latu = nc.createVariable('lat_u','f4',('yc','x'))
        nc_lonv = nc.createVariable('lon_v','f4',('y','xc'))
        nc_latv = nc.createVariable('lat_v','f4',('y','xc'))
        nc_angle = nc.createVariable('angle','f4',('yc','xc'))

        if is3d:
            pass
        else:
            nc_u = nc.createVariable('u','f4',('ocean_time','yc','x'), \
                fill_value=ufill)
            nc_v = nc.createVariable('v','f4',('ocean_time','y','xc'), \
                fill_value=vfill)
        
        # add data
        nc_lonp[:] = self.data['lon_psi']
        nc_latp[:] = self.data['lat_psi']
        nc_lonr[:] = self.data['lon_rho']
        nc_latr[:] = self.data['lat_rho']
        nc_lonu[:] = self.data['lon_u']
        nc_latu[:] = self.data['lat_u']
        nc_lonv[:] = self.data['lon_v']
        nc_latv[:] = self.data['lat_v']
        nc_angle[:] = self.grid['angle']
        
        #!!!!!!!!!!!!Add 3d
        if len(self.data['u'].shape) == 2:
            nc_time[0] = self.data[t_key]
            nc_u[0,:] = self.data['u']
            nc_v[0,:] = self.data['v']
        else:
            nc_time[:] = self.data[t_key]
            nc_u[:] = self.data['u']
            nc_v[:] = self.data['v']
            
        if 'mask_rho' in self.grid:
            nc_mask = nc.createVariable('mask_rho','f4',('yc','xc'))
            nc_mask[:] = self.grid['mask_rho']

        # add variable attributes from 'atts' (nested dict object)
        for key,val in self.atts['time'].items():
            if not key.startswith('_'):
                setattr(nc_time,key,val)
        
        self.atts['u']['coordinates'] = "lon_u lat_u ocean_time"        
        for key,val in self.atts['u'].items():
            if not key.startswith('_'):
                setattr(nc_u,key,val)
        
        self.atts['v']['coordinates'] = "lon_v lat_v ocean_time"
        for key,val in self.atts['v'].items():
            if not key.startswith('_'):
                setattr(nc_v,key,val)
        
    
        nc.close()
        
          
