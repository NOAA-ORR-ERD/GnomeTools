#!/usr/bin/env python
import numpy as np
from netCDF4 import Dataset, MFDataset
from matplotlib import path
from libgoods import nctools
reload(nctools)

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
            
    def get_dimensions(self,var_map):
        
        self.time = self.Dataset.variables[var_map['time']]  
        self.atts['time'] = {}
        for an_att in self.time.ncattrs():
            self.atts['time'][an_att] = getattr(self.time,an_att) 
        self.data['time'] = self.time[:]
    
        lon = self.Dataset.variables[var_map['lon']]
        self.atts['lon'] = {}
        for an_att in lon.ncattrs():
            self.atts['lon'][an_att] = getattr(lon,an_att)
        self.data['lon'] = lon[:]
        
        lat = self.Dataset.variables[var_map['lat']]
        self.atts['lat'] = {}
        for an_att in lat.ncattrs():
            self.atts['lat'][an_att] = getattr(lat,an_att)
        self.data['lat'] = lat[:]        
    
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
        ii,jj = np.meshgrid(xrange(m),xrange(n))
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
            self.x = [[0,np.size(glat,1),1]]
        else: #do subset
            
            if dl == 0:
                [yvec,xvec] = np.where(np.logical_and(np.logical_and(glat>=sl,glat<=nl),np.logical_and(glon>=wl,glon<=el)))
            else:
                [yvec,xvec] = np.where(np.logical_and(np.logical_and(glat>=sl,glat<=nl),np.logical_or(glon>=wl,glon<=el)))
            
            if len(yvec) > 2 and len(xvec) > 2:
                y1 = min(yvec)
                y2 = max(yvec)+1
                x1 = min(xvec)
                x2 = max(xvec)+1
                self.y = [y1,y2,stride]
                self.x = [x1,x2,stride]           
            else:
                self.y = [0,np.size(glat,0),1]
                self.x = [[0,np.size(glat,1),1]]
               
    def get_grid_info(self,grid_vars=['mask'],yindex=None,xindex=None):
    
        if xindex is None and yindex is None:
            x1 = 0; x2 = self.data['lon'].shape[1]; step = 1
            y1 = 0; y2 = self.data['lon'].shape[0]
        else:
            y1 = yindex[0]; y2 = yindex[1]; step = yindex[2]
            x1 = xindex[0]; x2 = xindex[1]
        
        for var in grid_vars: 
        #these are 2D spatial vars (need to add vertical info) for 3d
            try:
                self.grid[var] = self.Dataset.variables[var][y1:y2:step,x1:x2:step]
            except KeyError:
                pass

    
    def get_data(self,var_map,tindex=None,yindex=None,xindex=None,zindex=0,is3d=False):
    
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
        
        u = self.Dataset.variables[var_map['u']]
        self.atts['u'] = {}
        for an_att in u.ncattrs():
            self.atts['u'][an_att] = getattr(u,an_att) 

        v = self.Dataset.variables[var_map['v']]
        self.atts['v'] = {}
        for an_att in v.ncattrs():
            self.atts['v'][an_att] = getattr(v,an_att) 
        
        self.data['u'] = u[t1:t2:ts,zindex,y1:y2:step,x1:x2:step]
        self.data['v'] = v[t1:t2:ts,zindex,y1:y2:step,x1:x2:step]
        
    def write_nc(self,ofn,is3d=False):
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
                
        lon_key = 'lon'; lat_key = 'lat'
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
        
        if self.data[lon_key].shape == self.data['u'].shape[-2:]:
            xc = x
            yc = y
        else:  
            xc = x-1
            yc = y-1

        # add Dimensions
        nc.createDimension('x',x)
        nc.createDimension('y',y)
        nc.createDimension('xc',xc)
        nc.createDimension('yc',yc)
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
        nc_time = nc.createVariable('time','f4',('time',))
        nc_lon = nc.createVariable('lon','f4',('y','x'))
        nc_lat = nc.createVariable('lat','f4',('y','x'))
        
        if self.atts.has_key('wind'):
            nc_u = nc.createVariable('air_u','f4',('time','yc','xc'), \
                fill_value=ufill)
            nc_v = nc.createVariable('air_v','f4',('time','yc','xc'), \
                fill_value=vfill)
        elif is3d:
            pass
        else:
            nc_u = nc.createVariable('water_u','f4',('time','yc','xc'), \
                fill_value=ufill)
            nc_v = nc.createVariable('water_v','f4',('time','yc','xc'), \
                fill_value=vfill)
        
         # add data
        nc_lon[:] = self.data[lon_key]
        nc_lat[:] = self.data[lat_key]
        
        #!!!!!!!!!!!!Add 3d
        if len(self.data['u'].shape) == 2:
            nc_time[0] = self.data[t_key]
            nc_u[0,:] = self.data['u']
            nc_v[0,:] = self.data['v']
        else:
            nc_time[:] = self.data[t_key]
            nc_u[:] = self.data['u']
            nc_v[:] = self.data['v']
            
        if self.grid.has_key('mask'):
            nc_mask = nc.createVariable('mask','f4',('yc','xc'))
            nc_mask[:] = self.grid['mask']

        # add variable attributes from 'atts' (nested dict object)
        for key,val in self.atts['time'].iteritems():
            if not key.startswith('_'):
                setattr(nc_time,key,val)
            
        for key,val in self.atts['u'].iteritems():
            if not key.startswith('_'):
                setattr(nc_u,key,val)
    
        for an_att in self.atts['v'].iteritems():
            if not key.startswith('_'):
                setattr(nc_v,key,val)
    
        nc.close()
    
class roms(cgrid):
    """
    A class for dealing with curvilinear grid model output and converting to GNOME format
    Requires passing in a var_map dict so the 
    variable names can be customized for different models or datasets
            
   """
         
    def get_dimensions(self,var_map):
        
        self.time = self.Dataset.variables[var_map['time']]
        self.atts['time'] = {}
        for an_att in self.time.ncattrs():
            self.atts['time'][an_att] = getattr(self.time,an_att) 
        self.data['time'] = self.time[:]
        
        #load lat/lon for rho, u, and v grids
        for var in ['lat_rho','lat_u','lat_v','lon_rho','lon_u','lon_v']:
            ds_var = self.Dataset.variables[var]
            self.atts[var] = {}
            for an_att in ds_var.ncattrs():
                self.atts[var][an_att] = getattr(ds_var,an_att)
            self.data[var] = ds_var[:]
        
        #Now load or create P grid lat/lon (sometimes not included in ROMS output)
        try:
            lon_psi = self.Dataset.variables['lon_psi']
            self.atts['lon_psi'] = {}            
            for an_att in lon_psi.ncattrs():
                self.atts['lon_psi'][an_att] = getattr(lon_psi,an_att) 
            self.data['lon_psi'] = lon_psi[:]
            lat_psi = self.Dataset.variables['lat_psi']
            self.atts['lat_psi'] = {}
            for an_att in lat_psi.ncattrs():
                self.atts['lat_psi'][an_att] = getattr(lat_psi,an_att) 
            self.data['lat_psi'] = lat_psi[:]
        except KeyError:
            self.data['lon_psi'] = (self.data['lon_rho'][0:-1,0:-1]+self.data['lon_rho'][1:,1:])*0.5
            self.atts['lon_psi'] = self.atts['lon_rho']
            self.atts['lon_psi']['long_name'] = 'longitude of PSI-points'
            self.data['lat_psi'] = (self.data['lat_rho'][0:-1,0:-1]+self.data['lat_rho'][1:,1:])*0.5
            self.atts['lat_psi'] = self.atts['lat_rho']
            self.atts['lat_psi']['long_name'] = 'latitude of PSI-points'
                

    def get_grid_info(self,yindex=None,xindex=None,is3d=False):
        
        if xindex is None and yindex is None:
            x1 = 0; x2 = self.data['lon_psi'].shape[1]
            y1 = 0; y2 = self.data['lon_psi'].shape[0]
        else:
            y1 = yindex[0]; y2 = yindex[1]
            x1 = xindex[0]; x2 = xindex[1]
        

        self.grid['mask'] = self.Dataset.variables['mask_rho'][y1:y2+1,x1:x2+1] 
        self.grid['angle'] = self.Dataset.variables['angle'][y1:y2+1,x1:x2+1] 
        
        if is3d:
            self.grid['h'] = self.Dataset.variables['h'][y1:y2+1,x1:x2+1] 
            self.grid['hc'] = self.Dataset.variables['hc'][:]
            self.grid['Cs_r'] = self.Dataset.variables['Cs_r'][:]
            self.grid['sc_r'] = self.Dataset.variables['s_rho'][:]
    
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
            self.data['lon_psi_ss'] = self.data['lon_psi'][y1:y2+1:step,x1:x2+1:step]
            self.data['lat_psi_ss'] = self.data['lat_psi'][y1:y2+1:step,x1:x2+1:step]
        
        u = self.Dataset.variables['u']
        self.atts['u'] = {}
        for an_att in u.ncattrs():
            self.atts['u'][an_att] = getattr(u,an_att) 
        v = self.Dataset.variables['v']
        self.atts['v'] = {}
        for an_att in v.ncattrs():
            self.atts['v'][an_att] = getattr(v,an_att) 

        if is3d: 
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
        else:
            self.data['u'] = u_on_upts
            self.data['v'] = v_on_vpts
            
    
    def interp_and_rotate(self,u,v,is3d=False):
        '''Calculate u/v on rho points -- we lose exterior most u/v values
        Then rotate to north/east
        '''
        self.grid['angle'] = self.grid['angle'][1:-1,1:-1]
        self.grid['mask'] = self.grid['mask'][1:-1,1:-1]
        cosa = (np.cos(self.grid['angle']) * self.grid['mask'])
        sina = (np.sin(self.grid['angle']) * self.grid['mask'])
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
              
        if self.data.has_key('lon_psi_ss'): #subset
            self.data['lon_psi_ss'] = self.data['lon_psi_ss'][:-1,:-1]
            self.data['lat_psi_ss'] = self.data['lat_psi_ss'][:-1,:-1]   
        else:
            self.data['lon_psi'] = self.data['lon_psi'][:-1,:-1]
            self.data['lat_psi'] = self.data['lat_psi'][:-1,:-1]

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
        nc_time = nc.createVariable('time','f4',('time',))
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
            nc_u = nc.createVariable('water_u','f4',('time','yc','x'), \
                fill_value=ufill)
            nc_v = nc.createVariable('water_v','f4',('time','y','xc'), \
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
            
        if self.grid.has_key('mask'):
            nc_mask = nc.createVariable('mask_rho','f4',('yc','xc'))
            nc_mask[:] = self.grid['mask']

        # add variable attributes from 'atts' (nested dict object)
        for key,val in self.atts['time'].iteritems():
            if not key.startswith('_'):
                setattr(nc_time,key,val)
            
        for key,val in self.atts['u'].iteritems():
            if not key.startswith('_'):
                setattr(nc_u,key,val)
    
        for an_att in self.atts['v'].iteritems():
            if not key.startswith('_'):
                setattr(nc_v,key,val)
    
    
        nc.close()
        
          