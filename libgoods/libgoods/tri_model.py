#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from netCDF4 import Dataset, MFDataset, num2date
from libgoods import nctools, base

class tri(base.nc):
    """
    A class for dealing with triangular grid model output and converting to GNOME format
    Although I use variable names consistent with FVCOM, by passing in a var_map dict the 
    variable names can be customized for SELFE or ADCIRC
    
    Right now the attribute specifying whether the elements are orderd clockwise or counter
    clockwise needs to be manually added before writing to GNOME format (GNOME requres this, 
    but its not often specified in the model output)
        
    """
    
            
    def get_dimensions(self,var_map,get_time=True,get_xy=True):
        
        if get_time:
            self.time_varname = var_map['time']
            self.time = self.Dataset.variables[self.time_varname][:]
            self.time_units = self.Dataset[self.time_varname].units
            self.time_dimension = self.Dataset.variables[self.time_varname].dimensions
        
        if get_xy:
            if self.GridDataset is None:
                self.GridDataset = self.Dataset
                                
            self.lon = self.GridDataset.variables[var_map['longitude']][:]
            self.lat = self.GridDataset.variables[var_map['latitude']][:]               
            self.nv = self.Dataset.variables[var_map['nodes_surrounding_ele']]
            if self.nv.shape[0] > self.nv.shape[1]: #[num_faces,num_eles]
                self.nv = self.nv.transpose()
            
            try: #should there be a check here for nbe starting at 0?
                self.nbe = self.GridDataset.variables[var_map['eles_surrounding_ele']]
                     
                if self.nbe.shape[0] > self.nbe.shape[1]: #[num_faces,num_eles]
                    self.nbe = self.nbe[:].transpose()

            except KeyError:
                print('Building face-face connectivity')
                self.build_face_face_connectivity()
        
    def get_grid_topo(self,var_map):
        
        self.atts = dict()
        self.data = dict()
        nv = self.Dataset.variables[var_map['nodes_surrounding_ele']]
        self.atts['nv'] = dict()
        for an_att in nv.ncattrs():
            self.atts['nv'][an_att] = getattr(nv,an_att)
        nv = nv[:]
#        if nv.min() == 0:
#            nv = nv + 1
        if nv.shape[0] > nv.shape[1]: #nv should be nv[num_faces,num_eles]
            self.data['nv'] = nv.transpose()
        else:
            self.data['nv'] = nv
            
        try: #should there be a check here for nbe starting at 0?
            nbe = self.Dataset.variables[var_map['eles_surrounding_ele']]
            self.atts['nbe'] = dict()
            for an_att in nbe.ncattrs():
                self.atts['nbe'][an_att] = getattr(nbe,an_att)               
            if nbe.shape[0] > nbe.shape[1]: #nbe should be nbe[num_faces,num_eles]
                self.data['nbe'] = nbe[:].transpose()
            else:
                self.data['nbe'] = nbe[:]
        except KeyError:
            print('Building face-face connectivity')
            self.build_face_face_connectivity()
            
        try:
            #this returns ALL the edges (boundary and interior)
            edges = self.Dataset.variables[var_map['edge_node_connectivity']]
            self.atts['edges'] = dict()
            for an_att in edges.ncattrs():
                self.atts['edges'][an_att] = getattr(edges,an_att)  
            edges = edges[:]
#            if edges.min() == 0:
#                edges = edges + 1
            self.data['edges'] = edges
        except KeyError:
            print('No edge information')
            #TODO: call code to determine edges here
            pass
        
        for opt_var in ['depth','sigma','lonc','latc','a1u','a2u']:
            if opt_var in var_map:
                theVar = self.Dataset.variables[var_map[opt_var]]
                self.data[opt_var] = theVar[:]
                self.atts[opt_var] = dict()
                for an_att in theVar.ncattrs():
                    self.atts[opt_var][an_att] = getattr(theVar,an_att)  
            
#        if var_map.has_key('depth'):
#            self.data['depth'] = self.Dataset.variables[var_map['depth']][:]
#        if var_map.has_key('sigma'):
#            self.data['sigma'] = self.Dataset.variables[var_map['sigma']][:]

            
    def get_data(self,var_map,tindex=None,nindex=None,zindex=0):
    
        ''' 
        var_map is a dict mapping model variable names to common names
        tindex can be used to subset in time --> tindex = [start,stop,step]
        nindex is for subsetting grid -- list of nodes/elements to get u/v on
        zindex is for z layer (FVCOM has surface at zindex = 0, SELFE zindex = -1)
        '''
    
        if tindex:
            self.data['time_ss'] = self.data['time'][tindex[0]:tindex[1]:tindex[2]]
        else:
            tindex = [0,len(self.data['time']),1]
                    
        u = self.Dataset.variables[var_map['u_velocity']]
        u.set_auto_maskandscale(False)
        self.atts['u'] = dict()
        for an_att in u.ncattrs():
            self.atts['u'][an_att] = getattr(u,an_att)
        v = self.Dataset.variables[var_map['v_velocity']]
        v.set_auto_maskandscale(False)
        self.atts['v'] = dict()
        for an_att in v.ncattrs():
            self.atts['v'][an_att] = getattr(v,an_att)
        
        if nindex is None:
            if len(u.shape)==3:
                if zindex > 0:
                    self.data['u'] = u[tindex[0]:tindex[1]:tindex[2],0:zindex,:]
                    self.data['v'] = v[tindex[0]:tindex[1]:tindex[2],0:zindex,:]
                else:
                    self.data['u'] = u[tindex[0]:tindex[1]:tindex[2],zindex,:]
                    self.data['v'] = v[tindex[0]:tindex[1]:tindex[2],zindex,:]    
            elif len(u.shape)==2:
                self.data['u'] = u[tindex[0]:tindex[1]:tindex[2],:]
                self.data['v'] = v[tindex[0]:tindex[1]:tindex[2],:]
            else:
                print("Error:velocity is not 2 or 3 dimensional")
                raise
        else: #Spatial subset -- under development but *mostly* working (TODO: add 3D)
            if u.shape[-1] == max(self.data['nbe'].shape):
                # velocities on elements
                id = np.where(np.diff(self.eles_in_ss) > 1)[0]
                indexer= self.eles_in_ss
            else:
                # velocities on nodes
                id = np.where(np.diff(self.nodes_in_ss) > 1)[0]
                indexer = self.nodes_in_ss
            id2 = [-1]; id2.extend(id)

            #print 'Number of contiguous segments: ', len(id2)+1 
            firsttime = True
            for ii in range(len(id2)):
                #if not np.mod(ii,20):
                    #print ii,'of ', len(id2)
                sid = id2[ii]+1
                try:
                    fid = id2[ii+1]
                except IndexError:
                    fid = -1
                if  len(u.shape) == 3:
                    this_u = u[tindex[0]:tindex[1]:tindex[2],zindex,indexer[sid]-1:indexer[fid]]
                    this_v = v[tindex[0]:tindex[1]:tindex[2],zindex,indexer[sid]-1:indexer[fid]]
                elif len(u.shape)==2:
                    this_u = u[tindex[0]:tindex[1]:tindex[2],indexer[sid]-1:indexer[fid]]
                    this_v = v[tindex[0]:tindex[1]:tindex[2],indexer[sid]-1:indexer[fid]]
                else:
                    print("Error:velocity is not 2 or 3 dimensional")
                    raise
                    
                if firsttime:
                    self.data['u'] = this_u.copy()
                    self.data['v'] = this_v.copy()
                    firsttime = False
                    ax = this_u.ndim-1
                else:  
                    self.data['u'] = np.concatenate((self.data['u'],this_u),axis=ax)
                    self.data['v'] = np.concatenate((self.data['v'],this_v),axis=ax)
        
        # sometimes these are numpy masked arrays -- GNOME can't deal
        if type(self.data['u']) is np.ma.core.MaskedArray:
            self.data['u'] = self.data['u'].data
        if type(self.data['v']) is np.ma.core.MaskedArray:
            self.data['v'] = self.data['v'].data
        
        #self.atts['v']['fill_value'] = self.atts['v']['missing_value']         
        #self.atts['u']['fill_value'] = self.atts['u']['missing_value']       
        
    def find_bndry_segs(self,subset=False):
    
        if subset:
            nv = self.data['nv_ss']
            nbe = self.data['nbe_ss']
        else:
            nv = self.data['nv']
            nbe = self.data['nbe']
            
        bnd = []
        for i, face in enumerate(nbe.transpose()):
            for j, neighbor in enumerate(face):
                if neighbor == 0:
                    if j == 0:
                        bound = [nv[1,i], nv[2,i]]
                    elif j == 1:
                        bound = [nv[2,i], nv[0,i]]
                    elif j == 2:
                        bound = [nv[0,i], nv[1,i]]
#                    if j == num_vertices-1:
#                        bound = [nv[-1,i], nv[0,i]]
#                    else:
#                        bound = [nv[j,i], nv[j+1,i]]
                    bnd.append(bound)
        
        return bnd

    def order_boundary(self,b,seg_types=None):
    
        if seg_types is None:
            seg_types = [0]*len(b)
            
        obnd = dict()
        otype = dict()
        start_point = 0
        bnd_number = 0
        obnd[bnd_number] = [b.pop(start_point),]
        otype[bnd_number] = [seg_types.pop(start_point),]
        while len(b)>0:
            idx = [i for i,edge in enumerate(b) if edge[0]==obnd[bnd_number][-1][-1]]
            if len(idx) == 1:
                obnd[bnd_number].append(b.pop(idx[0]))
                otype[bnd_number].append(seg_types.pop(idx[0]))
            else:
                bnd_number = bnd_number + 1
                obnd[bnd_number] = [b.pop(0),]
                otype[bnd_number] = [seg_types.pop(0),]
                
        #format for GNOME ([node1,node2,bnd_num,bnd_type] - bnd_type=1 for open, 2 for closed)
        boundary = []
        for i, a_bnd in obnd.items():
            for j, seg in enumerate(a_bnd):
                #TODO -- need to make separate method for adding 1 to nv,boundary...
                boundary.append([seg[0], seg[1], i, otype[i][j]])
            
        self.data['bnd'] = boundary
        self.atts['bnd'] = {'long_name':'Boundary segment information required for GNOME model'}  
    
    def write_bndry_file(self,bnd_file):
        
        np.savetxt(bnd_file,self.data['bnd'],fmt='%i')
    
    def read_bndry_file(self,bnd_file): 
 
        bnd = []
        f = open(bnd_file,'r')
        for line in f:
            vals = [int(val) for val in line.split()]
            bnd.append(vals)
        f.close()
        return np.array(bnd),{'long_name':'Boundary segment information required for GNOME model'}
        
    def write_unstruc_grid(self,ofn):
        
        """
        
        Write GNOME compatible netCDF file (netCDF3) from unstructured (triangular) grid data
        
        """  
        nc = Dataset(ofn,'w',format='NETCDF3_CLASSIC')
        
        # Global Attributes
        setattr(nc,'grid_type','Triangular')
        
        # test u/v dimensions
        if self.data['u'].shape != self.data['v'].shape:
            print('u/v dimensions differ')
            raise
        
        # determine if its a subset in time
        t_key = 'time'
        try:
            if self.data['u'].shape[0] == len(self.data['time_ss']):
                t_key = 'time_ss'
        except KeyError:
            if self.data['u'].shape[0] != len(self.data['time']):
                print('Dimensions of u/v do not match time variable')
                raise
                
        lon_key = 'lon'; lat_key = 'lat'
        nv_key = 'nv'; nbe_key = 'nbe'
        lonc_key = 'lonc'; latc_key = 'latc'
        a1u_key = 'a1u'; a2u_key = 'a2u'
        # determine if its a subset of the grid
        try:
            if self.data['u'].shape[-1] == len(self.data['lon_ss']) or \
                self.data['u'].shape[-1] == self.data['nbe_ss'].shape[-1]:
                lon_key = 'lon_ss'; lat_key = 'lat_ss'
                nv_key = 'nv_ss'; nbe_key = 'nbe_ss'
                lonc_key = 'lonc_ss'; latc_key = 'latc_ss'
                a1u_key = 'a1u_ss'; a2u_key = 'a2u_ss'
        except KeyError:
            if self.data['u'].shape[-1] != len(self.data['lon']) and \
                self.data['u'].shape[-1] != self.data['nbe'].shape[-1]:
                print('Dimensions of u/v do not match grid variables')
                raise 
                
        # add Dimensions
        nc.createDimension('time',None)
        nc.createDimension('node',len(self.data[lon_key]))
        nc.createDimension('nele',np.shape(self.data[nbe_key])[1])
        nc.createDimension('nbnd',len(self.data['bnd']))
        nc.createDimension('nbi',4)
        nc.createDimension('three',3)
        if 'a1u' in self.data:
            nc.createDimension('four',4)
        if 'sigma' in self.data:
            nc.createDimension('sigma',len(self.data['sigma'])) 
        
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
        nc_lon = nc.createVariable('lon','f4',('node'))
        nc_lat = nc.createVariable('lat','f4',('node'))
        nc_nbe = nc.createVariable('nbe','int32',('three','nele'))
        nc_nv = nc.createVariable('nv','int32',('three','nele'))
        nc_bnd = nc.createVariable('bnd','int32',('nbnd','nbi'))
        if 'sigma' in self.data:
            nc_sigma = nc.createVariable('sigma','f4',('sigma','node'))
            #nc_sigma = nc.createVariable('sigma','f4',('sigma'))
            nc_depth = nc.createVariable('depth','f4',('node'))
        if 'lonc' in self.data:
            nc_lonc = nc.createVariable('lonc','f4',('nele'))
            nc_latc = nc.createVariable('latc','f4',('nele'))
        if 'a1u' in self.data:
            nc_a1u = nc.createVariable('a1u','f4',('four','nele'))
            nc_a2u = nc.createVariable('a2u','f4',('four','nele'))
            
        if self.data['u'].shape[-1] == len(self.data[lon_key]): #velocities on nodes
            if len(self.data['u'].shape) == 3:
                nc_u = nc.createVariable('u','f4',('time','sigma','node'),fill_value=ufill)
                nc_v = nc.createVariable('v','f4',('time','sigma','node'),fill_value=vfill)  
            else:
                nc_u = nc.createVariable('u','f4',('time','node'),fill_value=ufill)
                nc_v = nc.createVariable('v','f4',('time','node'),fill_value=vfill)
        else: #velocities on elements
            if len(self.data['u'].shape) == 3:
                nc_u = nc.createVariable('u','f4',('time','sigma','nele'),fill_value=ufill)
                nc_v = nc.createVariable('v','f4',('time','sigma','nele'),fill_value=vfill)
            else:   
                nc_u = nc.createVariable('u','f4',('time','nele'),fill_value=ufill)
                nc_v = nc.createVariable('v','f4',('time','nele'),fill_value=vfill)
                
        #adjust time if necessary
        ref_time = self.atts['time']['units'].split(' since ')[1]
        ref_year = int(ref_time[0:4])
        if ref_year < 1970:
            print('Adjusting reference time')
            self.data[t_key],self.atts['time']['units'] = \
                nctools.adjust_time(self.data[t_key],self.atts['time']['units'])
        
        #add data to netcdf file
        nc_time[:] = self.data[t_key]
        nc_lon[:] = self.data[lon_key]
        nc_lat[:] = self.data[lat_key]
        nc_u[:] = self.data['u']
        nc_v[:] = self.data['v']
        nc_bnd[:] = self.data['bnd']
        nc_nbe[:] = self.data[nbe_key]
        nc_nv[:] = self.data[nv_key]
        if 'sigma' in self.data:
            nc_sigma[:] = self.data['sigma'][:]
            nc_depth[:] = self.data['depth']
        if 'lonc' in self.data:
            lonc = self.data[lonc_key]
            nc_lonc[:] = (lonc > 180).choose(lonc,lonc-360)
            nc_latc[:] = self.data[latc_key]
        if 'a1u' in self.data:
            print(nc_a1u.shape)
            print(self.data[a1u_key].shape)
            nc_a1u[:] = self.data[a1u_key]
            nc_a2u[:] = self.data[a2u_key]
            
        #add variable attributes to netcdf file
        for an_att in self.atts['time'].items():
           setattr(nc_time,an_att[0],an_att[1])
        
        for an_att in self.atts['lon'].items():
            setattr(nc_lon,an_att[0],an_att[1])
        
        for an_att in self.atts['lat'].items():
            setattr(nc_lat,an_att[0],an_att[1])
        
        for an_att in self.atts['bnd'].items():
            setattr(nc_bnd,an_att[0],an_att[1])

        for an_att in self.atts['nbe'].items():
            setattr(nc_nbe,an_att[0],an_att[1])

        for an_att in self.atts['nv'].items():
            setattr(nc_nv,an_att[0],an_att[1])
        
        for an_att in self.atts['u'].items():
            if an_att[0] != '_FillValue':
                setattr(nc_u,an_att[0],an_att[1])
        
        for an_att in self.atts['v'].items():
            if an_att[0] != '_FillValue':
                setattr(nc_v,an_att[0],an_att[1])
        
        nc.close()
    
    def write_unstruc_grid_only(self,ofn):
        
        """
        
        Write netCDF file (netCDF3) of grid variables for unstructured (triangular) grid
        
        """  
        nc = Dataset(ofn,'w',format='NETCDF3_CLASSIC')
        
        # Global Attributes
        setattr(nc,'grid_type','Triangular')

        # add Dimensions
        nc.createDimension('node',len(self.data['lon']))
        nc.createDimension('nele',np.shape(self.data['nbe'])[1])
        nc.createDimension('nbnd',len(self.data['bnd']))
        nc.createDimension('nbi',4)
        nc.createDimension('three',3)
        #nc.createDimension('sigma',1) #coming soon?
        
        # create variables
        nc_lon = nc.createVariable('lon','f4',('node'))
        nc_lat = nc.createVariable('lat','f4',('node'))
        nc_nbe = nc.createVariable('nbe','int32',('three','nele'))
        nc_nv = nc.createVariable('nv','int32',('three','nele'))
        nc_bnd = nc.createVariable('bnd','int32',('nbnd','nbi'))
        
        #add data to netcdf file
        nc_lon[:] = self.data['lon']
        nc_lat[:] = self.data['lat']
        nc_bnd[:] = self.data['bnd']
        nc_nbe[:] = self.data['nbe']
        nc_nv[:] = self.data['nv']
        
        for an_att in self.atts['lon'].items():
            setattr(nc_lon,an_att[0],an_att[1])   
        for an_att in self.atts['lat'].items():
            setattr(nc_lat,an_att[0],an_att[1])       
        for an_att in self.atts['bnd'].items():
            setattr(nc_bnd,an_att[0],an_att[1])
        for an_att in self.atts['nbe'].items():
            setattr(nc_nbe,an_att[0],an_att[1])
        for an_att in self.atts['nv'].items():
            setattr(nc_nv,an_att[0],an_att[1])
        
        nc.close()
        
    def find_nodes_eles_in_ss(self,nl,sl,wl,el):
       
        print('Total number of eles: ', self.data['nbe'].shape[1])
        print('Total number of nodes: ', self.data['lon'].shape[0])
        
        #returns lists of eles and nodes, plus truncated and edited topology arrays (nbe, nv)
        subset_lat = np.nonzero(np.logical_and(self.data['lat']>=sl,self.data['lat']<=nl))[0]
        subset_lon = np.nonzero(np.logical_and(self.data['lon']>=wl,self.data['lon']<=el))[0]
        self.nodes_in_ss = np.intersect1d(subset_lat,subset_lon) + 1 #node numbering starts at 1 (subtract one for indexing lon/lat)
        
        self.eles_in_ss = []
        nv_ss = []; nbe_ss = []; 
        
        #determine which nodes are in subset boundary and elements with all nodes in ss
        print('Finding nodes and entire elements in ss')
        for ii, ele in enumerate(self.data['nv'].transpose()):
            #if all of the nodes are in subset domain keep -- otherwise get rid of it
            if (self.nodes_in_ss == ele[0]).any() and (self.nodes_in_ss == ele[1]).any() \
             and (self.nodes_in_ss == ele[2]).any():
                nv_ss.append(ele)
                nbe_ss.append(self.data['nbe'][:,ii])
                self.eles_in_ss.append(ii+1) #ele numbering starts at 1
            else:
                pass
            
        print('Number of eles in ss: ', len(self.eles_in_ss))
        print('Number of nodes in ss: ', len(self.nodes_in_ss))
              
        nbe_ss = np.array(nbe_ss).transpose()
        nv_ss = np.array(nv_ss).transpose()
        self.eles_in_ss = np.array(self.eles_in_ss)

        print('Remapping nodes and elements')
        #now remap nbe_ss, nv_ss to number of remaining nodes, and elements
        nv_ssr = nv_ss.copy()
        nbe_ssr = nbe_ss.copy()
        for ii in range(len(self.eles_in_ss)):
            for jj in range(3):
                nid = np.searchsorted(self.nodes_in_ss, nv_ss[jj,ii], side='left')
                nv_ssr[jj,ii] = nid+1
                if nbe_ss[jj,ii] != 0:
                    eid = np.searchsorted(self.eles_in_ss, nbe_ss[jj,ii], side='left')
                    if eid >= len(self.eles_in_ss) or self.eles_in_ss[eid] != nbe_ss[jj,ii]:
                        nbe_ssr[jj,ii] = 0
                    else:
                        nbe_ssr[jj,ii] = eid+1
                                   
        self.data['nbe_ss'] = nbe_ssr
        self.data['nv_ss'] = nv_ssr
        self.data['lon_ss'] = self.data['lon'][self.nodes_in_ss-1]
        self.data['lat_ss'] = self.data['lat'][self.nodes_in_ss-1]
        if 'lonc' in self.data:
            self.data['lonc_ss'] = self.data['lonc'][self.eles_in_ss-1]
            self.data['latc_ss'] = self.data['latc'][self.eles_in_ss-1]
        if 'a1u' in self.data:
            self.data['a1u_ss'] = self.data['a1u'][:,self.eles_in_ss-1]
            self.data['a2u_ss'] = self.data['a2u'][:,self.eles_in_ss-1]
        
    def find_subset_land_nodes(self,bndry_file):
      
        #find all the land segments and re-number to match new subset node numbers
        
        print('Remapping boundary segs to new subset node numbers')
        f = open(bndry_file)
        nodes = []
        for line in f:
            node1,node2,bnumber,flag = (int(l) for l in line.split())
            node1_id = np.where(self.nodes_in_ss == node1)[0]
            node2_id = np.where(self.nodes_in_ss == node2)[0]
            if len(node1_id) > 0 and len(node2_id) > 0 and flag == 0:
                nodes.extend([node1_id[0]+1,node2_id[0]+1])
             
        f.close()
        
        return np.unique(nodes)
    
    
    def build_face_face_connectivity(self):
        """
        builds the triangular connectivity array (nbe)
        essentially giving the neighbors of each triangle (face)
        """       
                
        num_vertices = 3
        num_faces = max(self.nv.shape)
        face_face = np.zeros( (num_faces, num_vertices), dtype=np.int32  )
        face_face += -1 # fill with -1

        # loop through all the triangles to find the matching edges:
        edges = {} # dict to store the edges in 
        for i, face in enumerate(self.nv.transpose()):
            # loop through edges:
            for j in range(num_vertices):
                edge = (face[j-1], face[j])
                if edge[0] > edge[1]: # sort the node numbers
                    edge = (edge[1], edge[0]) 
                # see if it is already in there
                prev_edge = edges.pop(edge, None)
                if prev_edge is not None:
                    face_num, edge_num = prev_edge
                    face_face[i,j] = face_num
                    face_face[face_num, edge_num] = i
                else:
                    edges[edge] = (i, j)
        nbe = (face_face + 1)
        nbe = nbe[:,[2,0,1]].transpose()
        self.nbe = nbe
        #self.atts['nbe'] = {'long_name': 'elements surrounding element'}
