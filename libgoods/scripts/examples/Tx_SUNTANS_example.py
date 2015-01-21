#!/usr/bin/env python
from libgoods import utools, nctools, data_files_dir
import os 
import numpy as np

'''
Sample script to retrieve data from unstructured grid netcdf "file" (can be
OPeNDAP url), generate necessary grid topology (boundary info), and write 
GNOME compatible output.

The boundary file is saved to the data files directory so it only needs 
to be generated once (unless you are subsetting the grid).

To process multiple files (urls) either
a) pass the filenames/urls in as a list -- this creates a netcdf4 MFDataset and is
a good option for not too many files (all output is written to one nc file for GNOME 
in this case)
b) add a file list loop -- in this case put it after the grid topo vars are loaded (as
this only has to be done once). See NGOFS_multifile_example.py

'''

# specify local file or opendap url 
data_file = os.path.join(data_files_dir,'GalvCoarse_2014_0000.nc')

# the utools class requires a mapping of specific model variable names (values)
# to common names (keys) so that the class methods can work with FVCOM, SELFE,
# and ADCIRC which have different variable names
# (This seemed easier than finding them by CF long_names etc)
#!!!!!!!!txsuntans output on server does not include eles_surrounding_ele info
#I have it saved as a netcdf file included in libgoods data_files directory
var_map = { 
            'time':'time',\
            'u_velocity':'uc', \
            'v_velocity':'vc', \
            'nodes_surrounding_ele':'cells',\
            'eles_surrounding_ele':'nbe',\
            'edge_node_connectivity':'edges',\
          }  

# class instantiation creates a netCDF Dataset object as an attribute
txsuntans = utools.ugrid(data_file)

# get longitude, latitude, and time variables
print 'Downloading data dimensions'
txsuntans.get_dimensions(var_map)

# UTM coordinates -- calculate lat/lon
x = txsuntans.Dataset.variables['xp'][:]
y = txsuntans.Dataset.variables['yp'][:]
lon = np.ones_like(x); lat = np.ones_like(x)
for ii in range(len(x)):
    lat[ii], lon[ii] = nctools.utmToLatLng(14,x[ii],y[ii])
txsuntans.data['lon'] = lon
txsuntans.data['lat'] = lat
txsuntans.atts['lon'] = {'long_name': 'longitude'}
txsuntans.atts['lat'] = {'long_name': 'latitude'}

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
txsuntans.get_grid_topo(var_map)

#TODO -- put this in utools
def order_boundary(b,seg_types):

    obnd = dict()
    otype = dict()
    bnd_number = 0
    obnd[bnd_number] = [b.pop(0),]
    otype[bnd_number] = [seg_types.pop(0),]
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
    for i, a_bnd in obnd.iteritems():
        for j, seg in enumerate(a_bnd):
            #TODO -- need to make separate method for adding 1 to nv,boundary...
            boundary.append([seg[0]+1, seg[1]+1, i, otype[i][j]])
        
    return np.array(boundary)
    

edge_types = txsuntans.Dataset.variables['mark'][:].tolist()
#"0 - computational; 1 - closed; 2 flux BC; 3 - stage BC; 4 - other BC; 5 - interproc; 6 - ghost
#Based on personal communication we use 1,2, and 3
bound_id, bound_type = zip(*[(i,x) for i,x in enumerate(edge_types) if x>0 and x<4])
bound_segs = txsuntans.data['edges'][bound_id,:]
bound_type_gnome = (np.array(bound_type)==1).choose(1,0)
#obnd,otype = order_boundary(bound_segs.tolist(),list(bound_type))

txsuntans.data['nv'] = txsuntans.data['nv'] + 1
txsuntans.data['bnd'] = order_boundary(bound_segs.tolist(),list(bound_type_gnome))
txsuntans.atts['bnd'] = {'long_name':'Boundary segment information required for GNOME model'} 


## GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
txsuntans.atts['nbe']['order'] = 'ccw'

'''
!!!!!!!!!!!!!All the stuff above here only has to be done once -- if you want to 
process multiple files, I'd put a loop here and just keep overwriting 
txsuntans.data['u'] and ['v'] and incrementing the output file name
Also need to change txsuntans.data['time'] appropriately
'''
# get the data
print 'Loading u/v'
txsuntans.get_data(var_map,zindex=0) #First time step only

  
print 'Writing to GNOME file'
txsuntans.write_unstruc_grid(os.path.join(data_files_dir, 'txsuntans_example.nc'))