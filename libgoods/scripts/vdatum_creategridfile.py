# -*- coding: utf-8 -*-
from libgoods import tri_grid
import os

vdatum_dir = 'C:\\Users\\amy.macfadyen\\PyProjects\\GOODS\\trunk\\static\\ocean_models\\VDATUM'
files =     {
            #'EC2001':'EC2001_NOS_euv.nc',
            #'DelChes':'DEdelches01_adcirc54.nc',
            'FlSAB':'FLsab_adcirc54.nc',
            'GulfME':'MEgulfme01_adcirc54.nc',
            'NY':'NYsndbght02_adcirc54.nc',
            'PR':'PuertoRico01_adcirc54.nc',
            'TX':'TXcoast01_adcirc54.nc',
                }

for ds,fn in files.iteritems():         
    print ds
    ncfile = os.path.join(vdatum_dir,fn)
    var_map = {'latitude':'lat','longitude':'lon','nodes_surrounding_ele':'ele'}
    adcirc = tri_grid.ugrid(ncfile)
    adcirc.get_dimensions(var_map,get_time=False)

    adcirc.get_grid_topo(var_map)

    # find and order the boundary
    print 'Finding boundary'
    bnd = adcirc.find_bndry_segs()
    seg_types = [0] * len(bnd)
    print 'Ordering boundary'
    adcirc.order_boundary(bnd,seg_types)
    adcirc.write_bndry_file(os.path.join(vdatum_dir,ds + '.bry'))
