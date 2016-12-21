# -*- coding: utf-8 -*-
from __future__ import print_function
from libgoods import tri_grid
from libgoods import data_files_dir
import os

ncfile = os.path.join(data_files_dir,'vdatum','vdatum_fl_sab_adcirc54.nc')
var_map = {'latitude':'lat','longitude':'lon','nodes_surrounding_ele':'ele'}
adcirc = tri_grid.ugrid(ncfile)
adcirc.get_dimensions(var_map,get_time=False)

adcirc.get_grid_topo(var_map)

# find and order the boundary
print('Finding boundary')
bnd = adcirc.find_bndry_segs()
seg_types = [0] * len(bnd)
print('Ordering boundary')
adcirc.order_boundary(bnd,seg_types)

adcirc.write_unstruc_grid_only(os.path.join(data_files_dir,'vdatum','fl_sab_grid.nc'))
