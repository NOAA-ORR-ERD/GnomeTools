from libgoods import curv_grid, data_files_dir
import os

'''
Sample script to retrieve data from Arakawa c-grid type model

'''
grid_file = 'grid_Arctic_1.nc'
data_file = 'arctic_avg2_0004.nc'
bbox = [65,170,75,240]
subset = False

var_map = { 'time':'ocean_time',
           }  
         
arcroms = curv_grid.roms(os.path.join(data_files_dir,grid_file))
arcroms.get_dimensions(var_map,get_time=False)

if subset:
    arcroms.subset(bbox,lat='lat_psi',lon='lon_psi')
    arcroms.get_grid_info(xindex=arcroms.x,yindex=arcroms.y)
    print arcroms.x, arcroms.y
else:
    arcroms.get_grid_info()



arcroms.update(os.path.join(data_files_dir,data_file))
arcroms.get_dimensions(var_map,get_xy=False)

if subset:
    arcroms.get_data(var_map,xindex=arcroms.x,yindex=arcroms.y)
else:
    arcroms.get_data(var_map)

ofn = os.path.join(data_files_dir,data_file[:-3] + '_gnome.nc')

'''
Right now we need to specify which are center lon/lat and which are stencil
This is because I am using the generic curv_grid cgrid.write_nc method
Eventually we will want a special format for ROMS preserving all the grids so this
is a bit of a kludge at this point...
Also, because we want to load into GUI gnome we have to reduce 
the size of lon/lat to match u/v
'''
if subset:
    arcroms.data['lon_ss'] = arcroms.data['lon_ss'][:-1,:-1]
    arcroms.data['lat_ss'] = arcroms.data['lat_ss'][:-1,:-1]
else:
    arcroms.data['lon'] = arcroms.data['lon_psi'][:-1,:-1]
    arcroms.data['lat'] = arcroms.data['lat_psi'][:-1,:-1]

    
arcroms.data['lonc'] = arcroms.data['lon_rho']
arcroms.data['latc'] = arcroms.data['lat_rho']
arcroms.grid['mask'] = arcroms.grid['mask_rho']

arcroms.write_nc(ofn,gui_gnome=True,is3d=False)



#    
    
    
    
    
    
    
