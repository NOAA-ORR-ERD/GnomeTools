from libgoods import curv_grid, data_files_dir
import os

'''
Sample script to retrieve data from Arakawa c-grid type model

'''

url = 'http://oos.soest.hawaii.edu/thredds/dodsC/hioos/roms_native/hiig/ROMS_Hawaii_Regional_Ocean_Model_Native_Grid_best.ncd'


var_map = { 'time':'time',
           }  
           
hiroms = curv_grid.roms(url)
hiroms.get_dimensions(var_map)
hiroms.data['lon'] = hiroms.data['lon_psi']
hiroms.data['lat'] = hiroms.data['lat_psi']

#Only download last five timesteps
ti=[len(hiroms.data['time'])-5,len(hiroms.data['time']),1]

#hiroms.subset([21.223,-158.387,21.647,-157.799],lat='lat_psi',lon='lon_psi')
hiroms.get_grid_info()
hiroms.get_data(var_map)

ofn = os.path.join(data_files_dir,'hiroms_example.nc')
#hiroms.data['lon_ss'] = hiroms.data['lon_psi_ss']
#hiroms.data['lat_ss'] = hiroms.data['lat_psi_ss']

hiroms.data['lonc'] = hiroms.data['lon_rho']
hiroms.data['latc'] = hiroms.data['lat_rho']
hiroms.grid['mask'] = hiroms.grid['mask_rho']

hiroms.write_nc(ofn,is3d=False)


# This was for testing different cases...
#if subset:
#    # this case interpolates to rho grid and reduces lon/lat size 
#    # (compatible with GNOME at present)
#    hiroms.subset([21.223,-158.387,21.647,-157.799],lat='lat_psi',lon='lon_psi')
#    hiroms.get_grid_info(xindex=hiroms.x,yindex=hiroms.y)
#    hiroms.get_data(var_map,tindex=ti,xindex=hiroms.x,yindex=hiroms.y)
#    hiroms.reduce_latlon_mesh_for_GNOME()
#    ofn = os.path.join(data_files_dir,'hiroms_ss_rho_reduced.nc')
##    hiroms.data['lon_ss'] = hiroms.data['lon_psi_ss']
##    hiroms.data['lat_ss'] = hiroms.data['lat_psi_ss']
#    hiroms.write_nc(ofn,is3d=False)
#else:
#    #u/v interpolated to rho grid, u/v and lat/lon the same size (works in current GNOME)
#    print 'interp and reduce'
#    hiroms.get_grid_info()
#    hiroms.get_data(var_map,tindex=ti)
#    hiroms.get_grid_info()
#    hiroms.reduce_latlon_mesh_for_GNOME()
#    ofn = os.path.join(data_files_dir,'hiroms_rho_reduced.nc')
#    hiroms.data['lon'] = hiroms.data['lon_psi']
#    hiroms.data['lat'] = hiroms.data['lat_psi']
#    hiroms.write_nc(ofn,is3d=False)
#    
#    #u/v interpolated to rho grid, lat/lon on psi grid larger than u/v
#    print 'interp only'
#    hiroms.get_dimensions(var_map)
#    hiroms.get_grid_info()
#    hiroms.get_data(var_map,tindex=ti)
#    ofn = os.path.join(data_files_dir,'hiroms_rho.nc')
#    hiroms.data['lon'] = hiroms.data['lon_psi']
#    hiroms.data['lat'] = hiroms.data['lat_psi']
#    hiroms.write_nc(ofn,is3d=False)
#    
#    #u/v on native grids
#    print 'native'
#    hiroms.get_dimensions(var_map)
#    hiroms.get_grid_info()
#    hiroms.get_data(var_map,tindex=ti,interp=False)
#    ofn = os.path.join(data_files_dir,'hiroms_native.nc')
#    hiroms.write_nc_native(ofn,is3d=False)
#    
    
    
    
    
    
    
