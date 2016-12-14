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


hiroms.write_nc(ofn,is3d=False)


