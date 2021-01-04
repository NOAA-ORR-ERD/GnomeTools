from libgoods import roms_model, data_files_dir
import os
reload(roms_model)
'''
Sample script to retrieve data from HIROMS aggreation

'''

url = 'http://oos.soest.hawaii.edu/thredds/dodsC/hioos/roms_native/hiig/ROMS_Hawaii_Regional_Ocean_Model_Native_Grid_best.ncd'


hiroms = roms_model.roms(url)
hiroms.get_dimensions(tvar='time') #its an aggregation so  this overrides "ocean_time" roms variable

hiroms.write_nc(ofn='hiroms_grid.nc',grid_only=True)

# Only download last five timesteps
# ti=[len(hiroms.time)-5,len(hiroms.time),1]

# hiroms.subset([21.223,-158.387,21.647,-157.799])


# ofn = os.path.join(data_files_dir,'hiroms_example.nc')

# hiroms.write_nc(ofn)


