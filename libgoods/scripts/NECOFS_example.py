from libgoods import utools
import datetime as dt
import netCDF4

data_url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/archives/necofs_mb'
bndry_file = 'MassB.bry'

#I use a mapping of FVCOM variable names to common names so that the class methods can also
#work with SELFE and ADCIRC which have different var names
#This seemed easier than finding them by CF long_names etc
var_map = { 'longitude':'lon', \
            'latitude':'lat', \
            'time':'time', \
            'u_velocity':'ua', \
            'v_velocity':'va', \
            'nodes_surrounding_ele':'nv',\
            'eles_surrounding_ele':'nbe',\
          }  

necofs = utools.ugrid(data_url)
time_var = necofs.Dataset.variables['time']
dtime = netCDF4.num2date(time_var[:],time_var.units)
start = dt.datetime(2011,4,13,6,0,0)
stop = dt.datetime(2011,4,14,6,0,0)
istart = netCDF4.date2index(start,time_var,select='nearest')
istop = netCDF4.date2index(stop,time_var,select='nearest')
print istart,istop

print 'Downloading data'
#necofs.get_data(var_map,tindex=[0,1,1]) #First time step only
#necofs.get_data(var_map) #All time steps in file
necofs.get_data(var_map,tindex=[istart,istop,1])

necofs.adjust_time() #GNOME can't handle pre 1980 start dates (in units)

necofs.get_bndry(bndry_file) 
#This file was pre-generated for this grid (somewhat manually as open water/land boundaries
#are not specified in the model output

necofs.atts['nbe']['order'] = 'cw'
#GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)

print 'Writing to GNOME file'
necofs.write_unstruc_grid('TestforGNOME.nc')