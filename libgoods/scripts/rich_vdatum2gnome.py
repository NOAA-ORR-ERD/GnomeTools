from __future__ import print_function
from libgoods import tri_grid
import numpy as np
from utide import _ut_constants_fname
from utide.utilities import loadbunch
import pytz
from datetime import datetime
from pandas import date_range

import iris
import warnings

from matplotlib.dates import date2num
from utide.harmonics import FUV

# vdatum file location and grid variable names
#nc_file = ('http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/tides/'
#              'vdatum_fl_sab/adcirc54.nc')
nc_file = 'http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/tides/DEdelches01_adcirc54.nc'
var_map = {'latitude':'lat','longitude':'lon','nodes_surrounding_ele':'ele','eles_surrounding_ele':''}

adcirc = tri_grid.ugrid(nc_file)
adcirc.get_dimensions(var_map,get_time=False)
adcirc.get_grid_topo(var_map)

# define a subset
#nl = 30.6; sl = 29
#wl = -82; el = -80.6
nl = 39.83; sl = 38.428
wl = -76.009; el = -74.224

adcirc.find_nodes_eles_in_ss(nl,sl,wl,el)

#date range
start = datetime.strptime('6-Jul-2017 00:00',
                          '%d-%b-%Y %H:%M').replace(tzinfo=pytz.utc)
stop = datetime.strptime('7-Jul-2017 00:00',  # '18-Sep-2015 18:00'
                         '%d-%b-%Y %H:%M').replace(tzinfo=pytz.utc)
dt = 1.0  # Hours.

#find and order grid boundary
bnd = adcirc.find_bndry_segs(subset=True)
print('Size of boundary: ', len(bnd))
seg_types = [0] * len(bnd)
adcirc.order_boundary(bnd,seg_types)



def parse_string(name):
    lista = [e.decode().strip() for e in name.tolist()]
    return ''.join(lista)

names = []
const = adcirc.Dataset.variables['tidenames'][:]
for name in const:
    names.append(parse_string(name)) #changed from name.data

con_info = loadbunch(_ut_constants_fname)['const']


k = 0
ind_nc, ind_ttide = [], []

const_name = [e.strip() for e in con_info['name'].tolist()]

consts = ['STEADY', 'M2', 'S2', 'N2', 'K1', 'O1', 'P1', 'M4', 'M6']
for name in consts:
    try:
        if name == 'STEADY':
            indx = const_name.index('Z0')
        else:
            indx = const_name.index(name)
        k += 1
        ind_ttide.append(indx)
        ind_nc.append(names.index(name))
    except ValueError:
        pass  # `const` not found.
 

glocals = date_range(start, stop, freq='1H').to_pydatetime()
ntimes = len(glocals)



lon = adcirc.data['lon'][:]
lat = adcirc.data['lat'][:]



inbox = np.logical_and(np.logical_and(lon >= wl,
                                      lon <= el),
                       np.logical_and(lat >= sl,
                                      lat <= nl))
lon = lon[inbox]
lat = lat[inbox]


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    cubes = iris.load_raw(nc_file)


ua = cubes.extract_strict('Eastward Water Velocity Amplitude')
up = cubes.extract_strict('Eastward Water Velocity Phase')
va = cubes.extract_strict('Northward Water Velocity Amplitude')
vp = cubes.extract_strict('Northward Water Velocity Phase')

uamp = ua.data[0, inbox, :][:, ind_nc]
vamp = va.data[0, inbox, :][:, ind_nc]
upha = up.data[0, inbox, :][:, ind_nc]
vpha = vp.data[0, inbox, :][:, ind_nc]


# this takes 12 minutes instead of Iris 14 seconds
#uamp = adcirc.Dataset.variables['u_amp'][0,adcirc.nodes_in_ss,:][:,ind_nc]
#vamp = adcirc.Dataset.variables['v_amp'][0,adcirc.nodes_in_ss,:][:,ind_nc]
#upha = adcirc.Dataset.variables['u_phase'][0,adcirc.nodes_in_ss,:][:,ind_nc]
#vpha = adcirc.Dataset.variables['v_phase'][0,adcirc.nodes_in_ss,:][:,ind_nc]


freq_nc = adcirc.Dataset.variables['tidefreqs'][:][ind_nc] 
freq_ttide = con_info['freq'][ind_ttide]
t_tide_names = np.array(const_name)[ind_ttide]
omega_ttide = 2*np.pi * freq_ttide  # Convert from radians/s to radians/hour.

omega = freq_nc * 3600

v, u, f = FUV(t=np.array([date2num(start)]), tref=np.array([0]),
              lind=np.array([ind_ttide]),
              lat=55, ngflgs=[0, 0, 0, 0])
              
# Convert phase in radians.
v, u, f = (np.squeeze(i) for i in (v, u, f))
v = v * 2 * np.pi
u = u * 2 * np.pi

thours = np.array([d.total_seconds() for d in
                   (glocals - glocals[0])]) / 60 / 60.
   
adcirc.data['u'] = np.ones((len(thours),len(adcirc.nodes_in_ss)),)   
adcirc.data['v'] = np.ones((len(thours),len(adcirc.nodes_in_ss)),)

k = 0
for k in range(len(thours)):          
    U = (f * uamp * np.cos(v + thours[k] * omega + u - upha * np.pi/180)).sum(axis=1)
    V = (f * vamp * np.cos(v + thours[k] * omega + u - vpha * np.pi/180)).sum(axis=1) 
    adcirc.data['u'][k,:] = U
    adcirc.data['v'][k,:] = V

adcirc.data['time'] = thours
adcirc.atts['time'] = {'units':'hours since ' + start.strftime('%Y-%m-%d %H:%M')}
adcirc.atts['u'] = {'units':'m/s','_FillValue':999}
adcirc.atts['v'] = {'units':'m/s','_FillValue':999}
adcirc.atts['nbe'] = {'order':'ccw'}
adcirc.write_unstruc_grid('test.nc')                           






