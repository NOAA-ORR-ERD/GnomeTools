#!/usr/bin/env python

"""
This script processes a GNOME-generated nc_particles file to re-mark
elements that end up on the beach

Each element that has the "beached" status flag set at the end
of the run is re-set to have status code 7: off-maps

This makes them show up in MapRoom as a bright magenta :-)

It will only work if there have been no elements added or removed
during the run

This should be re-written to be more flexible, but it gets the
job done.

"""

# short status_codes(data) ;
#     status_codes:long_name = "particle status code" ;
#     status_codes:flag_values = 0LL, 2LL, 3LL, 7LL, 10LL, 12LL, 32LL ;
#     status_codes:flag_meanings = "0:not_released 2:in_water 3:on_land 7:off_maps 10:evaporated 12:to_be_removed 32:on_tideflat" ;

import sys
import shutil
import numpy as np
from netCDF4 import Dataset

try:
    infilename = sys.argv[1]
except IndexError:
    print("You need to pass in the netcdf file to be processed")


out_filename = infilename.rsplit(".", maxsplit=1)[0] + "_adjusted.nc"

shutil.copy(infilename, out_filename)

ds = Dataset(out_filename, mode='a')

particle_count = ds.variables['particle_count']
status_codes = ds.variables['status_codes'][:]
ids = ds.variables['id']
# check to make sure there aren't any lost elements:
assert np.alltrue(np.diff(particle_count) == 0), "some elements were lost"

# this is now assuming that the number of particles doesn't change.
num_particles = particle_count[0]

# first set the ones that started on land to not_released code.
# print(status_codes[:num_particles])
start_on_land = status_codes[num_particles : 2 * num_particles] == 3

# reshape the status codes:
status_codes = np.reshape(status_codes, (-1, num_particles))
ids = np.reshape(ids, (-1, num_particles))


for particle in range(num_particles):
    history = status_codes[:, particle]
    if status_codes[-1, particle] == 3:
        status_codes[:, particle] = 7

# reset the ones that started on land
status_codes[:2 * num_particles][..., start_on_land] = 0

# reshape it back:
status_codes.shape = (-1,)

ds.variables['status_codes'][:] = status_codes

ds.close()


