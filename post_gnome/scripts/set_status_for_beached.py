#!/usr/bin/env python

"""
This script processes a GNOME-generated nc_particles file to re-mark
elements that end up on the beach

Each element that has the "beached" status flag set at the end
of the run is re-set to have status code 7: off-maps

This makes them show up in MapRoom as a bright magenta :-)

It will only work if there hve been no elements added or removed
during the run

This should be re-written to be more flexible, but it gets the
job done.

"""


import sys
import shutil
import numpy as np
from netCDF4 import Dataset

infilename = sys.argv[1]
out_filename = infilename[:-3] + "_adjusted.nc"

shutil.copy(infilename, out_filename)

ds = Dataset(out_filename, mode='a')

particle_count = ds.variables['particle_count']
status_codes = ds.variables['status_codes']
ids = ds.variables['id']
# check to make sure there aren't any lost elements:
assert np.alltrue(np.diff(particle_count) == 0), "some elements were lost"

num_particles = particle_count[0]
# reshape the status codes:
status_codes = np.reshape(status_codes, (-1, num_particles))
ids = np.reshape(ids, (-1, num_particles))

for particle in range(num_particles):
    history = status_codes[:, particle]
    if status_codes[-1, particle] == 3:
        status_codes[:, particle] = 7

# reshape it back:
status_codes.shape = (-1,)
ds.variables['status_codes'][:] = status_codes

ds.close()

