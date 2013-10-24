#!/usr/bin/env python

"""
script to convert nc_particle-sytle netcdf files (as exported by GNOME)
to kmz for viewing in Google Earth

This is very simple -- it would be nice to expand it with many options!

"""
import sys

from post_gnome import nc_particles
from post_gnome.kml_stuff.write_kmz import write_kmz

if __name__ == "__main__":
    nc_filename = sys.argv[1]

    kmz_filename = nc_filename.rstrip('.nc') + ".kmz"

    print "processing:", nc_filename
    print "creating:", kmz_filename

    pf = nc_particles.Reader(nc_filename)
    data = pf.get_all_timesteps()

    write_kmz('kmz_filename', pf.times, data)


