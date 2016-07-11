#!/usr/bin/env python

"""
test code for kmz files
"""

from post_gnome import kmz_particles

def test_init():
    kmz_particles.Writer()


def test_copy():
    kmz_particles.nc2kmz('sample.nc')

    raise False
