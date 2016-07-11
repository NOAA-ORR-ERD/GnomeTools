#!/usr/bin/env python

"""
test code for kmz files
"""

from post_gnome import kmz_particles


def test_init():
    kmz_particles.Writer('junk')


def test_copy():
    name = kmz_particles.nc2kmz('sample.nc')
    assert name == 'sample.kmz'

    ## fixme -- should put some real tests in here!
    assert True

# def test_copy_real_gnome():
#     """
#     see if converting to kml works for a real GNOME file
#     """
#     name = kmz_particles.nc2kmz('chesapeake_bay_example.nc')
#     assert name == 'chesapeake_bay_example.kmz'

#     ## fixme -- should put some real tests in here!
#     assert False

