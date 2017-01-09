#!/usr/bin/env python

"""
Test code for nc_particles

Not very complete

Designed to be run with pytest

"""
# for py2/3 compatibility
from __future__ import absolute_import, division, print_function, unicode_literals

import datetime

import pytest
import netCDF4
from post_gnome import nc_particles


# test the writer...
def test_writer_with_ref_time():
    w = nc_particles.Writer('junk_file.nc', ref_time=datetime.datetime(2010, 2, 3, 0))
    del w
    nc = netCDF4.Dataset('junk_file.nc')
    units = nc.variables['time'].units
    nc.close()
    assert units == "seconds since 2010-02-03T00:00:00"


def test_write_timestep():
    """very simple version"""
    w = nc_particles.Writer('junk_file.nc',
                            ref_time=datetime.datetime(2010, 2, 3, 0),
                            nc_version=4)
    data = {"longitude": [43.2, 43.3, 43.4],
            "latitude": [31.0, 31.2, 31.3],
            "id": [1, 2, 3]
            }
    w.write_timestep(datetime.datetime(2010, 2, 3, 0), data)
    # and another (same data, but whatever..)
    w.write_timestep(datetime.datetime(2010, 2, 3, 0), data)
    del w
    ## read it back in


def test_write_timestep_wrong_size():
    """very simple version"""
    w = nc_particles.Writer('junk_file.nc',
                            nc_version=4)
    data = {"longitude": [43.2, 43.3, 43.4],
            "latitude": [31.0, 31.2, 31.3],
            "id": [1, 2]
            }
    with pytest.raises(ValueError):
        w.write_timestep(datetime.datetime(2010, 2, 3, 0), data)
    del w

# this one fails when run a part of the test suite
#  but runs fine on its own. -- bug in python-netcdf4?
# def test_multi_close():
#     w = nc_particles.Writer('junk_file.nc',
#                             nc_version=4)
#     w.close()
#     w.close()

if __name__ == "__main__":
    test_multi_close()
