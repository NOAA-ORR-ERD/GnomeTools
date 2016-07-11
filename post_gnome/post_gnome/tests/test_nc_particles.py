#!/usr/bin/env python

"""
Test code for nc_particles

Not very complete

Designed to be run with pytest

"""
# for py2/3 compatibility
from __future__ import absolute_import, division, print_function, unicode_literals

import os

import pytest
import netCDF4
from post_gnome import nc_particles


def test_init():
    """
    Can the classes be intitialized?
    """
    w = nc_particles.Writer('junk_file.nc')
    del w
    print(os.getcwd())
    nc_particles.Reader('sample.nc')


def test_3_unlimited():
    with pytest.raises(ValueError):
        nc_particles.Writer('junk_file2.nc', nc_version=3)


def test_netcdf3():
    w = nc_particles.Writer('junk_file3.nc', num_timesteps=10, nc_version='3')
    w.close()
    nc = netCDF4.Dataset('junk_file3.nc')
    assert nc.file_format == 'NETCDF3_CLASSIC'


def test_netcdf4():
    w = nc_particles.Writer('junk_file4.nc', num_timesteps=10, nc_version=4)
    w.close()
    nc = netCDF4.Dataset('junk_file4.nc')
    assert nc.file_format == 'NETCDF4'


def test_netcdf_wrong():
    with pytest.raises(ValueError):
        nc_particles.Writer('junk_file.nc', nc_version='nc4')


def test_netcdf_wrong_num():
    with pytest.raises(ValueError):
        nc_particles.Writer('junk_file.nc', nc_version='5')


def test_multi_close():
    w = nc_particles.Writer('junk_file5.nc',
                            nc_version=4)
    w.close()
    w.close()
