# for py2/3 compatibility
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import datetime

import pytest
import numpy as np
import netCDF4
from post_gnome import nc_particles

## test the Reader
def test_read_required():
    """
    Does it find the required variables and attributes
    Should be able to set up data_index
    """
    r = nc_particles.Reader('sample.nc')
    assert len(r.times) == 3
    assert np.array_equal(r.data_index, np.array([0, 3, 7, 9]))

def test_read_existing_dataset():
    nc = netCDF4.Dataset('sample.nc')
    r = nc_particles.Reader(nc)
    assert len(r.times) == 3

def test_str():
    r = nc_particles.Reader('sample.nc')
    print(r)
    r.close()
    assert True

## other tests fail (E   RuntimeError: NetCDF: Not a valid ID)
## if this test is here -- no idea why, but I think NetCDF4 isn't cleaning up after itself well
def test_read_variables():
    """
    does it find the data variables ?
    """
    r = nc_particles.Reader('sample.nc')
    # set(), because order doesn't matter
    varnames = set(r.variables)
    assert varnames == set(['latitude', 'depth', 'mass', 'id', 'longitude'])

def test_get_all_timesteps():
    r = nc_particles.Reader('sample.nc')
    data = r.get_all_timesteps(variables=['depth', 'mass', 'id'])
    print(data)
    assert 'depth' in data
    assert 'mass' in data
    assert 'id' in data
    ## better to check actual data, but what can you do?

def test_get_timestep():
    r = nc_particles.Reader('sample.nc')
    data = r.get_timestep(2, variables=['latitude', 'depth', 'mass', 'id', 'longitude'])

    # specific results from the sample file
    assert np.array_equal(data['longitude'], [-88.3, -88.1])
    assert np.array_equal(data['latitude'], [28.1, 28.0])
    assert np.array_equal(data['depth'], [0.0, 0.1])
    assert np.array_equal(data['mass'], [0.05, 0.06])
    assert np.array_equal(data['id'], [1, 3])

def test_get_individual_trajectory():
    r = nc_particles.Reader('sample.nc')

    path = r.get_individual_trajectory(1)
    assert np.array_equal(path['latitude'],  [28.0, 28.05, 28.1])
    assert np.array_equal(path['longitude'], [-88.1, -88.2, -88.3])

def test_get_units():
    r = nc_particles.Reader('sample.nc')

    assert r.get_units('depth') == 'meters'
    assert r.get_units('longitude') == 'degrees_east'

def test_get_attributes():
    r = nc_particles.Reader('sample.nc')

    assert  r.get_attributes('depth') == {'units' : "meters",
                                          'long_name' : "particle depth below sea surface",
                                          'standard_name' : "depth",
                                          'axis' : "z positive down",
                                          }
