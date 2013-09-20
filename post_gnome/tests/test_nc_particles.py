#!/usr/bin/env python

"""
test code for reading an example netcdf file for particle trajectories

designed to be run with pytest

"""

import pytest

import datetime

import numpy as np

import netCDF4

from post_gnome import nc_particles

test_filename = 'test_particles.nc'

## write the test file
import nc_particle_build_sample
nc_particle_build_sample.write_sample_file(test_filename)

## NOTE: nc_particle_build_sample provides teh core testing of writting.
def test_write_missmatch():
    """
    test to see that an exception is raised if a miss-match in array sizes occurs
    """
    writer = nc_particles.Writer("junk.nc",
                                 num_timesteps=1,
                                 ref_time=datetime.datetime.now(),
                                 )
    data = {'longitude': np.array( [ -88.0,
                                     -88.1,
                                    ], dtype=np.float64 ),
            'latitude': np.array( [ 28.0, 
                                    28.0,
                                    28.1,
                                  ], dtype=np.float64 ),
            }
    with pytest.raises(ValueError):
        writer.write_timestep(datetime.datetime.now(), data)
    writer.close()


def test_open_with_filename():
    pf = nc_particles.Reader(test_filename)
    assert True # it would have failed already...

def test_open_with_filename_fail():
    """
    a non-existant file
    """
    with pytest.raises(RuntimeError):
        pf = nc_particles.Reader('random_name')

def test_open_with_ncfile():
    nc_file = netCDF4.Dataset(test_filename)
    pf = nc_particles.Reader(nc_file)
    assert True # it would have failed already...


def test_read_first_timestep():
    pf = nc_particles.Reader(test_filename)
    particles = pf.get_timestep(0, ['latitude','longitude', 'mass', 'status_code'])
    print particles
    ## checking against data in "nc_particle_build_sample.py"
    assert np.array_equal(particles['latitude'], (28.0, 28.0, 28.1))
    assert np.array_equal(particles['longitude'], (-88.0, -88.1, -88.1))
    assert np.array_equal(particles['mass'], (0.01, 0.005, 0.007))
    assert np.array_equal(particles['status_code'], (1, 2, 3) )

def test_read_third_timestep():
    pf = nc_particles.Reader(test_filename)
    particles = pf.get_timestep(2, ['latitude','status_code'])
    print particles
    ## checking against data in "nc_particle_build_sample.py"
    assert np.array_equal(particles['latitude'], (28.0, 28.0))
    assert np.array_equal(particles['status_code'], (2, 3) )


def test_data_not_there():
    pf = nc_particles.Reader(test_filename)
    
    with pytest.raises(KeyError):
        particles = pf.get_timestep(0, ['random_variable_name'])

def test_individual_trajectory():
    pf = nc_particles.Reader(test_filename)
    traj = pf.get_individual_trajectory(particle_id=1)
    assert np.array_equal(traj['latitude'], ( 28.,  28.,  28.) )
    assert np.array_equal(traj['longitude'], ( -88.1, -88.1, -88.) )

def test_get_units():
    pf = nc_particles.Reader(test_filename)
    print pf.get_units('depth')
    assert pf.get_units('depth') == 'meters'
    assert pf.get_units('mass') == 'grams'



                          