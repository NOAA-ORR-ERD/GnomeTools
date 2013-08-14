#!/usr/bin/env python

"""
tests of the netcdf file reading code

designed to be run with nose

"""

import sys

sys.path.append("../lib")

import numpy as np
import netCDF4
import nose

## assorted info about test file
data_file = "micro_test.nc"
num_LEs = 10


import nc_particles

def test_can_open():
    trajectory = nc_particles.nc_particle_file(data_file)

def test_num_LEs(): # should be ten at the end...
    trajectory = nc_particles.nc_particle_file(data_file)
    lat = trajectory.get_all_timesteps(variables=['latitude'])['latitude']
    assert len(lat[-1]) == 10
    
def test_flag():
    trajectory = nc_particles.nc_particle_file(data_file)
    flag = trajectory.get_all_timesteps(variables=['flag'])['flag']
    print flag
    #assert false
    

def test_all_vars():
    vars = ['longitude','latitude','mass','age','flag','id']
    trajectory = nc_particles.nc_particle_file(data_file)
    particle_count = trajectory.particle_count
    data = trajectory.get_all_timesteps(variables = vars)
    for key in vars:
        d = data[key]
        assert len(d) == trajectory.num_times
        for i in range(len(d)):
            assert len(d[i]) == particle_count[i]
            

         
