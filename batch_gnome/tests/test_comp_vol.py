#!/usr/bin/env python

"""
tests of the comp_volume code

designed to be run with nose

"""

import sys

sys.path.append("../lib")

import numpy as np
import nose

## now to test the volume computation on the grid:
#from tap_comp_volume import comp_volume
from cy_tap_comp_volume import comp_volume

import TAP_mod

# set up a grid:
grid = TAP_mod.Grid(min_long=-10, max_long=10, min_lat=-5, max_lat=5,num_lat=5,num_long=10)

def test_one_particle():
    positions = np.array(( (1.5, 2.5), ), np.float32)
    mass = np.array( (1.0,), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert mass_grid[5, 3] == mass[0] # should be full mass in the one grid box 
    assert mass_grid.sum() == mass.sum() # mass conserved

def test_two_particles():
    positions = np.array(( (1.5, 2.5), (-9.0, -4.0)), np.float32)
    mass = np.array( (1.0, 2.0), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert mass_grid[5, 3] == mass[0] # should be full mass in the one grid box 
    assert mass_grid[0, 0] == mass[1] # should be full mass in the one grid box 
    assert mass_grid.sum() == mass.sum() # mass conserved

def test_two_particles_one_cell():
    positions = np.array(( (1.5, 2.5), (1.0, 2.0)), np.float32)
    mass = np.array( (1.0, 2.0), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert mass_grid[5, 3] == mass[0] + mass[1] # should be full mass in the one grid box 
    assert mass_grid.sum() == mass.sum() # mass conserved

def test_off_grid_up():
    positions = np.array(( (5, 6), ), np.float32)
    mass = np.array( (1.0,), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert np.alltrue(mass_grid == 0) # nothing on the grid

def test_off_grid_down():
    positions = np.array(( (5, -6), ), np.float32)
    mass = np.array( (1.0,), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert np.alltrue(mass_grid == 0) # nothing on the grid

def test_off_grid_right():
    positions = np.array(( (11, 3), ), np.float32)
    mass = np.array( (1.0,), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert np.alltrue(mass_grid == 0) # nothing on the grid

def test_off_grid_left():
    positions = np.array(( (-11, 3), ), np.float32)
    mass = np.array( (1.0,), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert np.alltrue(mass_grid == 0) # nothing on the grid

def test_off_grid_upleft():
    positions = np.array(( (-11, 200), ), np.float32)
    mass = np.array( (1.0,), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert np.alltrue(mass_grid == 0) # nothing on the grid

def test_off_grid_downright():
    positions = np.array(( (1100, -200), ), np.float32)
    mass = np.array( (1.0,), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert np.alltrue(mass_grid == 0) # nothing on the grid

def test_on_edge_top(): # top edge is off grid
    positions = np.array(( (0.0, 5.0), ), np.float32)
    mass = np.array( (1.0,), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert np.alltrue(mass_grid == 0) # nothing on the grid
    
def test_on_edge_right(): # right edge is off grid
    positions = np.array(( (10.0, 0.0), ), np.float32)
    mass = np.array( (1.0,), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert np.alltrue(mass_grid == 0) # nothing on the grid

def test_upper_righ_corner(): #  off grid
    positions = np.array(( (10.0, 5.0), ), np.float32)
    mass = np.array( (1.0,), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert np.alltrue(mass_grid == 0) # nothing on the grid

def test_on_edge_bottom(): # bottom edge is on grid
    positions = np.array(( (1.0, -5.0), ), np.float32)
    mass = np.array( (1.0,), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert mass_grid[5,0] == mass[0]
    assert mass_grid.sum() == mass.sum() # conserve mass

def test_on_edge_left(): # left edge is on grid
    positions = np.array(( (-10.0, 0.0), ), np.float32)
    mass = np.array( (1.0,), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert mass_grid[0, 2] == mass[0]
    assert mass_grid.sum() == mass.sum() # conserve mass

def test_lower_left_corner(): #  on grid
    positions = np.array(( (-10.0, -5.0), ), np.float32)
    mass = np.array( (1.0,), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert mass_grid[0, 0] == mass[0]
    assert mass_grid.sum() == mass.sum() # conserve mass

## test the flags
        
def test_beached_flag():
    
    positions = np.array(( (1.5, 2.5),
                           (-9.0, -4.0)
                           ), np.float32)
    mass = np.array( (1.0, 2.0), np.float32 )
    flags = np.array( (0, 2) , dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert mass_grid[5, 3] == mass[0] # should be full mass in the one grid box 
    assert mass_grid[0, 0] == mass[1] # should be full mass in the one grid box 
    assert mass_grid.sum() == mass.sum() # mass conserved

def test_beached_flag2():
    # don't ignore it this time
    positions = np.array(( (1.5, 2.5),
                           (-9.0, -4.0)
                           ), np.float32)
    mass = np.array( (1.0, 2.0), np.float32 )
    flags = np.array( (0, 2) , dtype=np.uint8)
        
    mass_grid = comp_volume(positions, mass, flags, grid, flag_bitmask_to_ignore = 2)
    assert mass_grid[5, 3] == mass[0] # should be full mass in the one grid box 
    assert mass_grid[0, 0] == 0 # should be full mass in the one grid box 
    assert mass_grid.sum() == mass[0] # mass conserved (w/out beached oil)


def test_flag_notReleased():
    # don't ignore it this time
    positions = np.array(( (1.5, 2.5),
                           (-9.0, -4.0)
                           ), np.float32)
    mass = np.array( (1.0, 2.0), np.float32 )
    flags = np.array( (0, 1) , dtype=np.uint8) # 1 is the notReleased flag
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert mass_grid[5, 3] == mass[0] # should be full mass in the one grid box 
    assert mass_grid[0, 0] == 0 # should be full mass in the one grid box 
    assert mass_grid.sum() == mass[0] # mass conserved (w/out flagged oil)

def test_flag_offMap():
    # don't ignore it this time
    positions = np.array(( (1.5, 2.5),
                           (-9.0, -4.0)
                           ), np.float32)
    mass = np.array( (1.0, 2.0), np.float32 )
    flags = np.array( (0, 4) , dtype=np.uint8) # 4 is the offMap flag
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert mass_grid[5, 3] == mass[0] # should be full mass in the one grid box 
    assert mass_grid[0, 0] == 0 # should be full mass in the one grid box 
    assert mass_grid.sum() == mass[0] # mass conserved (w/out flagged oil)

def test_flag_evaporated():
    # don't ignore it this time
    positions = np.array(( (1.5, 2.5),
                           (-9.0, -4.0),
                           ), np.float32)
    mass = np.array( (1.0, 2.0), np.float32 )
    flags = np.array( (8, 0) , dtype=np.uint8) # 8 is the evaporated flag
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert mass_grid[5, 3] == 0 # should be full mass in the one grid box 
    assert mass_grid[0, 0] == mass[1] # should be full mass in the one grid box 
    assert mass_grid.sum() == mass[1] # mass conserved (w/out flagged oil)

def test_flag_noOnSurface():
    # don't ignore it this time
    positions = np.array(( (1.5, 2.5),
                           (-9.0, -4.0),
                           ), np.float32)
    mass = np.array( (1.0, 2.0), np.float32 )
    flags = np.array( (16, 0) , dtype=np.uint8) # 16 is the notOnSurface: flag
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert mass_grid[5, 3] == 0 # should be full mass in the one grid box 
    assert mass_grid[0, 0] == mass[1] # should be full mass in the one grid box 
    assert mass_grid.sum() == mass[1] # mass conserved (w/out flagged oil)

def test_multiple_LEs():
    # don't ignore it this time
    positions = np.array(( ( 1.5,  2.5),
                           (-9.0, -4.0),
                           ( 9.0,  4.0),
                           ( 0.0,  0.0),
                           ( 0.0, -4.0),
                           ( 5.0,  2.0),
                           (-6.0,  1.0)
                           
                           ), np.float32)
    mass = np.array( (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 5.0 ), np.float32 )
    flags = np.zeros_like(mass, dtype=np.uint8)
#    flags = np.array( (16, 0) , dtype=np.uint8) # 16 is the notOnSurface: flag
        
    mass_grid = comp_volume(positions, mass, flags, grid)
    assert mass_grid[5, 3] == mass[0] # should be full mass in the one grid box 
    assert mass_grid[0, 0] == mass[1] # should be full mass in the one grid box 
    assert mass_grid[9, 4] == mass[2] # should be full mass in the one grid box 
    assert mass_grid[5, 2] == mass[3] # should be full mass in the one grid box 
    assert mass_grid[5, 0] == mass[4] # should be full mass in the one grid box 
    assert mass_grid[7, 3] == mass[5] # should be full mass in the one grid box 
    assert mass_grid[2, 3] == mass[6] # should be full mass in the one grid box 
    assert mass_grid.sum() == mass.sum() # mass conserved (w/out flagged oil)

def test_multiple_flags():
    # don't ignore it this time
    positions = np.array(( ( 1.5,  2.5),
                           (-9.0, -4.0),
                           ( 9.0,  4.0),
                           ( 0.0,  0.0),
                           ( 0.0, -4.0),
                           ( 5.0,  2.0),
                           (-6.0,  1.0)
                           
                           ), np.float32)
    mass = np.array(  (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 5.0 ), np.float32 )
    flags = np.array( (  0,   1,   2,   4,   8,  16,   0 ), dtype=np.uint8)
    bitmask = 1 + 4 + 16
        
    mass_grid = comp_volume(positions, mass, flags, grid, flag_bitmask_to_ignore = bitmask)
    assert mass_grid[5, 3] == mass[0] # should be full mass in the one grid box 
    assert mass_grid[0, 0] == 0       # should be full mass in the one grid box 
    assert mass_grid[9, 4] == mass[2] # should be full mass in the one grid box 
    assert mass_grid[5, 2] == 0       # should be full mass in the one grid box 
    assert mass_grid[5, 0] == mass[4] # should be full mass in the one grid box 
    assert mass_grid[7, 3] == 0       # should be full mass in the one grid box 
    assert mass_grid[2, 3] == mass[6] # should be full mass in the one grid box 
    print mass_grid.sum()
    assert mass_grid.sum() == 14.0   # mass conserved (w/out flagged oil)

def test_multiple_flags_off_map():
    # don't ignore it this time
    positions = np.array(( ( 1.5,  2.5),
                           ( 100, 10  ),
                           ( 9.0,  4.0),
                           ( 0.0,  0.0),
                           ( 0.0, -4.0),
                           ( 5.0,  2.0),
                           (-6.0,  1.0)
                           
                           ), np.float32)
    mass = np.array(  (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 5.0 ), np.float32 )
    flags = np.array( (  0,   1,   2,   4,   8,  16,   0 ), dtype=np.uint8)
    bitmask = 1 + 4 + 16
        
    mass_grid = comp_volume(positions, mass, flags, grid, flag_bitmask_to_ignore = bitmask)
    assert mass_grid[5, 3] == mass[0] # should be full mass in the one grid box 
    assert mass_grid[0, 0] == 0       # should be full mass in the one grid box 
    assert mass_grid[9, 4] == mass[2] # should be full mass in the one grid box 
    assert mass_grid[5, 2] == 0       # should be full mass in the one grid box 
    assert mass_grid[5, 0] == mass[4] # should be full mass in the one grid box 
    assert mass_grid[7, 3] == 0       # should be full mass in the one grid box 
    assert mass_grid[2, 3] == mass[6] # should be full mass in the one grid box 
    print mass_grid.sum()
    assert mass_grid.sum() == 14.0   # mass conserved (w/out flagged oil)

