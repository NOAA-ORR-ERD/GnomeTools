#!/usr/bin/env python

"""
simple memory tet script

This simple allocates more an more arrays (about 10MB each) 

you can then see if/when it barfs with a memory error
"""


array_size = 10 * 1048576 #(10MB)

import numpy as np

all_arrays = []
    
count = 0    
while  True:
    print "created %i arrays"%count
    print "approx %i MB"%(10 * count)
    all_arrays.append(np.ones((array_size,), np.uint8) )
    count +=1
    
