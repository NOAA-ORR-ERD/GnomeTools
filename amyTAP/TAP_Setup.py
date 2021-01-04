#!/usr/bin/env python
"""
A way to import a particular setup by only changing one thing
"""

# getting the root dir from the command line:

import sys, os

## Add the "lib" dir to sys.path
##  fixme: this is really a bit of a kludge...

#print __file__
#lib_dir = os.path.join(os.path.split(__file__)[0], '../lib')
#print lib_dir
#sys.path.insert(0, lib_dir)

try:
    RootDir = sys.argv[1]
except IndexError:
    #raise Exception("You must pass in the RootDir on the command line")
    RootDir = os.getcwd()

if not os.path.exists(RootDir):
    raise Exception("RootDir: %s Doesn't exist"%RootDir)

sys.path.insert(0, RootDir)
setup = __import__('TAP_params')
setup.RootDir = RootDir


