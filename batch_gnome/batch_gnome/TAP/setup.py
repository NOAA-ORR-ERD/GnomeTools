#!/usr/bin/env python2.6

"""
The distutils setup.py script for the Hazmat package

"""

import numpy
from distutils.core import setup, Extension


setup(name="check_receptors",
      description="Module for TAP",
      author="Christopher Barker",
      author_email="Chris.Barker@noaa.gov",
#      packages=['TAP'],
      ext_modules=[Extension('check_receptors',
                             ['check_receptors.c']),
                   ],
      include_dirs = [ numpy.get_include() ],
      )


