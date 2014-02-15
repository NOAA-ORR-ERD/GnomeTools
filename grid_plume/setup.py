#!/usr/bin/env python

from __future__ import division #Change the / operator to ensure true division throughout (Zelenke).
import os
import sys
import glob
import shutil

#from setuptools import setup, find_packages
from distutils.core import setup

from setuptools import find_packages
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy


if "cleanall" in "".join(sys.argv[1:]):
    rm_files = ['gridplume/filter.so',
                'gridplume/filter.pyd',
                'gridplume/filter.c',
                ]

    for files_ in rm_files:
        for file_ in glob.glob(files_):
            print "Deleting auto-generated files: {0}".format(file_)
            os.remove(file_)

    rm_dir = ['grid_plume.egg-info', 'build']
    for dir_ in rm_dir:
        print "Deleting auto-generated directory: {0}".format(dir_)
        try:
            shutil.rmtree(dir_)
        except OSError as e:
            print e

    # this is what distutils understands
    sys.argv[1] = 'clean'


# JS: Set ARCHFLAGS so we're not trying to build for powerpc
if sys.maxsize <= 2 ** 32:
    architecture = 'i386'
else:
    architecture = 'x86_64'


if sys.platform == 'darwin':
    # for the mac -- decide whether we are 32 bit build
    if architecture == 'i386':
        #Setting this should force only 32 bit intel build
        os.environ['ARCHFLAGS'] = "-arch i386"
    else:
        os.environ['ARCHFLAGS'] = "-arch x86_64"

ext_modules=[
    Extension("gridplume.filter",
              ["gridplume/filter.pyx"],
              include_dirs = [ numpy.get_include() ],
              )
    ]


print "running setup(...) function"

setup(
    name = "grid_plume",
    description = "Filters plume data into a grid",
    author = "NOAA",
    packages = find_packages(),
    ext_modules = cythonize(ext_modules),
    requires=['numpy'],
    scripts = ['scripts/grid_plume.py',
               'scripts/merge_grids.py',
               'scripts/swept_area.py',
               ]

#    options = dict(
#        build = dict(
#            compiler = "msvc",
#        ) if sys.platform.startswith( "win" ) else {},
#    )

)
