#!/usr/bin/env python

import os
import sys
import glob
import shutil

from distutils.core import setup

from distutils.extension import Extension
from Cython.Build import cythonize

import numpy

if "cleanall" in "".join(sys.argv[1:]):
    rm_files = [
        'gridplume/filter.so',
        'gridplume/filter.pyd',
        'gridplume/filter.c',
    ]

    for files_ in rm_files:
        for file_ in glob.glob(files_):
            print("Deleting auto-generated files: {0}".format(file_))
            os.remove(file_)

    rm_dir = ['grid_plume.egg-info', 'build']
    for dir_ in rm_dir:
        print("Deleting auto-generated directory: {0}".format(dir_))
        try:
            shutil.rmtree(dir_)
        except OSError as e:
            print(e)

    # this is what distutils understands
    sys.argv[1] = 'clean'


ext_modules = [
    Extension(
        "gridplume.filter",
        ["gridplume/filter.pyx"],
        include_dirs=[numpy.get_include()],
    )
]

setup(
    ext_modules=cythonize(ext_modules),
)
