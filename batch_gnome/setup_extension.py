#!/usr/bin/env cython

"""
simple setup.py to build cython extension
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("cy_tap_comp_volume",
                             ["cy_tap_comp_volume.pyx"],
                             include_dirs=[np.get_include()]
                             )]
)

