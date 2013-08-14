#!/usr/bin/env python

"""
setup.py for the BatchGNOME package
"""

# This setup is suitable for "python setup.py develop".
# but it does snot install the scripts

from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np


# setup(
#     cmdclass = {'build_ext': build_ext},
#     ext_modules = [Extension("cy_tap_comp_volume",
#                              ["cy_tap_comp_volume.pyx"],
#                              include_dirs=[np.get_include()]
#                              )]
# )



setup(
    name = "batch_gnome",
    version = "0.2",
    packages = ["batch_gnome",],
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("batch_gnome.cy_tap_comp_volume",
                             ["batch_gnome/cy_tap_comp_volume.pyx"],
                             include_dirs=[np.get_include()],
                             )]

# should really install the scripts...
#    scripts = ["scripts/Verdat2Poly.py",
#               "scripts/Poly2Verdat.py",
#               "scripts/thin_bna.py"
#               ],
    )


