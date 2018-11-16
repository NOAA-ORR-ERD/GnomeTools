#!/usr/bin/env python

"""
setup.py for the libgoods package
"""

import os
# This setup is suitable for "python setup.py develop".
from setuptools import setup


def get_version():
    """
    get version from __init__.py
    """
    with open(os.path.join("libgoods", "__init__.py")) as initfile:
        for line in initfile:
            line = line.strip()
            if line.startswith("__version__"):
                version = line.split("=")[1].strip(' "')
                return version


setup(
    name="libgoods",
    version=get_version(),
    description="utilities for pre-processing input data for GNOME",
    long_description=open("README.md").read(),
    packages=["libgoods"],
    scripts=["scripts/hycom2bna.py",
             "scripts/nc_time_shift.py",
             "scripts/shape2bna",
             ],
    author="Amy MacFadyen, Christopher H. Barker",
    author_email="amy.macfadyen@noaa.gov, chris.barker@noaa.gov",
    url="http://www.gnome.orr.noaa.gov",
    license="LICENSE.txt",
    keywords="GNOME GOODS ORR",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Utilities",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Meteorology/Oceanography",
        "License :: Public Domain",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2 :: Only",
        "Programming Language :: Python :: Implementation :: CPython",
                  ],
     )
