#!/usr/bin/env python

"""
setup.py for the libgoods package
"""

# This setup is suitable for "python setup.py develop".

from setuptools import setup

setup(
    name = "libgoods",
    version = "0.1.0",
    description = "utilities for pre and post processing data for GNOME",
    long_description=open("README.md").read(),
    packages = ["libgoods",],
    scripts = ["scripts/sample_script.py",
              ],
    author = "Amy MacFadyen",
    author_email = "amy.macfadyen@noaa.gov",
    url="http://www.response.restoration.noaa.gov/gnome",
    license = "LICENSE.txt",
    keywords = "graphics cython drawing",
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


