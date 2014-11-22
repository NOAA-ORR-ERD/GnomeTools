#!/usr/bin/env python

"""
__init__.py for libgoods package
"""

__version__="0.1.1"

# find data_files
import os

data_files_dir = os.path.join(os.path.split(__file__)[0],"data_files")
if not os.path.exists(data_files_dir):
    os.mkdir(data_files_dir)

