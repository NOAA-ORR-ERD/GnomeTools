#!/usr/bin/env python2.3

from distutils.core import setup, Extension

setup (name = "TAP_ext",
       version = "1.0",
       maintainer = "Chris Barker",
       maintainer_email = "Chris.Barker@noaa.gov",
       description = "Multimodule package version of extentions for TAP",

       ext_modules = [Extension('TAP_ext',sources=['TAP_ext.c'])]
)

