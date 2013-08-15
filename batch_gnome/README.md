# batch_gnome:

A set of scripts and utilities for running GNOME in batch mode

This is a set of python scripts an utilities that can be used to run the NOAA GNOME model in a batch mode, particularly to generate data for TAP. For more details about that, see the Office of Response and Restoration web site:

<http://response.restoration.noaa.gov/gnome>

<http://response.restoration.noaa.gov/tap>


# Requirements:

(NOTE: all this should work in Windows or OS-X -- it will work on Linux too, but GNOME 1.3 itself doesn't, so not much use.)

## Python

Python is a very high level dynamic lanuage that is seeing wide adoption in the automation and scientific software communities.

We've been testing on version 2.7 -- older versions should work, but 3.* will not. You can get python from:

<http://www.python.org>

### The numpy package:

<http://numpy.scipy.org/>

### The netcdf4 package:

<http://code.google.com/p/netcdf4-python/>

(Windows binaries available from: <http://www.lfd.uci.edu/~gohlke/pythonlibs/> )

### Compiled code:

For some compiled code, you'll need:

#### Cython:

<http://cython.org/>

And a compiler set up to compile numpy extensions to Python:

- The system gcc on OS-X (and Linux)
- MS Visual studio 2008 (I think) -- the "Express" version works fine, at least for 32 bit python. It can be a bit hard to find, but MS should still distribute it.

Binaries for Windows and/or OS-X may be distributed with this package.

# What's here

With this package, you will find the following directories:

`batch_gnome` -- a python package of assorted modules used by the other scripts

`scripts` -- The various scripts you'll want to run

`tests` -- Some test code.

`sample_run` -- A complete sample set up for doing TAP.

# Running the Scripts

## Preparation 

The batch_gnome package has to be bulit and installed. I like to install it in "develop' mode:
  
    python setup.py develop

This should build the TAP extension required,and install eveythogn in batch_gnome where Python can find it.

## Setting up your GNOME, etc.
 
Create a directory to put all of your setup in -- see the `sample_run` dir for an example -- you may want to make of copy of that, and edit from there.
 
In your setup dir, create/Edit a `Setup_TAP.py` file -- this holds all configuration information for your setup -- the sample is commented, hopefully you can figure it out!

## Running the scripts

The scripts have not been set as a proper python package, so you need to run them from the Scripts directory.

`BuildAll.py` is a script that runs all the other scripts, in order, to do the whole procedure. In practice, you may want to run the scripts one at a time, so you can look in BuildAll.py to see which ones are run, in which order.

Most of the scripts are run by passing in the path to your Setup Data:

    python BuildAll.py /Users/chris.barker/Temp/BatchGNOME/Tests/Sample_run

(on your system, of course that path will be different. It may or may not work with a relative path -- using an absolute path is safer.

Good luck!




