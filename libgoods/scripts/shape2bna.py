#!/usr/bin/env python

"""
Script to translate a Shape file to BNA for GNOME
Usage:
  shape2bna the_shapefile [outputfile]


"""


from __future__ import print_function
import sys
import fiona

USAGE = """ shape2bna.py -- convert shapefiles withland polygons to a BNA file for GNOME

shape2bna shapefile [outputfile]

shapefile -- the shape file to convert (the *.shp file in the collections)

outputfile [optional] -- The name you want for the BNA file. If omitted a file
                         will be gernated with the same name as the input,
                         and ".bna" as the extension.

The script assumes that there are a number of polygons in the file, each surrounding land.

It does not currently support MapBounds or SpillableArea

This script was tested with one particular example

The vertices must be in geo-coordiantes (Ideally WGS84)

It may not work on others.

"""

# get the filename from the input
try:
    infilename = sys.argv[1]
except IndexError:
    print("you must provide a file to convert")
    sys.exit()
try:
    outfilename = sys.argv[2]
except IndexError:
    # use the default
    outfilename = infilename[:-4] + ".bna"

print("reading:", infilename)

# start the BNA file
print("writing out:", outfilename)
with open(outfilename, 'w') as bna:

    for feature in fiona.open(infilename):
        
        coords = feature['geometry']['coordinates']
        
        for count, points in enumerate(coords):
            points = points[0]
            bna.write('"%i","1", %i\n' % (count, len(points)))
            for p in points:
                bna.write("%f, %f\n" % (p))

print("Wrote %i polygons to the BNA file" % (count + 1))
