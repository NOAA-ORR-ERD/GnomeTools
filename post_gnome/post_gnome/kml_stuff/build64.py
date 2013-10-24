#!/usr/bin/env python

"""
build64.py

genereates a text file with the base64encoded contentes of a file

"""

import sys, base64

infilename = sys.argv[1]


outfilename = infilename + ".b64"

data = file(infilename, 'r').read()

data = base64.b64encode(data)

file(outfilename,'w').write(data)

file('junk.png','w').write(base64.b64decode(data))

print data

