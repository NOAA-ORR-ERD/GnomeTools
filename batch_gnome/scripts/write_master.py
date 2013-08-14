"""
writes a GNOME "master" file for multiple necdf files
"""


import sys, os, glob

try:
    path = sys.argv[1]
except IndexError:
    raise Exception("pass in the directory of all the files on the command line")


filenames = [os.path.split(f)[-1] for f in glob.glob(os.path.join(path, "*.nc"))]
filenames.sort()

print filenames

f = file(os.path.join(path, "master.txt"), 'w') 
f.write("NetCDF Files\n")
for name in filenames:
    f.write("[FILE] %s\n"%name)
f.close

