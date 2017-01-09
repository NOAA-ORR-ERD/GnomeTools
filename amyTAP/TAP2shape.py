#!/usr/bin/env python

"""
script to process TAP output, and create a shapefile with the data

requires the py shapefile library:

http://pypi.python.org/pypi/pyshp

fixme: We should add polygon ID to this file -- see example from Jill
( the issue is that TAP uses indexing from 1, Arc shape files use
  indexing from 0, so explicitly setting an ID make sit clear )

"""
import sys
import shapefile

USAGE = """ TAP2shape.py
Script to convert the TAP tab-delimited output to a shapefile.

Usage:

TAP2shape.py site.txt_file [data_file] output_file

site.txt_file:  the path to the site.txt file delivered with TAP.
                It canbe found in the TAPDATA directory.

data_file: is the path to the data file exported from TAP

output_file: path of the desired output shape file (and dbf file, etc.)

If no data_file is provided, an empty shape file will be produced.

"""



def IsClockwise(poly):
    """
    returns True if the polygon is clockwise ordered, false if not
    
    expects a sequence of tuples, or something like it (Nx2 array for instance),
    of the points:
    
    [ (x1, y1), (x2, y2), (x3, y3), ...(xi, yi) ]
    
    See: http://paulbourke.net/geometry/clockwise/
    """
    
    total = poly[-1][0] * poly[0][1] - poly[0][0]*poly[-1][1] # last point to first point
    for i in xrange(len(poly)-1):
        total += poly[i][0] * poly[i+1][1] - poly[i+1][0]*poly[i][1]
        
    if total <= 0:
        return True
    else:
        return False


def WriteSites(filename, polys, data, field_name):
    w = shapefile.Writer(shapefile.POLYGON) # shape type 5 is a polygon
    w.autoBalance = True
    
    # add the polygon index field:
    w.field("polygon", 'N', 10, 0)
    if data is not None:
        # add the data field:
        print "adding the data field"
        w.field(field_name,'N', 10, 4)
    for i, poly in enumerate(polys):
        #print "processing polygon:", i
        if not IsClockwise(poly):
            ## need to reverse the order of the points:
            ##  shape files expect clockwise for outer polygon
            ##   TAP used ccw for the grid polygons
            poly.reverse()
        w.poly(parts=[poly,])
        if data is None:
            #print "adding only polygon id"
            w.record(i+1)
        else:
            #print "adding polygon id and data"
            w.record(i+1, `data[i]`)
    w.save(filename)


def ReadSites(filename):
    """
    reads the site polygons from the SITE.txt file
    
    filename is the path to the SITE.txt file
    """
    
    infile = file(filename, 'U')
    
    # scan for sites line:
    print "Reading SITE.TXT file for site polygons"
    while True:
        line  = infile.readline()
        if line.split()[1].upper() == "SITES":
            break
    num_sites = int(line.split()[0])
    # read the polygons
    polys = []
    for i in xrange(num_sites):
        header = infile.readline().strip().split(",")
        num_points = int(header[2])
        #read the points
        poly = []
        for j in xrange(num_points):
            point = infile.readline().strip().split(",")
            point = (float(point[0]), float(point[1]))
            poly.append(point)
        polys.append(poly)
    print "Read %i polygons"%len(polys)
    return polys

def ReadData(infilename):
    """
    read the data from theTAP output file
    """
    # look for the "SITE" line:
    
    infile = file(infilename, 'rU')
    while True:
        line = infile.readline()
        if line.startswith("Site"):
            field_name = line.split('\t')[1].strip()
            if field_name == "% spills exceeding the LOC":
                field_name = "% spills"
            elif field_name == "% of amount released hitting this site":
                field_name = "% of total"
            elif field_name == "Response Time in Hours":
                field_name = "Time (Hrs)"
            else:
                field_name = field_name[:10]
            break
    
    # read the data:
    data = []
    for line in infile:
        d = line.strip().split("\t")[1]
        if d.startswith('<'):
            data.append(0.0)
        elif d.startswith ('>'):
            data.append(float(d[1:]))
        else:
            data.append(float(d))

    return data, field_name

def WritePrjFile(filename, proj='WGS84'):
    if proj <> 'WGS84':
        raise ValueError("WritePrjFile only suporrts the WGS84 coordinates system for now")
    
    # create the PRJ file
    if filename[-4:] == ".shp":
        filename = filename[:-4]
    prj = open("%s.prj" % filename, "w")
    wkt =  'GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]'
    
    #'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
    prj.write(wkt)
    prj.close()

    return None

if __name__ == "__main__":
    if len(sys.argv) == 3:
        try:
            site_file = sys.argv[1]
            output_file = sys.argv[2]
            data_file = None
            data = None
            field_name = None
        except IndexError:
            print USAGE
            sys.exit(1)
    else:
        try:
            site_file = sys.argv[1]
            data_file = sys.argv[2]
            output_file = sys.argv[3]
        except IndexError:
            print USAGE
            sys.exit(1)

    try:
        polys = ReadSites(site_file)
    except:
        print 'There is something wrong with the "SITE.TXT" file'
        raise

    if data_file is not None:
        try:
            data, field_name = ReadData(data_file)
        except:
            print 'There is something wrong with the data file'
            raise
        if len(data) <> len(polys):
            raise ValueError("""The SITE.TXT file doesn't match the TAP export file:
The SITE.TXT file has %i polygons.
The data file hasdata for %i polygons""")

        
    print "Writing shapefile(s)"
    WriteSites(output_file, polys, data, field_name)
    WritePrjFile(output_file)

    