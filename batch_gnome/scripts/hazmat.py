"""
hazmatPy module

This module contains assorted functions and classes that are useful
for NOAA/HAZMAT stuff.

It currently contains:

read_bna(filename,polytype = "list"):

        returns the data in a bna file, given in filename, in twop possible forms:

        If polytype is set to "list" (the default) what is returned is a
        list of two-element tuples. The first element is the polygon
        type, taken from the second field in the header of each polygon
        in the bna. The second element is the polygon. Each polygon is a
        Nx2 NumPy array of floats, such that polygons[n][m,0] is the
        latitude of the mth point of the nth polygon.  polygon[n][m,1]
        is the longitude (in decimal degrees)

        If polytype is set to "PolygonSet", a polygon set, as defined in
        the Geometry module is returned. A polygonset is a set of N
        polygons, such that P[n] returns a NX2 NumPy array of Floats, as
        defined above.


sort_by_other_list(list_to_sort,list_to_sort_by):

        returns a list of the elements of "list_to_sort", sorted by the
        elements of "list_to_sort_by".

        Example:
        >>> hazmat.sort_by_other_list(['a,','b','c','d'],[4,1,3,2])
        ['b', 'd', 'c', 'a,']

"""
# from Numeric import *
from TextFile import open

print "about to define read_bna"

def read_bna(filename,polytype = "list"):
    import string
    file = open(filename,'rt')

    if polytype == "list":
        polygons = []

        line = file.readline()
        while line:
            num_points = int(string.split(line,',')[2])
            polygon_type = string.replace(string.split(line,',')[1],'"','')
            polygon = zeros((num_points,2),Float)
            for i in range(num_points):
                polygon[i,:] = map(float,string.split(file.readline(),','))
            polygons.append((polygon_type,polygon))
            line = string.strip(file.readline())

    elif polytype == "PolygonSet":
        import Geometry

        polygons = Geometry.PolygonSet()

        while 1:
            line = file.readline()
            if not line:
                break
            num_points = int(line.split(',')[2])
            polygon = zeros((num_points,2),Float)
            for i in range(num_points):
                polygon[i,:] = map(float,file.readline().split(','))
            polygons.append(polygon)

    else:
        raise('polytype must be either "list" or "PolygonSet')
    

    file.close()
    return polygons



    #Sorting routine:
def sort_by_other_list(list_to_sort,list_to_sort_by):
    """
    sort_by_other list(list_to_sort,list_to_sort_by)
    
    function that sorts one list by the contents of another list.
    
    the list that is being sorted does not have to be sortable
    """
    pairs = map(None, list_to_sort_by,range(len(list_to_sort_by)))
    pairs.sort()
    out_list = []
    for i in map(lambda x: x[1],pairs):
        out_list.append(list_to_sort[i])
    return out_list
    
    
    
    
    
    
    
    
    
    
    
    
    
    
