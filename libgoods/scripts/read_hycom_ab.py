# -*- coding: utf-8 -*-

def ReadHYCOMMercatorParameters(gridAFile, gridBFile, gridCode):

    # Read the number of rows and columns from regional.grid.b.

    f = open(gridBFile, 'rU')
    rows = int(f.readline().split()[0])
    cols = int(f.readline().split()[0])
    f.close()

    Logger.Info('This HYCOM grid has %i rows and %i columns.' % (rows, cols))

    # Compute the number of bytes of a variable stored in a HYCOM .a
    # file for a grid of this size.
    
    varByteLength = int(math.ceil(float(cols*rows*4) / 16384) * 16384)

    Logger.Info('Each 2D variable in a HYCOM .a file takes %i bytes, including the padding to a 16 KB boundary.' % varByteLength)

    # Read the longitude and latitude grids as numpy arrays.

    import numpy

    Logger.Info('Reading the longitude and latitude grids from %s.' % gridAFile)

    if gridCode == 'p':
        varsToSkip = 0
    elif gridCode == 'q':
        varsToSkip = 2
    elif gridCode == 'u':
        varsToSkip = 4
    elif gridCode == 'v':
        varsToSkip = 6

    f = open(gridAFile, 'rb')
    f.seek(varsToSkip * varByteLength)
    lon = numpy.ascontiguousarray(numpy.flipud(numpy.cast['float'](numpy.fromfile(f, '>f4', rows*cols)).reshape((cols,rows))))
    f.seek((varsToSkip+1) * varByteLength)
    lat = numpy.ascontiguousarray(numpy.flipud(numpy.cast['float'](numpy.fromfile(f, '>f4', rows*cols)).reshape((cols,rows))))
    f.close()

    # I was told by Alan Wallcraft that the grid coordinates are only
    # accurate to four decimal digits (Fortran REAL*4 precision).
    # Round the coordinates to four digits.

    lon = lon.round(4)
    lat = lat.round(4)
    
    return lon,lat
    

if "__name__ == __main__":
    
    
    
    # Calculate the cell size. Alan said that the Mercator portion of
    # HYCOM's hybrid grid is based on a sphere with radius 6371001.0
    # meters.

#    p = Proj(proj='merc', ellps='sphere', a=6371001.0, b=6371001.0)
#    equatorialCellWidthInDegrees = (lon.max() - lon.min())/lon.shape[1]
#    cellSize = p(equatorialCellWidthInDegrees, 0)[0]
#
#    Logger.Info('The cell size of the Mercator section of the grid is %s m.' % repr(cellSize))
#
#    # The HYCOM grid is not centered on 0 longitude. The left-most
#    # edge is typically around 74 E. Determine the central meridian
#    # and build the ArcGIS WKT string for the projection. Keep in mind
#    # that the HYCOM grid coordinates are for the centers of the
#    # cells, not the corners.
#    #
#    # This calculation assumes the grid is global.
#
#    indexOfZeroLat = lat[:,0].tolist().index(0.0)
#    centralMeridian = -180.0 + lon[indexOfZeroLat,0] - equatorialCellWidthInDegrees/2
#    coordSys = "PROJCS['HYCOM_Mercator',GEOGCS['HYCOM_Sphere',DATUM['D_Sphere',SPHEROID['Sphere',6371001.0,0.0]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Mercator'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian'," + repr(centralMeridian) + "],PARAMETER['Standard_Parallel_1',0.0],UNIT['Meter',1.0]]"
#
#    Logger.Info('The central meridian is %s.' % repr(centralMeridian))
#
#    # Find the y index of the bottom of the Mercator section of the
#    # grid. Below this section is a section in which the latitude
#    # increases by a constant value, typically 0.032 degrees. Find
#    # where this section starts by looking for the y index at which
#    # latitude increases by a constant value.
#
#    indexOfBottom = (abs(lat[indexOfZeroLat:-1,0]-lat[indexOfZeroLat+1:,0] - (lat[-3,0]-lat[-2,0])) < 0.000001).tolist().index(True) + indexOfZeroLat
#
#    # Calculate the Mercator coordinates of the lower left corner. We
#    # will extract the entire x extent; assuming that extent is
#    # global, the x coordinate is just 1/2 way around the globe. The y
#    # coordinate is the projected latitude of the southmost cell of
#    # the Mercator section of the grid.
#
#    xLLC = p(-180.0, 0)[0]
#    yLLC = p(0, lat[indexOfBottom,0])[1] - cellSize/2
#
#    # Above the Mercator section of the grid is a section that has a
#    # bipoler projection. This typically starts around 47 N. Find the
#    # index of the top of the Mercator section by computing the
#    # difference in projected cell height from what is expected in the
#    # Mercator projection. In the bipolar section, the difference is
#    # very large (1000+ meters). Look for the cell at which the
#    # difference drops below 10 meters. It is a dramatic change, so
#    # this should be a reliable technique.
#
#    deltaY = numpy.array(map(lambda a, b: abs(p(0,a)[1] - p(0,b)[1] - cellSize), lat[:-1,0], lat[1:,0]))
#    indexOfTop = (deltaY < 10).tolist().index(True)
#
#    Logger.Info('The Mercator section of the grid spans rows %i through %i, inclusive, where the top row is 0.' % (indexOfTop, indexOfBottom))
#    Logger.Info('The center latitudes of the top and bottom rows of the Mercator section are %.4f and %.4f.' % (lat[indexOfTop,0], lat[indexOfBottom,0]))
#    Logger.Info('The projected x and y coordinates of the lower-left corner of the Mercator section are %s, %s.' % (repr(xLLC), repr(yLLC)))
#
#    # Return successfully.
#
#    return rows, cols, indexOfTop, indexOfBottom, xLLC, yLLC, cellSize, coordSys, lon, lat