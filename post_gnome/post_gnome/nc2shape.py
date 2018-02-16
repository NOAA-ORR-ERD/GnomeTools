#!/usr/bin/env python

from netCDF4 import Dataset
from post_gnome import nc_particles
import shapefile as shp
import os,glob
import zipfile
import numpy as np
import scipy.stats as st
from matplotlib.pyplot import contour
from shapely import geometry
import datetime, pytz

def write_proj_file(prj_filename):
    prj = open("%s.prj" % prj_filename, "w")
    epsg = 'GEOGCS["WGS 84",'
    epsg += 'DATUM["WGS_1984",'
    epsg += 'SPHEROID["WGS 84",6378137,298.257223563]]'
    epsg += ',PRIMEM["Greenwich",0],'
    epsg += 'UNIT["degree",0.0174532925199433]]'
    prj.write(epsg)
    prj.close()
    
def zip_shape_files(shp_fdir,shp_fname):
    print 'zipping'
    files = glob.glob(os.path.join(shp_fdir, shp_fname) + '*')
    print files
    zipfname = os.path.join(shp_fdir,shp_fname) + '.zip'
    zipf = zipfile.ZipFile(zipfname, 'w',compression=zipfile.ZIP_DEFLATED)
    for f in files:
        if f.split('.')[0] == os.path.join(shp_fdir, shp_fname) and f.split('.')[1] != 'zip':
            print f
            zipf.write(f, arcname=os.path.split(f)[-1])
            os.remove(f)
    zipf.close()
    return zipfname
    
def mpl_contours_2_shape(cs,levs,t,shp_fdir,shp_fname):
   
    w = shp.Writer(shp.POLYGON)
    w.autobalance = 1
    w.field('Time', 'C')
    w.field('Value', 'C')
    
    for c,col in enumerate(cs.collections):
        # Loop through all polygons that have the same intensity level
        for contour_path in col.get_paths(): 
            # Create the polygon for this intensity level
            # The first polygon in the path is the main one, the following ones are "holes"
            # Right now I am only using the exterior so maybe no holes?
            for ncp,cp in enumerate(contour_path.to_polygons()):
                x = cp[:,0]
                y = cp[:,1]
                
                new_shape = geometry.Polygon([(i[0], i[1]) for i in zip(x,y)])
                if ncp == 0:
                    poly = new_shape
                else:
                    # Remove the holes if there are any
                    poly = poly.difference(new_shape)

            # do something with polygon
            ex_x, ex_y = poly.exterior.xy
            coords = [[[ex_x[i], ex_y[i]] for i in range(len(ex_x))]]
            w.poly(shapeType=5, parts=coords)
            w.record(t.isoformat(), '> ' + str(levs[c]))

    w.save(os.path.join(shp_fdir, shp_fname))
        
    # create the PRJ file
    prj_filename = os.path.join(shp_fdir, shp_fname)
    write_proj_file(prj_filename)

    zip_shape_files(shp_fdir,shp_fname)

def points(fn, package_dir, t2convert):
    nc = Dataset(fn)
    particles = nc_particles.Reader(nc)
    times = particles.times

    timezone = pytz.timezone("America/Los_Angeles") #this should be passed in
    
    w = shp.Writer(shp.POINT)
    w.autobalance = 1
    w.field('Time', 'C')
    w.field('LE id', 'N')
    w.field('Depth', 'N')
    w.field('Mass', 'N')
    w.field('Age', 'N')
    w.field('status', 'N')
    
    if t2convert is None:
        ts = range(0,len(times))
    else:
        dt = [np.abs(((output_t - t2convert).total_seconds()) / 3600) for output_t in times]
        ts = [dt.index(min(dt))]
    
    for t in ts:
        print 'Converting output from: ', times[t]
        TheData = particles.get_timestep(t, variables=['latitude',
                                                       'longitude',
                                                       'id',
                                                       'depth',
                                                       'mass',
                                                       'age',
                                                       'status_codes'])
        #new_t = datetime.datetime(times[t].year,times[t].month,times[t].day,times[t].hour,times[t].minute,times[t].second,tzinfo=timezone)
        for k, p in enumerate(zip(TheData['longitude'], TheData['latitude'])):
            if TheData['status_codes'][k] != 10:
                w.point(p[0],p[1])
                w.record(times[t].strftime('%Y-%m-%dT%H:%M:%S -7:00'),
                         TheData['id'][k],
                         TheData['depth'][k],
                         TheData['mass'][k],
                         TheData['age'][k],
                         TheData['status_codes'][k])

    source_fdir = os.path.join(package_dir, 'source_files')
    #shapefile_name = os.path.split(fn)[-1].split('.')[0]
    if t2convert is None:
        shapefile_name = os.path.split(fn)[-1].split('.')[0] + '_all'
    else:
        shapefile_name = os.path.split(fn)[-1].split('.')[0] + '_' + times[t].strftime('%Y%b%d_%H%M')
    print 'sfn:', shapefile_name
    print os.path.join(source_fdir, shapefile_name)

    w.save(os.path.join(source_fdir, shapefile_name))

    nc.close()

    # create the PRJ file
    prj_filename = os.path.join(source_fdir, shapefile_name)
    write_proj_file(prj_filename)

    zipfname = zip_shape_files(source_fdir,shapefile_name)
    print 'zipfname', zipfname

    return zipfname


def contours(fn,
             package_dir,
             t2convert,
             levels=[0.1, 0.4, 0.8],
             names=['Light', 'Medium', 'Heavy'],
             ):

    print "contouring data in:", fn
    nc = Dataset(fn)
    particles = nc_particles.Reader(nc)
    times = particles.times
    dt = [np.abs(((output_t - t2convert).total_seconds()) / 3600) for output_t in times]
    t = dt.index(min(dt))
    print 'Converting output from: ', times[t]

    TheData = particles.get_timestep(t, variables=['latitude',
                                                   'longitude',
                                                   'id',
                                                   'depth',
                                                   'mass',
                                                   'age',
                                                   'status_codes'])

    # contouring
    status = TheData['status_codes']
    floating = np.where(status == 2)[0]
    x = TheData['longitude'][floating]
    y = TheData['latitude'][floating]
    
    # Peform the kernel density estimate
    xx, yy = np.mgrid[min(x) - .1:max(x) + .1:100j, min(y) - .1:max(y) + .1:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    max_density = f.max()

    levels.sort()
    particle_contours = [lev * max_density for lev in levels]

    cs = contour(xx, yy, f, particle_contours)

    w = shp.Writer(shp.POLYGON)
    w.autobalance = 1

    w.field('Time', 'C')
    w.field('Depth', 'N')
    w.field('Type', 'C')

    for c in range(len(cs.collections)):
        p = cs.collections[c].get_paths()[0]
        v = p.vertices
        coords = [[[i[0], i[1]] for i in v]]
        w.poly(shapeType=5, parts=coords)
        w.record(times[t].isoformat(),TheData['depth'][c], names[c])
        print names[c]

    source_fdir = os.path.join(package_dir, 'source_files')
    shapefile_name = os.path.split(fn)[-1].split('.')[0] + '_contours'
    w.save(os.path.join(source_fdir, shapefile_name))

    nc.close()

    # create the PRJ file
    prj_filename = os.path.join(source_fdir, shapefile_name)
    write_proj_file(prj_filename)

    files = os.listdir(source_fdir)
    zipfname = shapefile_name + '.zip'
    zipf = zipfile.ZipFile(os.path.join(source_fdir, zipfname), 'w')
    for f in files:
        if f.split('.')[0] == shapefile_name:
            zipf.write(os.path.join(source_fdir, f), arcname=f)
    zipf.close()

    return zipfname
