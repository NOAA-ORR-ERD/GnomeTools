#!/usr/bin/env python

from netCDF4 import Dataset
from post_gnome import nc_particles
import shapefile as shp
import os
import zipfile
import numpy as np
import scipy.stats as st
from matplotlib.pyplot import contour


def write_proj_file(prj_filename):
    prj = open("%s.prj" % prj_filename, "w")
    epsg = 'GEOGCS["WGS 84",'
    epsg += 'DATUM["WGS_1984",'
    epsg += 'SPHEROID["WGS 84",6378137,298.257223563]]'
    epsg += ',PRIMEM["Greenwich",0],'
    epsg += 'UNIT["degree",0.0174532925199433]]'
    prj.write(epsg)
    prj.close()


def points(fn, package_dir, t2convert):

    nc = Dataset(fn)
    particles = nc_particles.Reader(nc)
    times = particles.times
    dt = [np.abs(((output_t - t2convert).total_seconds()) / 3600) for output_t in times]
    t = dt.index(min(dt))
    print 'Converting output from: ', times[t]

    w = shp.Writer(shp.POINT)
    w.autobalance = 1

    w.field('Year', 'C')
    w.field('Month', 'C')
    w.field('Day', 'C')
    w.field('Hour', 'C')
    w.field('LE id', 'N')
    w.field('Depth', 'N')
    w.field('Mass', 'N')
    w.field('Age', 'N')
    w.field('Status_Code', 'N')

    TheData = particles.get_timestep(t, variables=['latitude',
                                                   'longitude',
                                                   'id',
                                                   'depth',
                                                   'mass',
                                                   'age',
                                                   'status_codes'])
    for k, p in enumerate(zip(TheData['longitude'], TheData['latitude'])):
        w.point(p[0],p[1])
        w.record(times[t].year,
                 times[t].month,
                 times[t].day,
                 times[t].hour,
                 TheData['id'][k],
                 TheData['depth'][k],
                 TheData['mass'][k],
                 TheData['age'][k],
                 TheData['status_codes'][k])

    source_fdir = os.path.join(package_dir, 'source_files')
    shapefile_name = os.path.split(fn)[-1].split('.')[0]
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

    w.field('Year', 'C')
    w.field('Month', 'C')
    w.field('Day', 'C')
    w.field('Hour', 'C')
    w.field('Depth', 'N')
    w.field('Type', 'C')

    for c in range(len(cs.collections)):
        p = cs.collections[c].get_paths()[0]
        v = p.vertices
        coords = [[[i[0], i[1]] for i in v]]
        w.poly(shapeType=3, parts=coords)
        w.record(times[t].year, times[t].month, times[t].day, times[t].hour,
                 TheData['depth'][c], names[c])
        print names[c]

    source_fdir = os.path.join(package_dir, 'source_files')
    shapefile_name = os.path.split(fn)[-1].split('.')[0]
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
