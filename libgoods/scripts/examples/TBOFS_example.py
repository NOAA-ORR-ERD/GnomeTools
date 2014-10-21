#!/usr/bin/env python
from libgoods import romstools, nctools, data_files_dir
import datetime as dt
import os 
import numpy as np

# specify local file or opendap url
#data_url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NGOFS/MODELS/201309/nos.ngofs.fields.f000.20130901.t03z.nc'
sdate = dt.datetime(2014,5,16,0,0)
edate = dt.datetime(2014,5,18,18,0)

flist = []
url_stem = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/TBOFS/MODELS/'
while sdate < edate:
    n = 1
    yr = str(sdate.year)
    mon = str(sdate.month).zfill(2)
    day = str(sdate.day).zfill(2)
    hr = str(sdate.hour).zfill(2)
    url_sub_dir = yr + mon
    while n < 7:
        fname = 'nos.tbofs.fields.n00' + str(n) + '.' + yr + mon + day + '.t' + hr + 'z.nc'
        n = n + 1
        flist.append(url_stem + url_sub_dir + '/' + fname)
    delta_t = dt.timedelta(0.25,0)
    sdate = sdate + delta_t

print "Number of files to aggregate: ", len(flist)
tbofs = romstools.romsgrid(flist)
tbofs.get_dimensions()
nctools.show_tbounds(tbofs.Dataset.variables['ocean_time'])

# subset
nl = 28.2; sl = 27.67
wl = -82.9; el = -82.2
tbofs.subset(nl,sl,wl,el)

print 'getting data'
tbofs.get_data(xindex=tbofs.x,yindex=tbofs.y)
#tbofs.get_data()

print 'writing data'
tbofs.write_it(os.path.join(data_files_dir, 'tbofs_example.nc'))