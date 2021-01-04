#!/usr/bin/env python

from __future__ import print_function
from netCDF4 import Dataset, num2date, date2num, date2index
import glob
import os
import math
import datetime as dt


def get_var_map(filename, var_list=['longitude',
                                    'latitude',
                                    'time',
                                    'u_velocity',
                                    'v_velocity',
                                    'air_u',
                                    'air_v']):

    var_map = dict()
    ncvars = Dataset(filename).variables
    long_names = {'longitude': 'longitude',
                  'latitude': 'latitude',
                  'time': 'time',
                  'u_velocity': 'eastward_sea_water_velocity',
                  'v_velocity': 'northward_sea_water_velocity',
                  'air_u': 'eastward_wind',
                  'air_v': 'northward_wind',
                  }

    for var in var_list:
        matches = []
        for varname in ncvars:
            try:
                if ncvars[varname].standard_name == long_names[var]:
                    matches.append(varname)
            except AttributeError:
                pass
        if not matches:
            var_map[var] = None
        elif len(matches) == 1:
            var_map[var] = matches[0]
        else:
            var_map[var] = matches

    return var_map

def sync_filelist_time_units(file_list,time_var='time'):
    '''
    Make all times relative to units in first file
    '''
    first_file = True
    for f in file_list:
        nc = Dataset(f,'r+')
        t = nc.variables[time_var]
        if first_file:
            first_file = False
            base_units = t.units
        else:
            ts = num2date(t[:],t.units)
            t[:] = date2num(ts,base_units)
            t.units = base_units
        print('Start date: ', num2date(t[0], t.units))
        nc.close()
    

def when(filename,tvar='time'):
        nc = Dataset(filename)
        t = nc.variables[tvar]
        print('Start date: ', num2date(t[0], t.units))
        try:
            print('End date: ', num2date(t[-1], t.units))
        except IndexError:
            print('End date: ', num2date(t[0], t.units))
        nc.close()
            
def shift_time(filename, tshift, tvar='time'):
    '''
    Shift time by hours in tshift
    for reference:
          [[0,'GMT'],
          [-10,'HST (-10)'],
          [-9,'AKST (-9)'],
          [-8,'AKDT (-8)'],
          [-8,'PST (-8)'],
          [-7,'PDT (-7)'],
          [-7,'MST (-7)'],
          [-6,'MDT (-6)'],       
          [-6,'CST (-6)'],
          [-5,'CDT (-5)'],
          [-5,'EST (-5)'],
          [-4,'EDT (-4)'],
          [-3,'ADT (-3)']]
    '''
    nc = Dataset(filename,'r+')
    t = nc.variables[tvar]
    
    offset = dt.timedelta(hours = tshift)
    oldtime = num2date(t[:],t.units)
    newtime = date2num(oldtime + offset, t.units)
    t[:] = newtime
    nc.close()
    
def fix_time_units(units):
    '''
    GNOME doesn't support units of this form: 'hours since 2014-12-12T18:00:00Z'
    Just replace T with space
    '''
    new_units = ' '.join(units.split('T'))
    return new_units


def show_ncfile_tbounds(filename, tvar='time'):
    '''
    pass in a netcdf filename to see time bounds
    assumes time variable is "time" unless specified.
    '''
    t = Dataset(filename).variables[tvar]
    print('Start date: ', num2date(t[0], t.units))
    try:
        print('End date: ', num2date(t[-1], t.units))
    except IndexError:
        print(num2date(t[0], t.units))


def show_tbounds(t):

    print('Start date: ', num2date(t[0], t.units))
    try:
        print('End date: ', num2date(t[-1], t.units))
    except IndexError:
        print(num2date(t[0], t.units))


def get_tindex(t, start_date, end_date, stride=None):

    tindex = []
    tindex.append(date2index(start_date, t, select='before'))
    tindex.append(date2index(end_date, t, select='after') + 1)
    if stride is None:
        tindex.append(1)
    else:
        tindex.append(stride)
    return tindex


def adjust_time(t, t_units):
    # GNOME can't handle pre-1970 date units
    # fixme: but using days ???

    dtime = num2date(t, t_units)
    new_units = 'days since 1980-1-1 00:00:00'
    new_time = date2num(dtime, units=new_units)

    return new_time, new_units


def round_time(datetime_in=None, roundto=60):
    """Round a datetime object to any time laps in seconds
    dt : datetime.datetime object
    roundTo : Closest number of seconds to round to, default 1 minute.
    Author: Thierry Husson 2012 - Use it as you want but don't blame me.
    """
    seconds = (datetime_in - datetime_in.min).seconds
    # // is a floor division, not a comment on following line:
    rounding = (seconds + roundto / 2) // roundto * roundto
    return datetime_in + dt.timedelta(0, rounding - seconds, -datetime_in.microsecond)


def make_filelist_for_GNOME(file_dir, file_match='*.*', outfilename='filelist.txt'):
    # used to load multiple files as one mover in GNOME
    flist = glob.glob(os.path.join(file_dir, file_match))
    f = open(os.path.join(file_dir, outfilename), 'w')
    f.write('NetCDF Files\n')
    # f.write('\n'.join(['[FILE] ' + os.path.split(file)[-1] for file in flist]))
    f.write('\n'.join(['[FILE] ' + file for file in flist]))
    f.close()


def get_time_var(ds):
    """
    get the time variable from a dataset

    looks for standard name

    :param ds: open Dataset object
    """
    var_list = ds.get_variables_by_attributes(standard_name='time')
    if len(var_list) == 0:
        # not found by standard name -- try plain old time name
        try:
            var = ds.variables['time']
        except KeyError:
            raise ValueError('no variable with name or standard_name: "time"')
    elif len(var_list) == 1:
        var = var_list[0]
    elif len(var_list) >= 1:
        raise ValueError("more than one time variable in dataset")
    return var


def utmToLatLng(zone, easting, northing, northernHemisphere=True):
    # Convert UTM coordinates to lat/lon
    # fixme: WOW! maybe use proj4 next time?

    if not northernHemisphere:
        northing = 10000000 - northing

    a = 6378137
    e = 0.081819191
    e1sq = 0.006739497
    k0 = 0.9996

    arc = northing / k0
    mu = arc / (a * (1 - math.pow(e, 2) / 4.0 - 3 *
                     math.pow(e, 4) / 64.0 - 5 * math.pow(e, 6) / 256.0))

    ei = (1 - math.pow((1 - e * e), (1 / 2.0))) / \
        (1 + math.pow((1 - e * e), (1 / 2.0)))

    ca = 3 * ei / 2 - 27 * math.pow(ei, 3) / 32.0

    cb = 21 * math.pow(ei, 2) / 16 - 55 * math.pow(ei, 4) / 32
    cc = 151 * math.pow(ei, 3) / 96
    cd = 1097 * math.pow(ei, 4) / 512
    phi1 = mu + ca * math.sin(2 * mu) + cb * math.sin(4 * mu) + \
        cc * math.sin(6 * mu) + cd * math.sin(8 * mu)

    n0 = a / math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (1 / 2.0))

    r0 = a * (1 - e * e) / \
        math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (3 / 2.0))
    fact1 = n0 * math.tan(phi1) / r0

    _a1 = 500000 - easting
    dd0 = _a1 / (n0 * k0)
    fact2 = dd0 * dd0 / 2

    t0 = math.pow(math.tan(phi1), 2)
    Q0 = e1sq * math.pow(math.cos(phi1), 2)
    fact3 = (5 + 3 * t0 + 10 * Q0 - 4 * Q0 * Q0 -
             9 * e1sq) * math.pow(dd0, 4) / 24

    fact4 = (61 + 90 * t0 + 298 * Q0 + 45 * t0 * t0 - 252 *
             e1sq - 3 * Q0 * Q0) * math.pow(dd0, 6) / 720

    lof1 = _a1 / (n0 * k0)
    lof2 = (1 + 2 * t0 + Q0) * math.pow(dd0, 3) / 6.0
    lof3 = (5 - 2 * Q0 + 28 * t0 - 3 * math.pow(Q0, 2) + 8 *
            e1sq + 24 * math.pow(t0, 2)) * math.pow(dd0, 5) / 120
    _a2 = (lof1 - lof2 + lof3) / math.cos(phi1)
    _a3 = _a2 * 180 / math.pi

    latitude = 180 * (phi1 - fact1 * (fact2 + fact3 + fact4)) / math.pi

    if not northernHemisphere:
        latitude = -latitude

    longitude = ((zone > 0) and (6 * zone - 183.0) or 3.0) - _a3

    return latitude, longitude
