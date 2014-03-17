#!/usr/bin/env python


usage = """
nc_time_shift(infilename, outfilename, amount [units=days])

Script to time-shift a CF-compliant netcdf file.

It does this by changing the base date of the time units.

units can be one of weeks, days, hours, minutes, or seconds 

If infilename and outfilename are the same, the file will be changed in-place.

"""

import sys
import shutil
import netCDF4

def shift_time(infilename, outfilename, amount, shift_units):
    """
    Do the actual time shift

    :param infilename: the name of the input file (full or relative path)
    :param outfilename: the name of the output file (full or relative path)

    :param amount: amount to shift the time
    :type amount: floating point number, of a string version

    :param shift_units: units of the time shift
    :type shift_units: string
    """

    units = str(shift_units).strip()

    # make a copy of the file if it's different
    if not infilename == outfilename:
        print "creating new file:", outfilename
        shutil.copyfile(infilename, outfilename)

    # pull the time units out of the source copy:
    nc = netCDF4.Dataset(infilename, mode='r')
    time_var = nc.variables['time']
    time_units = time_var.units
    print "original time units", time_units
    nc.close()

    units, epoch = time_units.split('since')
    # iso datetimes can be separated with whitesapce or "T"
    try:
        date, time = epoch.strip().split()
    except ValueError: 
        date, time = epoch.strip().split("T")
    year, month, day = [int(i) for i in date.split('-')]
    hour, minute, second = [int(i) for i in time.split(':')]
    
    dt = netCDF4.datetime(year, month, day, hour, minute, second)

    try:
        shift = netCDF4.timedelta(**{shift_units:float(amount)})
    except TypeError: 
        raise ValueError("Time shift must units options are: days, seconds, minutes, hours, or weeks")

    dt+=shift
    new_units = "%s since %s"%(units.strip(), dt.isoformat())
    print "new time units:", new_units

    ## write out to the new file:
    nc = netCDF4.Dataset(outfilename, mode='r+')
    time_var = nc.variables['time']
    time_var.units = new_units
    nc.close()

    print "Done"


if "__name__ == __main__":
    try:
        infilename, outfilename, amount = sys.argv[1:4]
    except ValueError:
         print "*ERROR*\nYou must pass in infilename, outfilename, and amount"
         print usage
         sys.exit(1)
    try: 
        units = sys.argv[4]
    except IndexError:
        units = 'days'

    print infilename, outfilename, amount, units 

    shift_time(infilename, outfilename, amount, units)


