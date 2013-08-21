#!/usr/bin/env python

from netCDF4 import num2date, date2num, date2index
import glob

def show_tbounds(t):
    
    print 'Start date: ', num2date(t[0],t.units)
    try:
        print 'End date: ', num2date(t[-1],t.units)
    except IndexError:
        print num2date(t[0],t.units)
        
def get_tindex(t,start_date,end_date,stride=None):
        
    tindex = []
    tindex.append(date2index(start_date,t,select='before'))
    tindex.append(date2index(end_date,t,select='after') + 1)
    if stride is None:
        tindex.append(1)
    else:
        tindex.append(stride)
    return tindex

def adjust_time(t,t_units):
    #GNOME can't handle pre-1970 date units
    
    dtime = num2date(t,t_units)
    new_units = 'days since 1980-1-1 00:00:00'
    new_time = date2num(dtime,units=new_units)
    
    return new_time,new_units

def make_filelist_for_GNOME(file_dir,file_match='*.*'):
    pass
    