#!/usr/bin/env python
import datetime as dt
import urllib, os
import numpy as np
from netCDF4 import Dataset, date2num

'''
!!!This script was written to be a stand alone without needing to 
install the libgoods package -- however, here I modify it to use
the libgoods data_files_dir 
'''
from libgoods import data_files_dir

'''
Sample script to download and reformat model output from the SLGO web service
Edit the main loop for the desired start and end dates or to point to local
files if output has already been downloaded (in .INI text format)

Scipt contains three methods that shouldn't need to be edited:
grid_latlon: returns lat/lon vectors based on grid dimensions
read_text_format: reads the INI format for time, lon, u and v plus descriptors (attributes)
write_reg_grid: writes the model data into a GNOME compatible netCDF file

'''
def grid_latlon(model):
    
    if model.lower() == 'stle':
        min_lon = -72.6
        max_lon = -68.0
        min_lat = 46.0
        max_lat = 49.5
        num_east_west = 1150
        num_north_south = 1400
    elif model.lower() == 'mogsl':
        min_lon = -71
        max_lon = -56.030000334605575
        min_lat = 45.5
        max_lat = 51.97999985516071
        num_east_west = 500
        num_north_south = 325
    else:
        return [],[]
    lon = np.linspace(float(min_lon),float(max_lon),int(num_east_west))
    lat = np.linspace(float(min_lat),float(max_lat),int(num_north_south))
    
    return lat,lon

def read_text_format(data,lat,lon):
    
    # Read and parse data from the "INI" (text) format downloaded from the SLGO web service
    # data is a file object (can be loaded from local file or url using urllib.urlopen
    atts = dict()
    keep_reading = True
    while keep_reading:
        line = data.readline()
        if line.startswith('Erreur'):
            break
        elif line.startswith('[Mask'):
            read_mask = 1
            mask = np.zeros((len(lat),len(lon)))
            skip = data.readline()
            skip = data.readline()
            while read_mask:
                row = data.readline()
                try:
                    row_num = int(row.split('=')[0].split('_')[1]) - 1
                    values = [int(x) for x in row.split('=')[1][1:-1]]
                    mask[row_num,:] = values
                except:
                    read_mask = 0
        elif line.startswith('[Current_def'):
            atts['u'] = {'units' : data.readline().split(':')[1].strip()}
            atts['v'] = {'units' : data.readline().split(':')[1].strip()}
        elif line.startswith('[Time_Def'):
            t_units = 'days since 2012-01-01 00:00:00'
            atts['t'] = {'units' : t_units}
            skip = data.readline()
            tlen = int(data.readline().split('=')[1])
            ts = []
            for ii in range(tlen):
                d = data.readline().split('= ')[1].split()
                d = [int(x) for x in d]
                ts.append(dt.datetime(d[0],d[1],d[2],d[3],d[4],d[5]))
            time = date2num(ts,t_units)
        elif line.startswith('[Current_Forecast'):
            read_currents = 1
            skip = data.readline()
            u = np.zeros((tlen,len(lat),len(lon)))
            v = np.zeros((tlen,len(lat),len(lon)))
            while read_currents:
                line = data.readline()
                try:
                    idxs,cdata = line.split('=')
                    col,row = idxs.split('_')
                    line_f = [float(x) for x in cdata.split()]
                    u[:,row,col] = line_f[::2]
                    v[:,row,col] = line_f[1::2]
                except:
                    read_currents = 0
            keep_reading = False
                 
    atts['mask'] = mask
    
    if atts['u']['units'] == 'mm/s': #not supported in GNOME
        u = u/1000
        v = v/1000
        atts['u']['units'] = 'm/s'
        atts['v']['units'] = 'm/s'
        
    return time, u, v, atts

def write_reg_grid(time,lon,lat,u,v,atts,ofn):
    """
  
    Write GNOME compatible netCDF file (netCDF3) from regular grid data
    
    
    subset(time,lon,lat,u,v,atts,dfn)
   
    Arguments: 
       *time*: 1-D time vector
       *lon* : 1-D in decimal degrees -180:180
       *lat* : 1-D in decimal degrees -90:90
       *u*   : 2-D eastward velocity component
       *v*   : 2-D northward velocity component
       *atts*: nested dict object of variable attributes with keys u, v, and t
       *ofn* : netCDF output file name -- string

    """
    
    nc = Dataset(ofn,'w',format='NETCDF3_CLASSIC')
    
    # Global Attributes
    setattr(nc,'grid_type','regular')

    # add Dimensions
    nc.createDimension('lon',len(lon))
    nc.createDimension('lat',len(lat))
    nc.createDimension('time',None)

    try:
        ufill = atts['u']['_FillValue']
        vfill = atts['v']['_FillValue']
    except KeyError:
        try:
            ufill = atts['u']['missing_value']
            vfill = atts['v']['missing_value']
        except KeyError:
            ufill = 999.
            vfill = 999.

    # create variables
    nc_time = nc.createVariable('time','f8',('time',))
    nc_lon = nc.createVariable('lon','f4',('lon',))
    nc_lat = nc.createVariable('lat','f4',('lat',))
    if atts.has_key('wind'):
        nc_u = nc.createVariable('air_u','f4',('time','lat','lon'), \
            fill_value=ufill)
        nc_v = nc.createVariable('air_v','f4',('time','lat','lon'), \
            fill_value=vfill)
    else:
        nc_u = nc.createVariable('water_u','f4',('time','lat','lon'), \
            fill_value=ufill)
        nc_v = nc.createVariable('water_v','f4',('time','lat','lon'), \
            fill_value=vfill)
    
    # add data
    nc_lon[:] = lon
    nc_lat[:] = lat
    if len(u.shape) == 2:
        nc_time[0] = time
        nc_u[0,:] = u
        nc_v[0,:] = v
    else:
        nc_time[:] = time
        nc_u[:] = u
        nc_v[:] = v

    # add variable attributes from 'atts' (nested dict object)
    for an_att in atts['t'].iteritems():
        setattr(nc_time,an_att[0],an_att[1])
        
    for an_att in atts['u'].iteritems():
        if an_att[0] != '_FillValue':
            setattr(nc_u,an_att[0],an_att[1])

    for an_att in atts['v'].iteritems():
        if an_att[0] != '_FillValue':
            setattr(nc_v,an_att[0],an_att[1])
    
            
    nc.close()
    
def main():

    # change model and/or bounding times (start=sdate, end=edate)
    model = 'mogsl'
    sdate = dt.date(2013,9,3)
    edate = dt.date(2013,9,6)
    
    # i just hard-coded the grid dimensions for the stle and mogsl models rather
    # than read the info from the file each time
    lat,lon = grid_latlon(model)
    
    # download data from web service in 1 day chunks and write to GNOME netCDF files
    # if data has already been downloaded, edit this loop to generate sequential filenames
    # then use data = open(filename) instead of data = urllib.urlopen(q_url)
    ofns = [] #this will hold the names of the output GNOME netCDF files -- used to generate a text file with the list of filenames
    while sdate < edate:
        sdate_plus1 = sdate + dt.timedelta(days=1)
        date_query_str1 = str(sdate.year) + str(sdate.month).zfill(2) + str(sdate.day).zfill(2) + '000000'
        date_query_str2 = str(sdate_plus1.year) + str(sdate_plus1.month).zfill(2) + str(sdate_plus1.day).zfill(2) + '000000'
        q_url = 'http://ws.ns-shc.qc.dfo-mpo.gc.ca/OO-CurrentsIceWeb/ExportData?model=' + model + \
            '&format=text&data=u,v&datemin=' + date_query_str1 + '&datemax=' + date_query_str2
        print q_url
        try:
            data = urllib.urlopen(q_url)
            # data = open('C:\\Users\\amy.macfadyen\\Downloads\\mogsldata_' + sdate.isoformat() + '.ini')
            time,u,v,atts = read_text_format(data,lat,lon)
            ofn = model + '_' + sdate.isoformat() + 'to' + sdate_plus1.isoformat() + '.nc'
            write_reg_grid(time,lon,lat,u,v,atts,os.path.join(data_files_dir,ofn))
            ofns.append(ofn)
        except:
            print 'Data download must have failed'
        finally:
            sdate = sdate_plus1
    
    # write a filelist of the output netCDF files -- can load this into GNOME as a pointer to all the files
    f = open(os.path.join(data_files_dir,model + '_filelist.txt'),'w')
    f.write('NetCDF Files\n')
    f.write('\n'.join(['[FILE] ' + file for file in ofns]))
    f.close()

    
if __name__ == "__main__":  main()



    
    
        

