from libgoods import reg_grid,data_files_dir, nctools
import os
from netCDF4 import date2num
import datetime

out_dir = os.path.join(data_files_dir,'CAROMS')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

var_map = { 'longitude':'lon', \
            'latitude':'lat', \
            'time':'time', \
            'u_velocity':'u', \
            'v_velocity':'v', \
            }

t_units = 'days since 2015-05-01 00:00:00'

#get last two days
sdate = datetime.date.today() - datetime.timedelta(days=2)
edate = datetime.date.today() + datetime.timedelta(days=1)

while sdate < edate:
    day = sdate.day
    month = sdate.month   
    year = sdate.year
    url_stem = 'http://west.rssoffice.com:8080/thredds/dodsC/pacific/CA3km-forecast/CA/ca_subCA_fcst_'

    model_time = str(year) + str(month).zfill(2) + str(day).zfill(2) + '03'
    
    try:
        url = url_stem + model_time + '.nc'
        print url
        caroms = reg_grid.rgrid(url)
        caroms.get_dimensions(var_map)
        
        subset = 1
        tlen = len(caroms.data['time'])
        
#        dt = datetime.datetime(2015,month,day,int(hour),0)
#        caroms.atts['time']['units'] = t_units
#        caroms.data['time'][:] = date2num(dt,t_units)
        
        ofn = os.path.join(out_dir,'scb_' + model_time + '.nc')
        if subset:
            caroms.subset([33,-121,34.8,-117.2])
            caroms.get_data(var_map,yindex=caroms.y,xindex=caroms.x)
            caroms.write_nc(ofn)
        else:
            caroms.get_data(var_map)
            caroms.write_nc(ofn)
    except:
        print 'File not found'
            
    sdate = sdate + datetime.timedelta(days=1)

nctools.make_filelist_for_GNOME(out_dir,'*.nc')