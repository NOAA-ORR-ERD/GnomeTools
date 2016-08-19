from libgoods import reg_grid,data_files_dir
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

sdate = datetime.date(2015,6,1)

while sdate < datetime.date(2015,6,10):
    day = sdate.day
    month = sdate.month
    #url_stem = 'http://west.rssoffice.com:8080/thredds/dodsC/pacific/CA3km-nowcast/SCB/ca_subSCB_das_2015' 
    #url_stem = 'http://west.rssoffice.com:8080/thredds/dodsC/pacific/CA3km-forecast'    
    url_stem = 'http://west.rssoffice.com:8080/thredds/dodsC/pacific/CA3km-forecast/CA/'
    times = ['03','09','15','21']

    for hour in times:
        model_time = str(month).zfill(2) + str(day).zfill(2) + hour
        print model_time
        
        try:
            url = url_stem + model_time + '.nc'
            caroms = reg_grid.rgrid(url)
            caroms.get_dimensions(var_map)
            
            subset = 1
            tlen = len(caroms.data['time'])
            
            dt = datetime.datetime(2015,month,day,int(hour),0)
            caroms.atts['time']['units'] = t_units
            caroms.data['time'][:] = date2num(dt,t_units)
            
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
            print url
            
    sdate = sdate + datetime.timedelta(days=1)
