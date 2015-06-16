from libgoods import reg_grid
import os

out_dir = 'C:\\Users\\amy.macfadyen\\Desktop\\CAROMS'


var_map = { 'longitude':'lon', \
            'latitude':'lat', \
            'time':'time', \
            'u_velocity':'u', \
            'v_velocity':'v', \
            }

url_stem = 'http://west.rssoffice.com:8080/thredds/dodsC/pacific/CA3km-nowcast/SCB/ca_subSCB_das_201505'
#url_stem = 'http://west.rssoffice.com:8080/thredds/dodsC/pacific/CA3km-nowcast/CA/ca_subCA_das_201505'
times = ['03','09','15','21']
for day in range(19,31):
    for hour in times:
        model_time = str(day).zfill(2) + hour
        print model_time
        
        url = url_stem + model_time + '.nc'
        caroms = reg_grid.rgrid(url)
        caroms.get_dimensions(var_map)
        
        subset = 1
        tlen = len(caroms.data['time'])
        
        caroms.atts['time']['units'] = 'days since 2015-05-' + str(day) + '-'
        caroms.data['time'][:] = 0.0
        
        ofn = os.path.join(out_dir,'scb_' + model_time + '.nc')
        if subset:
            caroms.subset([33,-121,34.8,-117.2])
            caroms.get_data(var_map,yindex=caroms.y,xindex=caroms.x)
            caroms.write_nc(ofn)
        else:
            caroms.get_data(var_map)
            caroms.write_nc(ofn)
