#!/usr/bin/env python
import datetime
import urllib2, requests
from netCDF4 import Dataset, num2date
import os, glob
'''
Methods for generating ordered filelist for a time series of CO-OPS data
(Nowcasts + Forecasts) based on user specified start and end dates. If end 
date is unspecified or greater than datetime.utcnow() the latest forecast
will be automatically be appended.

Notes on COOPS naming and aggregations:
Nowcast and forecast files are created four times a day. Output is hourly in 
individual files. So each update generates 6 nowcast files and 48 forecast files
The update cycle time will be the last model output timestep in the nowcast 
files and the first timestep in the forecast files

Example filenames from one update cycle (20141027.t15z):
Nowcast:
nos.ngofs.fields.n000.20141027.t15z.nc
nos.ngofs.fields.n001.20141027.t15z.nc 
...
nos.ngofs.fields.n006.20141027.t15z.nc 

Forecast:
nos.ngofs.fields.f000.20141027.t15z.nc
nos.ngofs.fields.f002.20141027.t15z.nc
...
nos.ngofs.fields.f048.20141027.t15z.nc

So to make a time series, use subsequent nowcasts updates strung together sequentially
by update date/time then by n0001-n005 (leave off the last one as it overlaps with
the next set of files)

Similarly append the forecast that is the same update cycle as the most recent nowcast

Confusing? Yes. Run the code and look at the output, the hard work is already done :)

!!!!!Important note: this is for newer ROMS and FVCOM models only. The POM models
still have old file structure with more than one time step per file

'''

def make_server_filelist(model,hour0,start,end=None,test_exist=False):
    '''
    Create a list of model file urls for an aggregated time series based on 
    user specified start and end dates
    
    Args:
        model (string): The COOPS OFS (e.g. NGOFS)
        hour0 (int): The first hour that the model is updated on
            For triangular grid models this is typically 3
            For ROMS models this is typically 0
        start (datetime.date): Look for model output beginning 
            on this date
        end (datetime.date): Look for model output ending before
            this date (if None or > datetime.utcnow() append latest forecast)
        test_exists(bool): Set to True when accessing files from COOPS server
            and want to check existence before operating on them       

    Returns:
        flist (list): List of urls                 
    '''    
    
    flist = []
    stem = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/' + model.upper() + '/MODELS/'
    sdate = datetime.datetime.combine(start,datetime.time(hour0,0))
    if end is None or end > datetime.datetime.utcnow().date():
        edate = datetime.datetime.utcnow()
        append_fc = 1
    else:
        edate = datetime.datetime.combine(end,datetime.time(hour0,0))
        append_fc = 0
        
    while sdate <= edate:
        ym = str(sdate.year) + str(sdate.month).zfill(2)
        ymd = ym + str(sdate.day).zfill(2)
        h = str(sdate.hour).zfill(2)
        fname = stem + ym + '/nos.' + model.lower() + '.fields.n000.' + ymd + '.t' + h + 'z.nc'
        
        agg = make_agg(fname,type='nc')
        flist.extend(agg)
        sdate = sdate + datetime.timedelta(days=.25) #nowcast files are 6 hourly
        
    #check files exist by looking for 404 error
    if test_exist:
        flist = [f for f in flist if test_server_existence(f + '.html')]

    if append_fc:
        last_nc = flist[-1].split('/')[-1].split('n005.')[-1]
        fc_file0 = stem + ym + '/nos.' + model.lower() + '.fields.f000.' + last_nc
        fc_flist = make_agg(fc_file0)
        flist.extend(fc_flist)
        
    return flist
    
def sort_local_files(local_dir,model):
    '''
    Create a filelist for an aggregated time series in local directory 
    
    
         

    Returns:
        flist (list): List of absolute filepaths                  
    '''    
    nc0_files = glob.glob(os.path.join(local_dir,'*n000*'))
    
    flist = []
    for f in nc0_files:
        nc_complete = True #Add nc files if all 6 hours are there
        agg = make_agg(f,'nc')
        for f in agg:
            if not os.path.exists(f):
                nc_complete = False
        if nc_complete:
            flist.extend(agg)
    
    fc0_file = flist[-1].replace('n005','f000')   
    fc_complete = True #Add nc files if all 6 hours are there
    agg = make_agg(fc0_file,'fc')
    for f in agg:
        if not os.path.exists(f):
            fc_complete = False
    if fc_complete:
        flist.extend(agg)

    return flist, nc_complete + fc_complete
    
    
def make_agg(fc_file0,type='fc'):
    if type == 'fc':
        num_files = 48
    elif type == 'nc':
        num_files = 5
        # here we leave off the last file in order to make best time series of nowcast files
        # there is a one hour overlap between adjacent nowcasts
    else:
        print 'Type must be fc or nc'
    a,b = fc_file0.split('000')
    agg = [fc_file0,]
    for h in range(1,num_files+1):
        agg.append(a + str(h).zfill(3) + b)
    return agg

def test_existence(url):
    req = urllib2.Request(url)
    try:
        urllib2.urlopen(req)
        exists = True
    except:
        print 'Not found: ', url
        exists = False
    return exists
    
def test_server_existence(url):
    resp = requests.get(url)
    if resp.status_code == 200:
        exists = True
    else:
        print 'Not found: ', url
        exists = False
    return exists
     
def download_and_save(url,output_dir):
    nc_in = Dataset(url)
    fname = url.split('/')[-1]
    nc_out = Dataset(os.path.join(output_dir,fname),'w')
    
    #Copy dimensions
    for dname, the_dim in nc_in.dimensions.iteritems():
        nc_out.createDimension(dname, len(the_dim))
        
    for var in ['time','u','v']:
        varin = nc_in.variables[var]
        varout = nc_out.createVariable(var, varin.datatype, varin.dimensions)
        varout[:] = varin[:]
        for name in varin.ncattrs():
            value = getattr(varin,name)
            setattr(varout, name, value)

    nc_in.close()
    nc_out.close()
    
if __name__ == "__main__":
    ofs = 'ngofs'
    hour0 = 3
    #sdate = datetime.date.today()-datetime.timedelta(days=14)
    sdate = datetime.date(2014,10,28)
    flist = make_server_filelist(ofs,3,sdate)
    output_dir = 'C:\\Users\\amy.macfadyen\\Documents\\Projects\\goods\\trunk\\static\\ocean_models\\COOPS\\NGOFS'
    for f in flist:
        nc = Dataset(f)
        t = nc.variables['time']
        ts = num2date(t[:],t.units)
        print ts, '...writing'
        download_and_save(f,output_dir)
        
    

#http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NWGOFS/MODELS/201403/nos.nwgofs.fields.nowcast.20140322.t03z.nc.html
#http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NWGOFS/MODELS/201403/nos.nwgofs.fields.n000.20140322.t03z.nc.html