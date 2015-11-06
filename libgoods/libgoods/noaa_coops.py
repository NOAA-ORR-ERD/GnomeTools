#!/usr/bin/env python
import datetime
import urllib2, requests
from netCDF4 import Dataset
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
def specify_bnd_types(grid,segs,ss_land_nodes=[]):
    '''
    The node values were determined by plotting grid, they
    are not included in the model output
    Land_bnd_segs are needed to get the boundary right for subset grids only
    They are obtained by tri_grid remap_bry_nodes method
    '''
    if grid.lower() == 'ngofs':
        ow = range(1,180)
    elif grid.lower() == 'nwgofs':
        ow = range(1,207)
    elif grid.lower() == 'negofs':
        ow = range(1,139)
    elif grid.lower() == 'creofs':
        ow = [68408,68409,68410,68411,68412,68414,68604,68605,68606,68607,68608,68791,68792,68793,68962,68963,68964,68965,69130,69131,69132,69133,69303,69304,69305,69479,69481,69669,69670,69671,69672,69674,69675,69866,69867,69868,69869,69870,70062,70063,70064,70065,70271,70272,70489,70490,70704,70705,70927,70928,71144,71346,71520,71683,71844,72001,72154,72281,72377,72462,72532,72583,72631,72676,72720,72765,72810,72851,72897,72939,72981,73023,73061,73099,73138,73178,73215,73251,73283,73313,73346,73381,73417,73453,73454,73481,73502,73523]
    elif grid.lower() == 'sfbofs':
        ow = [1,2,3,4,5,97,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,144,51,52,53,54,55,150,56,57,58,59,60,61,62,63,64,65,66,162,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91]
    
    seg_types= []
    
    if len(ss_land_nodes) > 0: #subset
        for seg in segs:
            if seg[0] in ss_land_nodes and seg[1] in ss_land_nodes:
                seg_types.append(0)
            else:
                seg_types.append(1)
    else:
        for seg in segs:
            if seg[0] in ow and seg[1] in ow:
                seg_types.append(1)
            else:
                seg_types.append(0)
            
    return seg_types
    
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
            this date (if None or > datetime.utcnow() append latest forecast,
            it will not be truncated so it may go beyond end date)
        test_exists(bool): Set to True when accessing files from COOPS server
            and want to check existence before operating on them       

    Returns:
        flist (list): List of urls                 
    '''    
    
    flist = []
    stem = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/' + model.upper() + '/MODELS/'
    sdate = datetime.datetime.combine(start,datetime.time(hour0,0))
    if end is None or end > datetime.datetime.utcnow().date() - datetime.timedelta(hours=6):
        edate = datetime.datetime.utcnow() - datetime.timedelta(hours=6)
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
    
    

    
#if __name__ == "__main__":
#    ofs = 'ngofs'
#    hour0 = 3
#    #sdate = datetime.date.today()-datetime.timedelta(days=14)
#    sdate = datetime.date(2014,10,28)
#    flist = make_server_filelist(ofs,3,sdate)
#    output_dir = 'C:\\Users\\amy.macfadyen\\Documents\\Projects\\goods\\trunk\\static\\ocean_models\\COOPS\\NGOFS'
#    for f in flist:
#        nc = Dataset(f)
#        t = nc.variables['time']
#        ts = num2date(t[:],t.units)
#        print ts, '...writing'
#        download_and_save(f,output_dir)
        
    