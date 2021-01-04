#!/usr/bin/env python
from __future__ import print_function
import datetime
from netCDF4 import Dataset, num2date
import os, glob
'''
Methods for generating ordered filelist for a time series of CO-OPS data
(Nowcasts + Forecasts) based on user specified start and end dates. If end 
date is unspecified or greater than datetime.utcnow() the latest forecast
will be automatically be appended.

Notes on COOPS naming and aggregations:

Nowcast and forecast files are created four times a day. Output is hourly in 
individual files. So each update generates some number of nowcast files and forecast files
The update cycle time will be the last model output timestep in the nowcast 
files and the first timestep in the forecast files

nos.OFS.fields.nHHH.YYYYMMDD.tCCz.nc or 
glos.OFS.fields.nHHH.YYYYMMDD.tCCz.nc (Great Lakes models)

OFS = model name, eg. "DBOFS"
nHHH = nowcast number, eg. n001-n006 or "nowcast" for old POM models
CC = update hour, e.g. 00,06...

Different models have different update cycles:
FVCOM (non-Great Lakes): 3,9,15,21
ROMS: 0,6,12,18
FVCOM (Great Lakes): 0,6,12,18
Great Lakes POM: 0,6,12,18
Other POM: 5,11,17,23

So to make a time series, use subsequent nowcasts updates strung together sequentially
by update date/time then by individual nowcast files.

The trick is, there is inconsistency between number of nowcast files among models.
Here's the lowdown:
FVCOM: nowcast numbering from 000-->006 with overlap between 6 and 0
ROMS (except GOMOFS): nowcast numbering from 001-->006 with no overlap
GOMOFS ROMS: nowcast 003 and 006 only - ouptut is 3 hourly, no overlap
POM models: one nowcast file "nowcast" in file name, 6 hours, no overlap between nowcasts

Then append the forecast that is the same update cycle as the most recent nowcast

Confusing? Yes. Run the code and look at the output, the hard work is already done :)


'''

    
def make_server_filelist(model,start,end=None,test_exist=False):
    '''
    Create a list of model file urls for an aggregated time series based on 
    user specified start and end dates
    
    Args:
        model (string): The COOPS OFS (e.g. NGOFS)
        
        start (datetime.date): Look for model output beginning 
            on this date
        end (datetime.date): Look for model output ending before
            this date (if None or > datetime.utcnow() append latest forecast,
            it will not be truncated so it may go beyond end date)
        test_exists(bool): Set to True when accessing files from COOPS server
            and want to check existence before operating on them - can take a long time      

    Returns:
        flist (list): List of urls                 
    '''    
    #A bunch of special casing for the variability among models (see docstring)
    
    #Update cycle - every six hours beginning from hour0
    if model.upper() in ['NGOFS','NWGOFS','NEGOFS','CREOFS','SFBOFS']:
        hour0 = 3
    elif model.upper() in ['NYOFS']:
        hour0 = 5
    else: 
        hour0 = 0
    
    #Variability in file names due to model originator
    if model.upper() in ['LOOFS','LSOFS']: #Great Lakes OFS except for LMHOFS/LEOFS which are labeled nos, nothing is consistent
        origin = 'glofs'
    else:
        origin = 'nos'
    
    #And a difference in number of nowcast/forecast files, ugh
    if model.upper() in ['NGOFS','NEGOFS','NWGOFS','CREOFS','SFBOFS','LEOFS','LMHOFS']:       
        nowcast_numbers = range(0,6)
        if model.upper() in ['LEOFS','LMHOFS']:
            forecast_numbers = range(0,121)
        else:
            forecast_numbers = range(0,49)
    elif model.upper() in ['CBOFS','DBOFS','TBOFS','CIOFS']:
        nowcast_numbers = range(1,7)
        forecast_numbers = range(1,49)
    elif model.upper() in ['GOMOFS']:
        nowcast_numbers = [3,6]
        forecast_numbers = range(0,73)
    else:
        nowcast_numbers = 'nowcast'
        forecast_numbers = 'forecast'
        
    #Now make the filelist based on specified dates and specific model info
    flist = []
    stem = 'https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/' + model.upper() + '/MODELS/'
    sdate = datetime.datetime.combine(start,datetime.time(hour0,0))
    if end is None or end > datetime.datetime.utcnow().date() - datetime.timedelta(hours=8):
        edate = datetime.datetime.utcnow() - datetime.timedelta(hours=8)
        append_fc = 1
    else:
        edate = datetime.datetime.combine(end,datetime.time(hour0,0))
        append_fc = 0
 
    while sdate <= edate:
        ym = str(sdate.year) + str(sdate.month).zfill(2)
        ymd = ym + str(sdate.day).zfill(2)
        h = str(sdate.hour).zfill(2)

        fname = stem + ym + '/' + origin + '.' + model.lower() + '.fields.n000.' + ymd + '.t' + h + 'z.nc'
        agg = make_agg(fname,nowcast_numbers)
        flist.extend(agg)

        sdate = sdate + datetime.timedelta(days=.25) #nowcast files are 6 hourly

    if append_fc:
        fname = stem + ym + '/' + origin + '.' + model.lower() + '.fields.f000.' + ymd + '.t' + h + 'z.nc'
        fc_flist = make_agg(fname,forecast_numbers)
        flist.extend(fc_flist)
        
    #check files exist by creating NetCDF Dataset
    if test_exist:
        flist = [f for f in flist if test_existence(f)]
        
    return flist
    
def make_agg(fname,num_files):
 
    a,b = fname.split('000')

    if isinstance(num_files,str):
        agg = [a[:-1] + num_files + b]
    else:
        agg = []
        for h in num_files:
            agg.append(a + str(h).zfill(3) + b)
      
    return agg
    
# def specify_bnd_types(grid,segs,ss_land_nodes=[]):
    # '''
    # The node values were determined by plotting grid, they
    # are not included in the model output
    # Land_bnd_segs are needed to get the boundary right for subset grids only
    # They are obtained by tri_grid remap_bry_nodes method
    # '''
    # if grid.lower() == 'ngofs':
        # ow = list(range(1,180))
    # elif grid.lower() == 'nwgofs':
        # ow = list(range(1,207))
    # elif grid.lower() == 'negofs':
        # ow = list(range(1,139))
    # elif grid.lower() == 'creofs':
        # ow = [68408,68409,68410,68411,68412,68414,68604,68605,68606,68607,68608,68791,68792,68793,68962,68963,68964,68965,69130,69131,69132,69133,69303,69304,69305,69479,69481,69669,69670,69671,69672,69674,69675,69866,69867,69868,69869,69870,70062,70063,70064,70065,70271,70272,70489,70490,70704,70705,70927,70928,71144,71346,71520,71683,71844,72001,72154,72281,72377,72462,72532,72583,72631,72676,72720,72765,72810,72851,72897,72939,72981,73023,73061,73099,73138,73178,73215,73251,73283,73313,73346,73381,73417,73453,73454,73481,73502,73523]
    # elif grid.lower() == 'sfbofs':
        # ow = [1,2,3,4,5,97,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,144,51,52,53,54,55,150,56,57,58,59,60,61,62,63,64,65,66,162,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91]
    # elif grid.lower() == 'gom3':
		# ow = list(range(1,120))
	# else:
        # ow = [1,10000]
        
    # seg_types= []
    
    # if len(ss_land_nodes) > 0: #subset
        # for seg in segs:
            # if seg[0] in ss_land_nodes and seg[1] in ss_land_nodes:
                # seg_types.append(0)
            # else:
                # seg_types.append(1)
    # else:
        # for seg in segs:
            # if seg[0] in ow and seg[1] in ow:
                # seg_types.append(1)
            # else:
                # seg_types.append(0)
            
    # return seg_types
    
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
    
    
def test_existence(url):
    #Also reporting time(s) in file for checking time series is ordered correctly and units same throughout
    try:
        nc = Dataset(url)
        exists = True
        try:
            t = nc.variables['ocean_time']           
        except KeyError:
            t = nc.variables['time']        
        print(num2date(t[:],t.units))
        print(t.units)
        nc.close()
    except:
        print('Not found: ', url)
        exists = False
    return exists
    
def fix_dbofs_mask(grid_file):
    #unclear why the native grid has an error in the rho mask
    nc = Dataset(grid_file,'r+')
    mask_rho = nc.variables['mask_rho']
    mask_rho[71,104:114] = 0  
    nc.close()
    
def download_and_save(url,output_dir):
    nc_in = Dataset(url)
    fname = url.split('/')[-1]
    nc_out = Dataset(os.path.join(output_dir,fname),'w')
    
    #Copy dimensions
    for dname, the_dim in nc_in.dimensions.items():
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
    
def ofs_info(ofs): #This is text to display in GOODS

    ofs = ofs.upper()

    if ofs == 'CREOFS':

        info = '''
        The <a href="http://tidesandcurrents.noaa.gov/ofs/creofs/creofs.html" target="_blank">
        Columbia River Estuary Operational Forecast System (CREOFS)</a> was 
        jointly developed by the <a href="http://www.ohsu.edu/xd/" target="_blank">
        Oregon Health & Science University (OHSU)</a>, 
        the <a href="http://www.nauticalcharts.noaa.gov/" target="_blank">
        NOAA/National Ocean Service's (NOS) Office of Coast Survey </a> and 
        <a href="http://tidesandcurrents.noaa.gov/" target="_blank">
        Center for Operational Oceanographic Products and Services (CO-OPS)</a>, 
        and the <a href="http://mag.ncep.noaa.gov" target="_blank">
        NOAA/National Weather Service's (NWS) National Centers 
        for Environmental Prediction (NCEP) Central Operations (NCO)</a>.
        The CREOFS model domain includes the upper and lower Columbia River and Estuary. 
        
        For detailed model information, visit the NOAA CO-OPS 
        <a href="http://tidesandcurrents.noaa.gov/ofs/creofs/creofs_info.html" target="_blank">
        model information page</a>.
        '''
        
    elif any(x in ofs for x in ['NGOFS','NEGOFS','NWGOFS']):
        
        info = '''
        A <a href="http://tidesandcurrents.noaa.gov/ofs/ngofs/ngofs.html" target="_blank">
        Northern Gulf of Mexico Operational Forecast System (NGOFS)</a>
        including two nested Northeast and Northwest Gulf of Mexico Operational 
        Forecast Systems (NEGOFS/NWGOFS)
        has been developed to serve the maritime user community. 
        NGOFS was developed in a joint project of the 
        <a href="http://www.nauticalcharts.noaa.gov/" target="_blank">
        NOAA/National Ocean Service's (NOS) Office of Coast Survey </a>, 
        the <a href="http://tidesandcurrents.noaa.gov/" target="_blank">
        NOAA/NOS Center for Operational Oceanographic Products and Services (CO-OPS)</a>, 
        the <a href="http://mag.ncep.noaa.gov" target="_blank">
        NOAA/National Weather Service's (NWS) National Centers 
        for Environmental Prediction (NCEP) Central Operations (NCO)</a>, and the 
        <a href="http://fvcom.smast.umassd.edu/" target="_blank">
        University of Massachusetts, Dartmouth </a> using the Finite Volume Coastal Ocean 
        Model (FVCOM). NGOFS generates water level, current, temperature and salinity 
        nowcast and forecast guidance four times per day.
        
        For detailed model information, visit the NOAA CO-OPS 
        <a href="http://tidesandcurrents.noaa.gov/ofs/ngofs/ngofs_info.html" target="_blank">
        model information page.</a>
        '''
        
    elif any(x in ofs for x in ['DBOFS','TBOFS','CBOFS']):
        
        info = '''
        The <a href="http://tidesandcurrents.noaa.gov/ofs/cbofs/cbofs.html" target="_blank">
        Chesapeake Bay Operational Forecast System (CBOFS)</a> was developed by 
        the <a href="http://www.nauticalcharts.noaa.gov/" target="_blank">
        NOAA/National Ocean Service/Office of Coast Survey</a> in a joint project 
        with the <a href="http://tidesandcurrents.noaa.gov/" target="_blank">
        NOAA/NOS/Center for Operational Oceanographic Products and Services 
        (CO-OPS)</a> and the <a href="http://mag.ncep.noaa.gov" target="_blank">
        NOAA/National Weather Service/National Centers for 
        Environmental Prediction (NCEP) Central Operations (NCO)</a> using 
        <a href="http://www.myroms.org/" target="_blank">Rutgers 
        University's Regional Ocean Modeling System (ROMS)</a>. 
        CBOFS generates water level, current, temperature and salinity nowcast 
        and forecast guidance four times per day.
    
        For detailed model information, visit the NOAA CO-OPS 
        <a href="http://tidesandcurrents.noaa.gov/ofs/cbofs/cbofs_info.html" target="_blank">
        model information page.</a>   
        '''
        if ofs == 'DBOFS':
            info = info.replace('cbofs','dbofs')
            info = info.replace('CBOFS','DBOFS')
            info = info.replace('Chesapeake','Delaware')
        elif ofs == 'TBOFS':
            info = info.replace('cbofs','tbofs')
            info = info.replace('CBOFS','TBOFS')
            info = info.replace('Chesapeake','Tampa')
            
    elif ofs == 'SFBOFS':
        
        info = '''
        A <a href="http://tidesandcurrents.noaa.gov/ofs/sfbofs/sfbofs.html" target="_blank">
        San Francisco Bay Operational Forecast System (SFBOFS)</a>
        has been developed to serve the San Francisco Bay maritime communities. 
        SFBOFS was jointly developed by <a href="http://www.nauticalcharts.noaa.gov/" target="_blank">
        NOAA/National Ocean Service's (NOS) Office of Coast Survey </a>, 
        the <a href="http://tidesandcurrents.noaa.gov/" target="_blank">
        NOAA/NOS Center for Operational Oceanographic Products and Services (CO-OPS)</a>, 
        the <a href="http://mag.ncep.noaa.gov" target="_blank">
        NOAA/National Weather Service's (NWS) National Centers 
        for Environmental Prediction (NCEP) Central Operations (NCO)</a>, and the 
        <a href="http://fvcom.smast.umassd.edu/" target="_blank">
        University of Massachusetts, Dartmouth </a> using the Finite Volume Coastal Ocean 
        Model (FVCOM).The NWS and NOS work together to run SFBOFS operationally.
        
        For detailed model information, visit the NOAA CO-OPS 
        <a href="http://tidesandcurrents.noaa.gov/ofs/sfbofs/sfbofs_info.html" target="_blank">
        model information page.</a>   
        '''
        
    elif ofs == 'LEOFS':
        
        info = '''
        
        The upgraded <a href="http://tidesandcurrents.noaa.gov/ofs/leofs/leofs.html" target="_blank">
        Lake Erie Operational Forecast System (LEOFS)</a> was jointly developed by the
        <a href="http://tidesandcurrents.noaa.gov/" target="_blank">
        NOAA/NOS Center for Operational Oceanographic Products and Services (CO-OPS)</a>
        and <a href="http://www.nauticalcharts.noaa.gov/" target="_blank">
        Office of Coast Survey</a>, <a href="http://www.glerl.noaa.gov/" target="_blank">
        the Great Lakes Environmental Research Laboratory (GLERL)</a>, 
        the <a href="http://mag.ncep.noaa.gov" target="_blank">
        NOAA/National Weather Service's (NWS) National Centers 
        for Environmental Prediction (NCEP) Central Operations (NCO)</a>, 
        and the<a href="http://fvcom.smast.umassd.edu/" target="_blank">
        University of Massachusetts, Dartmouth </a> using the Finite Volume Coastal Ocean 
        Model (FVCOM). The NWS and NOS work together to run LEOFS operationally.
        
        For detailed model information, visit the NOAA CO-OPS 
        <a href="http://tidesandcurrents.noaa.gov/ofs/leofs/leofs_info.html" target="_blank">
        model information page.</a>   
        
        '''
        
    else:
  
        return ''
        
        
    return info

    
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
        
    
