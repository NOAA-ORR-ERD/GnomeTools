import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from post_gnome import nc_particles
from netCDF4 import Dataset, num2date
import numpy as np

def add_map(bbox=None,bna=None):
    
    print 'Using mercator projection' #might want to extend to include other options?
    ax = plt.axes(projection=ccrs.Mercator())
    
    if bna is not None:
        #load bna and determine bounding box from bna map bounds
        #no informative error messaging if invalid bna
        coast_polys = {}
        
        with open(bna) as f:
            while True:
                try:
                    id, type, num_pts = f.readline().split(',')
                    points = np.zeros((int(num_pts), 2))
                    for i in range(int(num_pts)):
                        points[i,:] = [float(j) for j in f.readline().split(',')]
                    coast_polys[id] = points
                except ValueError:
                    if len(coast_polys.keys()) == 0:
                        print 'bna did not load correctly'
                    break
                    
        if coast_polys.has_key('"Map Bounds"'):
            mb = coast_polys.pop('"Map Bounds"') #need to get rid of this regardless for plotting
            if bbox is None:
                print 'Using map bounds from bna file'
                x0 = mb[:,0].min()
                x1 = mb[:,0].max()
                y0 = mb[:,1].min()
                y1 = mb[:,1].max()
                bbox = (x0,x1,y0,y1)
                print 'bbox:', bbox
    
    if bbox is None:
            bbox=(-180,180,-80,80)  
    ax.set_extent(bbox)
    
    if len(coast_polys.keys()) > 0:
        for poly in coast_polys.itervalues():
            ax.plot(poly[:,0],poly[:,1],'k',transform=ccrs.Geodetic())
    else:
        ax.coastlines(resolution='10m')
    
    ax.gridlines(draw_labels=True)
    
    return ax

def plot_particles(ax,filename,t,color='k',marker='.',markersize=4):
    '''
    plot all LEs at one time step
    ax: (matplotlib.axes object) the map on which the LEs will be plotted
    filename: (str) complete path filename of particle file
    t: (datetime obj) closest time to this will be plottted
    '''
    
    particles = nc_particles.Reader(filename)
    times = particles.times
    dt = [np.abs(((output_t - t).total_seconds())/3600) for output_t in times]
    tidx = dt.index(min(dt))
    try:
        TheData = particles.get_timestep(tidx,variables=['latitude','longitude','status_codes'])
    except: #GUI GNOME < 1.3.10
        TheData = particles.get_timestep(tidx,variables=['latitude','longitude','status'])
        TheData['status_codes'] = TheData['status']
    
    status = TheData['status_codes']
    label = t.isoformat()
    for sc in [2,3]:
        if sc==3:
            marker='x'
            label=None
        pid = np.where(status==sc)[0]
        if len(pid) > 0:
            ax.scatter(TheData['longitude'][pid],TheData['latitude'][pid],transform=ccrs.Geodetic(),\
                color=color,marker=marker,s=markersize,label=label)

    print 'Closest time found: ', times[tidx]
    
    return ax
    
def plot_single_trajectory(ax,filename,particle_id,color='k',addmarker=True,marker='.',markersize='4'):
    '''
    plot single particle trajectory
    ax: (matplotlib.axes object) the map on which the LEs will be plotted
    filename: (str) complete path filename of particle file
    particle_id: (int) particle id
    status_code: 2 for floating only, 3 for beached only
    '''
    
    particles = nc_particles.Reader(filename)
    try:
        TheData = particles.get_individual_trajectory(particle_id,variables=['latitude','longitude','status_codes'])
    except: #GUI GNOME < 1.3.10
        TheData = particles.get_individual_trajectory(particle_id,variables=['latitude','longitude','status'])
        TheData['status_codes'] = TheData['status']
    
    ax.plot(TheData['longitude'],TheData['latitude'],transform=ccrs.Geodetic(),color=color)
    if addmarker:
        status = TheData['status_codes']
        for sc in [2,3]:
            if sc==3:
                marker='x'
            pid = np.where(status==sc)[0]
            if len(pid) > 0:
                ax.plot(TheData['longitude'][pid],TheData['latitude'][pid],transform=ccrs.Geodetic(),\
                    color=color,marker=marker, markersize=markersize,linestyle='None')
    
    return ax
        
def plot_all_trajectories(ax,filename,color='slategray',addmarker=False,marker='.',markersize='4'):
    '''
    plot particle trajectories by ids
    ax: (matplotlib.axes object) the map on which the LEs will be plotted
    filename: (str) complete path filename of particle file
    '''
    
    particles = nc_particles.Reader(filename)
    try:
        TheData = particles.get_all_timesteps(variables=['latitude','longitude','id','status_codes'])
    except: #GUI GNOME < 1.3.10
        TheData = particles.get_all_timesteps(variables=['latitude','longitude','id','status'])
        TheData['status_codes'] = TheData['status']
    
    id = np.array(TheData['id'])
    lon = np.array(TheData['longitude'])
    lat = np.array(TheData['latitude'])
    status = np.array(TheData['status'])
    
    pids = np.unique(id)
    
    for pid in pids:
        x,y = np.where(id==pid)
    
        ax.plot(lon[x,y],lat[x,y],transform=ccrs.Geodetic(),color=color)
    
        if addmarker:
            le_marker=marker
            for sc in [2,3]:
                if sc==3:
                    le_marker='x'
                sid = np.where(status[x,y]==sc)[0]
                if len(sid) > 0:
                    ax.plot(lon[x,y][sid],lat[x,y][sid],transform=ccrs.Geodetic(),\
                        color=color,marker=le_marker, markersize=markersize,linestyle='None')                    
    
    return ax

def add_vectors(ax,filename,t,bbox=None,tvar='time',lonvar='lon',latvar='lat',uvar='water_u',vvar='water_v'):
    '''
    plot vectors from netCDF file - this is pretty much just a start that needs customizing
    to be useful in the general case
    
    ax: (matplotlib.axes object) the map on which the LEs will be plotted
    filename: (str) complete path filename of netCDF file
    t: (datetime obj) closest time to this will be plottted
    '''
    
    nc = Dataset(filename)
    nc_t = nc.variables[tvar]
    nc_dt = num2date(nc_t[:],nc_t.units)
    d = [np.abs(((output_t - t).total_seconds())/3600) for output_t in nc_dt]
    tidx = d.index(min(d))
    
    print 'Plotting vectors at: ', nc_dt[tidx]
    lon = nc.variables[lonvar][:]
    lat = nc.variables[latvar][:]
    u = nc.variables[uvar][tidx,:] #todo: add check for 3d
    v = nc.variables[vvar][tidx,:]
    ax.quiver(lon,lat,u,v,scale=2,transform=ccrs.PlateCarree())
    
    
    return ax
