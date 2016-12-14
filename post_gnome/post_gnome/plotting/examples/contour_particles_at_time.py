from post_gnome.plotting import geo_plots
reload(geo_plots)
import matplotlib.pyplot as plt
import datetime
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import cartopy.feature as cfeature
import os

plt.clf()

#ax = geo_plots.add_map(bbox=(-89,-85,27.5,31))
ax = geo_plots.add_map(bbox=(-90, -84, 27.5, 31))
#land_10m = cfeature.NaturalEarthFeature('physical','land','50m',\
#    edgecolor='face',facecolor='0.75')
##ax.add_feature(cfeature.LAND, facecolor='0.75')
#ax.add_feature(land_10m)
#if bbox not specified, this will use map bounds from bna

bathy_file = os.path.join('gom_bathy','Bathymetry.shp')
bathy = Reader(bathy_file)

for rec,geo in zip(bathy.records(),bathy.geometries()):
    if rec.attributes['DEPTH_METR'] in ['100m','500m','1000m']:
        shape_feature = ShapelyFeature(geo,ccrs.PlateCarree(), facecolor='none')
        ax.add_feature(shape_feature)

#add particles at one time
t = datetime.datetime(2016, 9, 21, 1)

filename = 'gulf_tamoc.nc'
ax = geo_plots.contour_particles(ax,filename,t,levels=[0.1, 0.4, 0.6, 1])
#ax = geo_plots.plot_particles(ax,filename,t,depth=0,color='b')


#add initial location


plt.show()