from post_gnome.plotting import geo_plots
reload(geo_plots)
import matplotlib.pyplot as plt
import datetime

plt.clf()

ax = geo_plots.add_map(bbox=(-89,-85,27.5,31)) 
#if bbox not specified, this will use map bounds from bna

#add particles at one time
t = datetime.datetime(2016,9,21,1)

filename = 'script_plume.nc'
ax = geo_plots.contour_particles(ax,filename,t)
#ax = geo_plots.plot_particles(ax,filename,t,depth=0,color='b')


#add initial location


plt.show()