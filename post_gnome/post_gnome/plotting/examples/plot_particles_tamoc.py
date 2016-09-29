from post_gnome.plotting import geo_plots
reload(geo_plots)
import matplotlib.pyplot as plt
import datetime

plt.clf()

ax = geo_plots.setup_3d(bbox=(-88,-87,27.5,28.5,0,2000)) 
#if bbox not specified, this will use map bounds from bna

#add particles at one time
t0 = datetime.datetime(2016,9,18,1)
#ax = geo_plots.plot_particles(ax,'script_plume.nc',t0,color='b')

#sz = 1000
#sv = 'droplet_diameter'
sz = 4
sv = None

t1 = t0 + datetime.timedelta(hours=48)
ax = geo_plots.plot_particles_3d(ax, 'script_plume.nc',t1,
                                 colormap='plasma',
                                 var = 'droplet_diameter',
                                 drop_size=sz,
                                 drop_scale_var=sv)

ax.legend()

#add initial location
#t = datetime.datetime(2016,9,18,1)
#geo_plots.plot_particles(ax,'script_plume.nc',t,marker='+',markersize=8)

plt.show()