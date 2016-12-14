from post_gnome.plotting import geo_plots
reload(geo_plots)

import matplotlib.pyplot as plt
import datetime
from cartopy.io.shapereader import Reader
import os
import numpy as np

plt.clf()

ax = geo_plots.setup_3d(bbox=(-88, -87, 27.5, 28.5, 0, 2000))
# if bbox not specified, this will use map bounds from bna

# add bathymetry
bathy_file = os.path.join('gom_bathy','Bathymetry.shp')
bathy = Reader(bathy_file)

for rec,geo in zip(bathy.records(), bathy.geometries()):
    if rec.attributes['DEPTH_METR'] in ['1000m']:
        print rec.attributes['DEPTH_METR']
        for g in geo:
            x = g.xy[0]
            y = g.xy[1]
            ax.plot3D(x, y, 1000 * np.ones_like(x), 'k')

ax.set_xlim(-87.5, -86.9)
ax.set_ylim(27.6, 28.7)

# add particles at one time
t0 = datetime.datetime(2016, 9, 18, 1)
# ax = geo_plots.plot_particles(ax,'script_plume.nc',t0,color='b')

# sz = 1000
# sv = 'droplet_diameter'
sz = 4
sv = None

t1 = t0 + datetime.timedelta(hours=48)
ax = geo_plots.plot_particles_3d(ax, 'gulf_tamoc.nc', t1,
                                     colormap='plasma',
                                     var='droplet_diameter',
                                     drop_size=sz,
                                     drop_scale_var=sv)

ax.legend()

# add initial location
# t = datetime.datetime(2016,9,18,1)
# geo_plots.plot_particles(ax,'script_plume.nc',t,marker='+',markersize=8)

plt.show()
