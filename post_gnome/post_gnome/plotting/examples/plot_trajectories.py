from post_gnome.plotting import geo_plots
reload(geo_plots)
import matplotlib.pyplot as plt

plt.clf()

ax = geo_plots.add_map(bbox=(-125.2,-124.2,47.7,48.3), bna='coast_wa.bna') 
#if bbox not specified, this will use map bounds from bna

geo_plots.plot_all_trajectories(ax,'WA_particles.nc',addmarker=True)
geo_plots.plot_single_trajectory(ax,'WA_particles.nc',1,color='r')

plt.show()