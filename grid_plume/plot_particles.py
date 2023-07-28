#!/usr/bin/env python

"""
Module that plots nc_particle_files
Dan Helfman
"""

from __future__ import division #Change the / operator to ensure true division throughout (Zelenke).
from gridplume import nc_particles as ncp
print ncp.__file__
#from gridplume import particles #Commented-out import of unused module (Zelenke).
from pylab import *
import matplotlib.pyplot as plt
import netCDF4 as nc
#import numpy as np #Commented-out import of unused module (Zelenke).
#import gdata  #Commented-out import of unused module (Zelenke).
#import googlemaps #Commented-out import of unused module (Zelenke).

 
 #experimental plot for GNOME2 
    #def gnome_plot():
        
        

#bring in relevant data for plotting
data = nc.Dataset('uh.nc', 'r')

#instantiate particle from data
particle = ncp.nc_particle_file(data)

#Tk Window Plot (quick)
"""quick plots: user clicks exit button to view next plot"""
def static_plots():
    for i in range(10):
        datadict = particle.get_timestep(i)
        plt.plot(datadict['longitude'],datadict['latitude'], '.')
        plt.show()

#Tk Window Plot (test)
"""pylab neutral animation plot"""
def animated_plot_test(): 
    ion()
    #dataaxisdict = particle.get_all_timesteps #Commented-out variable, as it's never used (Zelenke).
    #axisscale =[[np.amin(datadict['longitude']),np.amax(datadict['longitude']),np.amin(datadict['latitude']),np.max(datadict['latitude']]             
     
     
    for i in range(1000):
        #graph particles        
        datadict = particle.get_timestep(i)
        line, = plot(datadict['longitude'],datadict['latitude'], '.')
        #set labels, axies, update data
        xlabel('longitude' r"$(\circ)$")
        ylabel('latitude'  r"$(\circ)$")
        title(r"$UXO$")
        #ylim(amin(datadict['latitude']), amax(datadict['latitude']))
        #xlim(amin(datadict['longitude']), amax(datadict['longitude']))
        #line.set_xdata(datadict['latitude'])        
        #line.set_ydata(datadict['longitude'])
        #ylim(21.5,22)
        #xlim(-158,-158.3)
        v = [-158.275,-158.18,21.40,21.568]       
        axis(v)     
        draw()  
        clf()
#animated_plot_test()                       

#Tk Window Plot (production animated plot)
def animated_plot(*args):
    uxoplot.setydata()
    ax.draw_animated
    
    #Tk Window Plot (test)
"""pylab neutral animation plot"""
def animated_plot_test_2(): 
    ion()
    #axisscale =[[np.amin(datadict['longitude']),np.amax(datadict['longitude']),np.amin(datadict['latitude']),np.max(datadict['latitude']]             
    data = particle.get_all_timesteps(['longitude','latitude','depth','mass'])
    numtimes = len(particle.times)
    lon = data['longitude'].reshape((numtimes, -1))
    lat = data['latitude'].reshape((numtimes,-1))
  
    for i in range(1000):
        #graph particles        
        line, = plot(lon[i],lat[i], '.')
        #set labels, axies, update data
        xlabel('longitude' r"$(\circ)$")
        ylabel('latitude'  r"$(\circ)$")
        title(r"$UXO$")
        #ylim(amin(datadict['latitude']), amax(datadict['latitude']))
        #xlim(amin(datadict['longitude']), amax(datadict['longitude']))
        #line.set_xdata(datadict['latitude'])        
        #line.set_ydata(datadict['longitude'])
        #ylim(21.5,22)
        #xlim(-158,-158.3)
        v = [-158.275,-158.18,21.40,21.568]       
        axis(v)     
        draw()  
        clf()
animated_plot_test_2()
