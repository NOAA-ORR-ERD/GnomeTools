import numpy as np
#import netCDF4 as nc #Commented-out import of unused module (Zelenke).
import nc_particles as ncp


def load_particles(particle_dataset):
    """
    Load a netCDF4 dataset of particles and return it. Close the dataset when
    done.

    :param dataset: open readable NetCDF dataset
    :type dataset: netCDF.Dataset
    :return: (
                x coordinates of particles (longitude),
                y coordinates of particles (latitude),
                z coordinates of particles in meters,
                ages of particles in seconds,
                times converted to datetime objects,
                original time units,
                original time base date,
            )
    :rtype: tuple of (
                TxN array,
                TxN array,
                TxN array,
                TxN array,
                T array of datetime objects,
                str,
                array,
            )

    Where T is the number of time steps and N is the number of particles.
  
"""



    ncparticle = ncp.nc_particle_file(particle_dataset)
    
    
    
    NumLE = ncparticle.particle_count[:]
    times = ncparticle.times[:]    
    NumTimeSteps = len(times)
    data = ncparticle.get_all_timesteps(['longitude','latitude','depth','mass','age'])
    #print data.keys()
    
    
    #Check whether the number of particles remains constant in each time-step
    #(instaneous release) or changes over time (continuous release) (Zelenke).
    if len( np.unique(NumLE) ) == 1:
        
        #Reshape LE time-series vectors into rectangular matricies (NumTimeSteps x NumLE).
        lon = data['longitude'].reshape( (NumTimeSteps, -1) )
        lat = data['latitude'].reshape( (NumTimeSteps, -1) )
        depth = data['depth'].reshape( (NumTimeSteps, -1) )
        mass = data['mass'].reshape( (NumTimeSteps, -1) )
        age = data['age'].reshape( (NumTimeSteps, -1) )
        
    else:
        
        #Put varying number of LEs in each time-step into uniformly sized vectors
        #to allow for resizing into rectangular matricies (NumTimeSteps x NumLE).
        
        MaxNumLE = np.max(NumLE)
        lons = np.ones( MaxNumLE*NumTimeSteps, dtype = data['longitude'].dtype ) * np.array(np.NaN, dtype = data['longitude'].dtype ) #Preallocate vector for subsequent loop.
        lats = np.ones( MaxNumLE*NumTimeSteps, dtype = data['latitude'].dtype ) * np.array(np.NaN, dtype = data['latitude'].dtype )
        depths = np.ones( MaxNumLE*NumTimeSteps, dtype = data['depth'].dtype ) * np.array(np.NaN, dtype = data['depth'].dtype )
        masses = np.ones( MaxNumLE*NumTimeSteps, dtype = data['mass'].dtype ) * np.array(np.NaN, dtype = data['mass'].dtype )
        ages = np.copy(masses)
        
        CumSumLE = np.append(0,np.cumsum(NumLE))
        for i in range(NumTimeSteps):
            #Put however many LEs there were in each time-step into that time-step's (uniformly sized) portion of the vector.
            #Leave the rest of the chunk populated with NaNs, if there aren't enough LEs in that particular time-step to fill-up the chunk.
            lons[ i*MaxNumLE : (i*MaxNumLE)+NumLE[i] ] = data['longitude'][ CumSumLE[i] : CumSumLE[i+1] ]
            lats[ i*MaxNumLE : (i*MaxNumLE)+NumLE[i] ] = data['latitude'][ CumSumLE[i] : CumSumLE[i+1] ]
            depths[ i*MaxNumLE : (i*MaxNumLE)+NumLE[i] ] = data['depth'][ CumSumLE[i] : CumSumLE[i+1] ]
            masses[ i*MaxNumLE : (i*MaxNumLE)+NumLE[i] ] = data['mass'][ CumSumLE[i] : CumSumLE[i+1] ]
            ages[ i*MaxNumLE : (i*MaxNumLE)+NumLE[i] ] = data['age'][ CumSumLE[i] : CumSumLE[i+1] ]
        #Reshape (NumTimeSteps x NumLE).
        lon=lons.reshape( (NumTimeSteps, -1) )
        lat = lats.reshape( (NumTimeSteps, -1) )
        depth = depths.reshape( (NumTimeSteps, -1) )
        mass = masses.reshape( (NumTimeSteps, -1) )
        age = ages.reshape( (NumTimeSteps, -1) )
    
    
    
    time_units = ncparticle.time_units
    mass_units = ncparticle.mass_units #Included units of mass specified in the NetCDF file as an output variable (Zelenke).

    return (
        lon,
        lat,
        depth,
        age,
        times,
        time_units,
        mass,
        mass_units,
    )
