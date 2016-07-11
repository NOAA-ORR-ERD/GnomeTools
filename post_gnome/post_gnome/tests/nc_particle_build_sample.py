#!/usr/bin/env python

"""
A script to build a little example netcdf file for trajectories
"""

import numpy as np
import datetime

from post_gnome import nc_particles

def write_sample_file(filename):
    start_time = datetime.datetime(2010, 11, 3, 12)
    #timesteps = [T.start_time + datetime.timedelta(hour=1) for i in range(10)]
    timesteps = [datetime.datetime(2010, 11, 3, 12, 0),
                 datetime.datetime(2010, 11, 3, 12, 30),
                 datetime.datetime(2010, 11, 3, 13, 0),
                 ]

    trajectory = []

    # first timestep: three particles
    trajectory.append( {'longitude': np.array( [ -88.0,
                                           -88.1,
                                           -88.1,
                                           ], dtype=np.float64 ),
                        'latitude': np.array( [ 28.0, 
                                           28.0,
                                           28.1,
                                           ], dtype=np.float64 ),
                        'depth': np.array( [ 0.0, 
                                             0.1,
                                             0.2,
                                             ], dtype=np.float64 ),
                        'mass': np.array( [ 0.01, 
                                            0.005,
                                            0.007,
                                            ], dtype=np.float64 ),
                        'id': np.array( [ 0, 
                                          1,
                                          2,
                                          ], dtype=np.int32 ),
                        'status_code': np.array( [1,
                                                  2,
                                                  3,
                                                  ], dtype=np.int16,)
                        } )
                      
    # second timestep: four particles
    trajectory.append( {'longitude': np.array( [ -88.0,
                                           -88.1,
                                           -88.1,
                                           -87.9,
                                           ], dtype=np.float64 ),
                        'latitude': np.array( [ 28.0, 
                                           28.0,
                                           28.1,
                                           27.9
                                           ], dtype=np.float64 ),
                        'depth': np.array( [ 0.0, 
                                             0.1,
                                             0.2,
                                             0.1
                                             ], dtype=np.float64 ),
                        'mass': np.array( [ 0.01, 
                                            0.005,
                                            0.007,
                                            0.006,
                                            ], dtype=np.float64 ),
                        'id': np.array( [ 0, 
                                          1,
                                          2,
                                          3,
                                          ], dtype=np.int32 ),
                        'status_code': np.array( [1,
                                                  2,
                                                  3,
                                                  4,
                                                  ], dtype=np.int16,)
                        } )

    # third timestep: two particles
    trajectory.append( {'longitude': np.array( [ -88.0,
                                           -88.1,
                                           ], dtype=np.float64 ),
                        'latitude': np.array( [ 28.0, 
                                           28.0,
                                           ], dtype=np.float64 ),
                        'depth': np.array( [ 0.0, 
                                             0.1,
                                             ], dtype=np.float64 ),
                        'mass': np.array( [ 0.01, 
                                            0.005,
                                            ], dtype=np.float64 ),
                        'id': np.array( [ 1, 
                                          3,
                                          ], dtype=np.int32 ),
                        'status_code': np.array( [2,
                                                  3,
                                                  ], dtype=np.int16,)
                        } )
                      
    writer = nc_particles.Writer(filename,
                                 num_timesteps=len(timesteps),
                                 ref_time=timesteps[0],
                                 )
    for i, time in enumerate(timesteps):
        writer.write_timestep(time, trajectory[i])
    writer.close()

if __name__ == "__main__":
    write_sample_file('test_particles.nc')



                          