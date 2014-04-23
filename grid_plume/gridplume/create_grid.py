#!/usr/bin/env python

from __future__ import division #Change the / operator to ensure true division throughout (Zelenke).
import numpy as np

def sigma(C, t, t_0, alpha):
    """
    Given some constants and a particular time step in seconds, calculate the
    value of sigma at that time for a particular dimension.
    """
    return ( C * (t + t_0) )**alpha 


def create_grid( x,
                 y ,
                 z ,
                 times,
                 C_x,
                 C_y,
                 C_z,
                 t_0,
                 alpha,
                 min_bucket_count_x,
                 min_bucket_count_y,
                 min_bucket_count_z,
                 max_grids_size,
                 bucket_size,
                 ):
    """
    Create a new concentration grid for subsequent use with the
    calculate_concentration_timestep_agemassvars() function. The extent of this
    grid is buffered out a bit beyond the extent of the given particles.

    :param x: x coordinates of particles in meters
    :type x: array of np.double_t
    :param y: y coordinates of particles in meters
    :type y: array of np.double_t
    :param z: z coordinates of particles in meters (positive up)
    :type z: array of np.double_t
    :param times: times in seconds since first time step
    :type times: array of np.double_t
    :param sigma_x_0: initial x dimension of the plume in meters
    :type sigma_x_0: np.double_t
    :param sigma_y_0: initial y dimension of the plume in meters
    :type sigma_y_0: np.double_t
    :param sigma_z_0: initial z dimension of the plume in meters
    :type sigma_z_0: np.double_t
    :param alpha: speed of spreading in centimeters per second
    :type alpha: np.double_t
    :param min_bucket_count_x: minimum number of buckets in the x direction
    :type min_bucket_count_x: np.uint32_t
    :param min_bucket_count_y: minimum number of buckets in the y direction
    :type min_bucket_count_y: np.uint32_t
    :param min_bucket_count_z: minimum number of buckets in the z direction
    :type min_bucket_count_z: np.uint32_t
    :param max_grids_size: maximum size of all grids for all time steps in
                           bytes
    :type max_grids_size: np.uint64_t
    :param bucket_size: bucket size in bytes for purposes of calculating the
                        size of all grids
    :type bucket_size: np.uint32_t
    :return: newly created concentration grid for a single time step, zeroed
    :rtype: ZxYxX array of np.double_t
    :return: tuple of (
                minimum x value in meters,
                maximum x value in meters,
                minimum y value in meters,
                maximum y value in meters,
                minimum z value in meters,
                maximum z value in meters,
             )
    :rtype: tuple of (
                np.double_t,
                np.double_t,
                np.double_t,
                np.double_t,
                np.double_t,
                np.double_t,
            )
    """    
    time_count = times.shape[ 0 ]

    # Calculate the distance between buckets based on the sigma for the last
    # time step.
    bucket_x_delta = sigma( C_x, times[ time_count - 1 ], t_0, alpha )
    bucket_y_delta = sigma( C_y, times[ time_count - 1 ], t_0, alpha )
    bucket_z_delta = sigma( C_z, times[ time_count - 1 ], t_0, alpha )

    # Buffer extent out by the area of effect (3*sigma_final) on each side.
    grid_x_min = np.nanmin( x ) - 3 * bucket_x_delta
    grid_x_max = np.nanmax( x ) + 3 * bucket_x_delta
    grid_y_min = np.nanmin( y ) - 3 * bucket_y_delta
    grid_y_max = np.nanmax( y ) + 3 * bucket_y_delta

    # don't create a grid above water level: z <= 0
    grid_z_min = max( 0, np.nanmin( z ) - 3 * bucket_z_delta )
    grid_z_max = np.nanmax( z ) + 3 * bucket_z_delta

    bucket_count_z = int(max( ( grid_z_max - grid_z_min ) / bucket_z_delta + 0.5,
                              min_bucket_count_z,
                              ))

    bucket_count_y = int(max( ( grid_y_max - grid_y_min ) / bucket_y_delta + 0.5,
                              min_bucket_count_y,
                              ))
    
    bucket_count_x = int(max( ( grid_x_max - grid_x_min ) / bucket_x_delta + 0.5,
                              min_bucket_count_x,
                              ))

    grids_size = (time_count * bucket_count_z * 
                  bucket_count_y * bucket_count_x *
                  bucket_size)
    
    #cdef np.float64_t factor

    # If necessary, scale the grid dimensions down to fit within the maximum
    # requested grids size.
    if grids_size > max_grids_size:
        # To the 1/3 power because there are three dimensions being scaled.
        factor = ( max_grids_size / grids_size ) ** ( 1 / 3.0 )
        bucket_count_z = ( bucket_count_z * factor )
        bucket_count_y = ( bucket_count_y * factor )
        bucket_count_x = ( bucket_count_x * factor )

    grid = np.zeros(( bucket_count_z,
                      bucket_count_y,
                      bucket_count_x,
                      ),
                    np.float64,
                    )

    return ( grid,
             grid_x_min,
             grid_x_max,
             grid_y_min,
             grid_y_max,
             grid_z_min,
             grid_z_max,
             )
