from __future__ import division #Change the / operator to ensure true division throughout (Zelenke).
import cython
import numpy as np
cimport numpy as np


#This is a version that allows full flexibility:
#    
#    sigma = C * (t-t_0)^alpha
#    
#    alpha can be anything -- 0.5 for "classic" diffusion, 3/2 for the theoretical solution, etc.
#    
#    for "classic diffusion: alpha=0.5 C = 2D.
#
    


cdef extern from "math.h" nogil:
    double exp( double x )


@cython.boundscheck( False )
@cython.cdivision( True )
def create_grid(
    np.ndarray[ np.double_t, ndim = 2 ] x not None,
    np.ndarray[ np.double_t, ndim = 2 ] y not None,
    np.ndarray[ np.double_t, ndim = 2 ] z not None,
    np.ndarray[ np.double_t ] times not None,
    np.double_t C_x,
    np.double_t C_y,
    np.double_t C_z,
    np.double_t t_0,
    np.double_t alpha,
    np.uint32_t min_bucket_count_x,
    np.uint32_t min_bucket_count_y,
    np.uint32_t min_bucket_count_z,
    np.uint64_t max_grids_size,
    np.uint32_t bucket_size,
):
    """
    Create a new concentration grid for subsequent use with the
    calculate_concentration() function. The extent of this grid is buffered
    out a bit beyond the extent of the given particles.

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
    cdef np.uint32_t time_count = times.shape[ 0 ]

    # Calculate the distance between buckets based on the sigma for the last
    # time step.
    cdef np.double_t bucket_x_delta = sigma(
        C_x, times[ time_count - 1 ], t_0, alpha
    )
    cdef np.double_t bucket_y_delta = sigma(
        C_y, times[ time_count - 1 ], t_0, alpha
    )
    cdef np.double_t bucket_z_delta = sigma(
        C_z, times[ time_count - 1 ], t_0, alpha
    )


    # Buffer extent out by the area of effect (3*sigma_final) on each side.
    cdef np.double_t grid_x_min = np.nanmin( x ) - 3 * bucket_x_delta
    cdef np.double_t grid_x_max = np.nanmax( x ) + 3 * bucket_x_delta
    cdef np.double_t grid_y_min = np.nanmin( y ) - 3 * bucket_y_delta
    cdef np.double_t grid_y_max = np.nanmax( y ) + 3 * bucket_y_delta
    cdef np.double_t grid_z_min = np.nanmin( z ) - 3 * bucket_z_delta
    cdef np.double_t grid_z_max = np.nanmax( z ) + 3 * bucket_z_delta

    cdef np.uint32_t bucket_count_z = int(max(
        ( grid_z_max - grid_z_min ) / bucket_z_delta + 0.5,
        min_bucket_count_z,
        ))

    cdef np.uint32_t bucket_count_y = int(max(
        ( grid_y_max - grid_y_min ) / bucket_y_delta + 0.5,
        min_bucket_count_y,
        ))
    
    cdef np.uint32_t bucket_count_x = int(max(
        ( grid_x_max - grid_x_min ) / bucket_x_delta + 0.5,
        min_bucket_count_x,
        ))

    cdef np.uint64_t grids_size = \
        <np.uint64_t>time_count * <np.uint64_t>bucket_count_z * \
        <np.uint64_t>bucket_count_y * <np.uint64_t>bucket_count_x * \
        <np.uint64_t>bucket_size
    cdef np.float64_t factor

    # If necessary, scale the grid dimensions down to fit within the maximum
    # requested grids size.
    if grids_size > <np.uint64_t>max_grids_size:
        # To the 1/3 power because there are three dimensions being scaled.
        factor = ( <np.float64_t>max_grids_size / <np.float64_t>grids_size ) \
                 ** <np.float64_t>( 1 / 3.0 )
        bucket_count_z = <np.uint32_t>( <np.float64_t>bucket_count_z * factor )
        bucket_count_y = <np.uint32_t>( <np.float64_t>bucket_count_y * factor )
        bucket_count_x = <np.uint32_t>( <np.float64_t>bucket_count_x * factor )

    cdef np.ndarray[ np.double_t, ndim = 3 ] grid = \
        np.zeros(
            (
                bucket_count_z,
                bucket_count_y,
                bucket_count_x,
            ),
            np.double,
        )

    return (
        grid,
        grid_x_min,
        grid_x_max,
        grid_y_min,
        grid_y_max,
        grid_z_min,
        grid_z_max,
    )



#==============================================================================
# Commented-out original function as it's no longer used.  Retained for
# reference.  (Zelenke)
# 
# @cython.boundscheck( False )
# @cython.cdivision( True )
# def calculate_concentration(
#     np.ndarray[ np.double_t, ndim = 3 ] grid,
#     file_grid,
#     np.ndarray[ np.double_t, ndim = 2 ] x not None,
#     np.ndarray[ np.double_t, ndim = 2 ] y not None,
#     np.ndarray[ np.double_t, ndim = 2 ] z not None,
#     np.double_t grid_x_min,
#     np.double_t grid_x_max,
#     np.double_t grid_y_min,
#     np.double_t grid_y_max,
#     np.double_t grid_z_min,
#     np.double_t grid_z_max,
#     np.double_t mass,
#     np.ndarray[ np.double_t ] times not None,
#     np.double_t C_x,
#     np.double_t C_y,
#     np.double_t C_z,
#     np.double_t t_0,
#     np.double_t alpha,
#     np.double_t conversion_factor,
# ):
#     """
#     Calculate a concentration grid based on the given particles and save it to
#     file as each time step's grid is calculated.
# 
#     :param grid: concentration grid for a single time step, zeroed
#     :type grid: ZxYxZ array of np.double_t
#     :param file_grid: writable NetCDF4 concentration grid variable object
#     :type file_grid: netCDF4.Variable of TxZxYxZ
#     :param x: x coordinates of particles in meters
#     :type x: array of np.double_t
#     :param y: y coordinates of particles in meters
#     :type y: array of np.double_t
#     :param z: z coordinates of particles in meters (positive up)
#     :type z: array of np.double_t
#     :param grid_x_min: minimum grid x value in meters
#     :type grid_x_min: np.double_t
#     :param grid_x_max: maximum grid x value in meters
#     :type grid_x_max: np.double_t
#     :param grid_y_min: minimum grid y value in meters
#     :type grid_y_min: np.double_t
#     :param grid_y_max: maximum grid y value in meters
#     :type grid_y_max: np.double_t
#     :param grid_z_min: minimum grid z value in meters
#     :type grid_z_min: np.double_t
#     :param grid_z_max: maximum grid z value in meters
#     :type grid_z_max: np.double_t
#     :param mass: mass of each particle in kg
#     :type mass: np.double_t
#     :param times: times in seconds since first time step
#     :type times: array of np.double_t
#     :param sigma_x_0: initial x dimension of the plume in meters
#     :type sigma_x_0: np.double_t
#     :param sigma_y_0: initial y dimension of the plume in meters
#     :type sigma_y_0: np.double_t
#     :param sigma_z_0: initial z dimension of the plume in meters
#     :type sigma_z_0: np.double_t
#     :param alpha: speed of spreading in centimeters per second
#     :type alpha: np.double_t
#     :param conversion_factor: factor to multiply each concentration value by
#     :type conversion_factor: np.double_t
#     """
#     cdef np.uint32_t time_count = times.shape[ 0 ]
# 
# #==============================================================================
# #     Subtracted 1 from the denominator in the calculation of each of the three
# #     delta variables below so that grid dimensions calculated by filter.pyx
# #     would match those calculated by grid.py (Zelenke).
# #==============================================================================
#     cdef np.double_t bucket_x_delta = \
#         ( grid_x_max - grid_x_min ) / ( grid.shape[ 2 ] -1.0 )
#     cdef np.double_t bucket_y_delta = \
#         ( grid_y_max - grid_y_min ) / ( grid.shape[ 1 ] - 1.0 )
#     cdef np.double_t bucket_z_delta = \
#         ( grid_z_max - grid_z_min ) / ( grid.shape[ 0 ] - 1.0 )
# 
#     # Make one density grid per time step.
#     cdef np.uint32_t time_index, particle_index
#     cdef np.double_t sigma_x, sigma_y, sigma_z
#     cdef np.double_t effect_distance_x, effect_distance_y, effect_distance_z
#     cdef np.double_t initial_concentration
#     cdef np.uint32_t x_index_min, x_index_max
#     cdef np.uint32_t y_index_min, y_index_max
#     cdef np.uint32_t z_index_min, z_index_max
#     cdef np.uint32_t x_index, y_index, z_index
#     cdef np.double_t particle_x, particle_y, particle_z
#     cdef np.double_t bucket_x, bucket_y, bucket_z
#     cdef np.double_t bucket_y_min, bucket_z_min
# 
#     for time_index in range( time_count ):
#         # Clear the grid at the start of each time step.
#         if time_index != 0:
#             grid.fill( 0 )
# 
#         # Calculate how far a particle at this time step can affect the grid.
#         sigma_x = sigma(
#             C_x, times[ time_index ], t_0, alpha
#         )
#         sigma_y = sigma(
#             C_y, times[ time_index ], t_0, alpha
#         )
#         sigma_z = sigma(
#             C_z, times[ time_index ], t_0, alpha
#         )
# 
#         effect_distance_x = 3 * sigma_x
#         effect_distance_y = 3 * sigma_y
#         effect_distance_z = 3 * sigma_z
# 
#         initial_concentration = start_concentration(
#             sigma_x,
#             sigma_y,
#             sigma_z,
#             mass,
#         )
# 
#         for particle_index in range( x.shape[ 1 ] ):
#             particle_x = x[ time_index, particle_index ]
#             particle_y = y[ time_index, particle_index ]
#             particle_z = z[ time_index, particle_index ]
# 
#             # Skip NaN coordinates.
#             if particle_x != particle_x or particle_y != particle_y or \
#                particle_z != particle_z:
#                 continue
# 
#             # Calculate the window of effect of this particle. This defines a
#             # 3D box. Specifically, we're interested in the grid bucket
#             # indices within that box.
#             x_index_min = bucket_index(
#                 particle_x - effect_distance_x,
#                 grid_x_min,
#                 grid_x_max,
#                 bucket_x_delta,
#             )
#             x_index_max = bucket_index(
#                 particle_x + effect_distance_x,
#                 grid_x_min,
#                 grid_x_max,
#                 bucket_x_delta,
#             )
#             y_index_min = bucket_index(
#                 particle_y - effect_distance_y,
#                 grid_y_min,
#                 grid_y_max,
#                 bucket_y_delta,
#             )
#             y_index_max = bucket_index(
#                 particle_y + effect_distance_y,
#                 grid_y_min,
#                 grid_y_max,
#                 bucket_y_delta,
#             )
#             z_index_min = bucket_index(
#                 particle_z - effect_distance_z,
#                 grid_z_min,
#                 grid_z_max,
#                 bucket_z_delta,
#             )
#             z_index_max = bucket_index(
#                 particle_z + effect_distance_z,
#                 grid_z_min,
#                 grid_z_max,
#                 bucket_z_delta,
#             )
# 
#             # For each combination of those indices, calculate a
#             # concentration contribution and add it that grid bucket.
#             bucket_x = bucket_position(
#                 x_index_min, grid_x_min, bucket_x_delta,
#             )
#             bucket_y_min = bucket_position(
#                 y_index_min, grid_y_min, bucket_y_delta,
#             )
#             bucket_z_min = bucket_position(
#                 z_index_min, grid_z_min, bucket_z_delta,
#             )
# 
#             for x_index in range( x_index_min, x_index_max + 1 ):
#                 bucket_y = bucket_y_min
# 
#                 for y_index in range( y_index_min, y_index_max + 1 ):
#                     bucket_z = bucket_z_min
# 
#                     for z_index in range( z_index_min, z_index_max + 1 ):
#                         grid[ z_index, y_index, x_index ] += concentration( initial_concentration,
#                                                                             particle_x,
#                                                                             particle_y,
#                                                                             particle_z,
#                                                                             bucket_x,
#                                                                             bucket_y,
#                                                                             bucket_z,
#                                                                             sigma_x,
#                                                                             sigma_y,
#                                                                             sigma_z,
#                                                                             ) / conversion_factor
#                         bucket_z += bucket_z_delta
# 
#                     bucket_y += bucket_y_delta
# 
#                 bucket_x += bucket_x_delta
# 
#         # Write the grid for this time step out to file.
#         file_grid[ time_index ] = grid
#==============================================================================



#==============================================================================
# Commented-out previous version of function as it's no longer used.  Retained
# for reference.  (Zelenke)
# 
# #@cython.boundscheck( False )
# @cython.cdivision( True )
# def calculate_concentration_timestep(
#     np.ndarray[ np.double_t, ndim = 3 ] grid,
#     np.ndarray[ np.double_t, ndim = 1 ] x not None,
#     np.ndarray[ np.double_t, ndim = 1 ] y not None,
#     np.ndarray[ np.double_t, ndim = 1 ] z not None,
#     np.double_t time, 
#     np.double_t grid_x_min,
#     np.double_t grid_x_max,
#     np.double_t grid_y_min,
#     np.double_t grid_y_max,
#     np.double_t grid_z_min,
#     np.double_t grid_z_max,
#     np.double_t mass,
#     np.ndarray[ np.double_t ] times not None,
#     np.double_t C_x,
#     np.double_t C_y,
#     np.double_t C_z,
#     np.double_t t_0,
#     np.double_t alpha,
#     np.double_t conversion_factor,
# ):
#     """
#     Calculate a concentration grid based on the given particles and save it to
#     file as each time step's grid is calculated.
# 
#     :param grid: concentration grid for a single time step, zeroed
#     :type grid: ZxYxZ array of np.double_t
#     :param file_grid: writable NetCDF4 concentration grid variable object
#     :type file_grid: netCDF4.Variable of TxZxYxZ
#     :param x: x coordinates of particles in meters
#     :type x: array of np.double_t
#     :param y: y coordinates of particles in meters
#     :type y: array of np.double_t
#     :param z: z coordinates of particles in meters (positive up)
#     :type z: array of np.double_t
#     :param grid_x_min: minimum grid x value in meters
#     :type grid_x_min: np.double_t
#     :param grid_x_max: maximum grid x value in meters
#     :type grid_x_max: np.double_t
#     :param grid_y_min: minimum grid y value in meters
#     :type grid_y_min: np.double_t
#     :param grid_y_max: maximum grid y value in meters
#     :type grid_y_max: np.double_t
#     :param grid_z_min: minimum grid z value in meters
#     :type grid_z_min: np.double_t
#     :param grid_z_max: maximum grid z value in meters
#     :type grid_z_max: np.double_t
#     :param mass: mass of each particle in kg
#     :type mass: np.double_t
#     :param times: times in seconds since first time step
#     :type times: array of np.double_t
#     :param sigma_x_0: initial x dimension of the plume in meters
#     :type sigma_x_0: np.double_t
#     :param sigma_y_0: initial y dimension of the plume in meters
#     :type sigma_y_0: np.double_t
#     :param sigma_z_0: initial z dimension of the plume in meters
#     :type sigma_z_0: np.double_t
#     :param alpha: speed of spreading in centimeters per second
#     :type alpha: np.double_t
#     :param conversion_factor: factor to multiply each concentration value by
#     :type conversion_factor: np.double_t
#     """
# #==============================================================================
# #     Subtracted 1 from the denominator in the calculation of each of the three
# #     delta variables below so that grid dimensions calculated by filter.pyx
# #     would match those calculated by grid.py (Zelenke).
# #==============================================================================
#     cdef np.double_t bucket_x_delta = \
#         ( grid_x_max - grid_x_min ) /  ( grid.shape[ 2 ] - 1.0 )
#     cdef np.double_t bucket_y_delta = \
#         ( grid_y_max - grid_y_min ) / ( grid.shape[ 1 ] - 1.0 )
#     cdef np.double_t bucket_z_delta = \
#         ( grid_z_max - grid_z_min ) / ( grid.shape[ 0 ] - 1.0 )
# 
#     # Make one density grid per time step.
#     cdef np.uint32_t particle_index
#     cdef np.double_t sigma_x, sigma_y, sigma_z
#     cdef np.double_t effect_distance_x, effect_distance_y, effect_distance_z
#     cdef np.double_t initial_concentration
#     cdef np.uint32_t x_index_min, x_index_max
#     cdef np.uint32_t y_index_min, y_index_max
#     cdef np.uint32_t z_index_min, z_index_max
#     cdef np.uint32_t x_index, y_index, z_index
#     cdef np.double_t particle_x, particle_y, particle_z
#     cdef np.double_t bucket_x, bucket_y, bucket_z
#     cdef np.double_t bucket_y_min, bucket_z_min
# 
#     # Clear the grid at the start.
#     grid.fill( 0 )
# 
#     # Calculate how far a particle at this time can affect the grid.
#     sigma_x = sigma(
#         C_x, time, t_0, alpha
#     )
#     sigma_y = sigma(
#         C_y, time, t_0, alpha
#     )
#     sigma_z = sigma(
#         C_z, time, t_0, alpha
#     )
# 
#     effect_distance_x = 3 * sigma_x
#     effect_distance_y = 3 * sigma_y
#     effect_distance_z = 3 * sigma_z
# 
#     initial_concentration = start_concentration(
#         sigma_x,
#         sigma_y,
#         sigma_z,
#         mass,
#     )
# 
#     for particle_index in range( x.shape[ 0 ] ):
#         particle_x = x[ particle_index ]
#         particle_y = y[ particle_index ]
#         particle_z = z[ particle_index ]
# 
#         # Skip NaN coordinates.
#         if particle_x != particle_x or particle_y != particle_y or \
#            particle_z != particle_z:
#             continue
# 
#         # Calculate the window of effect of this particle. This defines a
#         # 3D box. Specifically, we're interested in the grid bucket
#         # indices within that box.
#         x_index_min = bucket_index(
#             particle_x - effect_distance_x,
#             grid_x_min,
#             grid_x_max,
#             bucket_x_delta,
#         )
#         x_index_max = bucket_index(
#             particle_x + effect_distance_x,
#             grid_x_min,
#             grid_x_max,
#             bucket_x_delta,
#         )
#         y_index_min = bucket_index(
#             particle_y - effect_distance_y,
#             grid_y_min,
#             grid_y_max,
#             bucket_y_delta,
#         )
#         y_index_max = bucket_index(
#             particle_y + effect_distance_y,
#             grid_y_min,
#             grid_y_max,
#             bucket_y_delta,
#         )
#         z_index_min = bucket_index(
#             particle_z - effect_distance_z,
#             grid_z_min,
#             grid_z_max,
#             bucket_z_delta,
#         )
#         z_index_max = bucket_index(
#             particle_z + effect_distance_z,
#             grid_z_min,
#             grid_z_max,
#             bucket_z_delta,
#         )
# 
#         # For each combination of those indices, calculate a
#         # concentration contribution and add it that grid bucket.
#         bucket_x = bucket_position(
#             x_index_min, grid_x_min, bucket_x_delta,
#         )
#         bucket_y_min = bucket_position(
#             y_index_min, grid_y_min, bucket_y_delta,
#         )
#         bucket_z_min = bucket_position(
#             z_index_min, grid_z_min, bucket_z_delta,
#         )
# 
#         for x_index in range( x_index_min, x_index_max + 1 ):
#             bucket_y = bucket_y_min
# 
#             for y_index in range( y_index_min, y_index_max + 1 ):
#                 bucket_z = bucket_z_min
# 
#                 for z_index in range( z_index_min, z_index_max + 1 ):
#                     grid[ z_index, y_index, x_index ] += concentration( initial_concentration,
#                                                                         particle_x,
#                                                                         particle_y,
#                                                                         particle_z,
#                                                                         bucket_x,
#                                                                         bucket_y,
#                                                                         bucket_z,
#                                                                         sigma_x,
#                                                                         sigma_y,
#                                                                         sigma_z,
#                                                                         ) / conversion_factor
#                     bucket_z += bucket_z_delta
# 
#                 bucket_y += bucket_y_delta
# 
#             bucket_x += bucket_x_delta
# 
#     return None
#==============================================================================



#==============================================================================
# Commented-out previous version of function as it's no longer used.  Retained
# for reference.  (Zelenke)
# #@cython.boundscheck( False )
# @cython.cdivision( True )
# def calculate_concentration_timestep_massvars(
#     np.ndarray[ np.double_t, ndim = 3 ] grid,
#     np.ndarray[ np.double_t, ndim = 1 ] x not None,
#     np.ndarray[ np.double_t, ndim = 1 ] y not None,
#     np.ndarray[ np.double_t, ndim = 1 ] z not None,
#     np.double_t time,
#     np.double_t grid_x_min,
#     np.double_t grid_x_max,
#     np.double_t grid_y_min,
#     np.double_t grid_y_max,
#     np.double_t grid_z_min,
#     np.double_t grid_z_max,
#     np.ndarray[ np.double_t, ndim = 1 ] mass not None,
#     np.ndarray[ np.double_t ] times not None,
#     np.double_t C_x,
#     np.double_t C_y,
#     np.double_t C_z,
#     np.double_t t_0,
#     np.double_t alpha,
#     np.double_t conversion_factor,
# ):
#     """
#     Calculate a concentration grid based on the given particles and save it to
#     file as each time step's grid is calculated. (Zelenke)
# 
#     :param grid: concentration grid for a single time step, zeroed
#     :type grid: ZxYxZ array of np.double_t
#     :param file_grid: writable NetCDF4 concentration grid variable object
#     :type file_grid: netCDF4.Variable of TxZxYxZ
#     :param x: x coordinates of particles in meters
#     :type x: array of np.double_t
#     :param y: y coordinates of particles in meters
#     :type y: array of np.double_t
#     :param z: z coordinates of particles in meters (positive up)
#     :type z: array of np.double_t
#     :param grid_x_min: minimum grid x value in meters
#     :type grid_x_min: np.double_t
#     :param grid_x_max: maximum grid x value in meters
#     :type grid_x_max: np.double_t
#     :param grid_y_min: minimum grid y value in meters
#     :type grid_y_min: np.double_t
#     :param grid_y_max: maximum grid y value in meters
#     :type grid_y_max: np.double_t
#     :param grid_z_min: minimum grid z value in meters
#     :type grid_z_min: np.double_t
#     :param grid_z_max: maximum grid z value in meters
#     :type grid_z_max: np.double_t
#     :param mass: mass of each particle in kg
#     :type mass: array of np.double_t
#     :param times: times in seconds since first time step
#     :type times: array of np.double_t
#     :param sigma_x_0: initial x dimension of the plume in meters
#     :type sigma_x_0: np.double_t
#     :param sigma_y_0: initial y dimension of the plume in meters
#     :type sigma_y_0: np.double_t
#     :param sigma_z_0: initial z dimension of the plume in meters
#     :type sigma_z_0: np.double_t
#     :param alpha: exponent of spreading (unitless)
#     :type alpha: np.double_t
#     :param conversion_factor: factor to multiply each concentration value by
#     :type conversion_factor: np.double_t
#     """
# #==============================================================================
# #     Subtracted 1 from the denominator in the calculation of each of the three
# #     delta variables below so that grid dimensions calculated by filter.pyx
# #     would match those calculated by grid.py (Zelenke).
# #==============================================================================
#     cdef np.double_t bucket_x_delta = \
#         ( grid_x_max - grid_x_min ) /  ( grid.shape[ 2 ] - 1.0 )
#     cdef np.double_t bucket_y_delta = \
#         ( grid_y_max - grid_y_min ) / ( grid.shape[ 1 ] - 1.0 )
#     cdef np.double_t bucket_z_delta = \
#         ( grid_z_max - grid_z_min ) / ( grid.shape[ 0 ] - 1.0 )
# 
#     # Make one density grid per time step.
#     cdef np.uint32_t particle_index
#     cdef np.double_t sigma_x, sigma_y, sigma_z
#     cdef np.double_t effect_distance_x, effect_distance_y, effect_distance_z
#     cdef np.ndarray[ np.double_t, ndim = 1 ] initial_concentration = np.zeros_like(mass)
#     cdef np.uint32_t x_index_min, x_index_max
#     cdef np.uint32_t y_index_min, y_index_max
#     cdef np.uint32_t z_index_min, z_index_max
#     cdef np.uint32_t x_index, y_index, z_index
#     cdef np.double_t particle_x, particle_y, particle_z
#     cdef np.double_t bucket_x, bucket_y, bucket_z
#     cdef np.double_t bucket_y_min, bucket_z_min
# 
#     # Clear the grid at the start.
#     grid.fill( 0 )
# 
#     # Calculate how far a particle at this time can affect the grid.
#     sigma_x = sigma(
#         C_x, time, t_0, alpha
#     )
#     sigma_y = sigma(
#         C_y, time, t_0, alpha
#     )
#     sigma_z = sigma(
#         C_z, time, t_0, alpha
#     )
# 
#     effect_distance_x = 3 * sigma_x
#     effect_distance_y = 3 * sigma_y
#     effect_distance_z = 3 * sigma_z
# 
# #==============================================================================
# #     initial_concentration = start_concentration(
# #         sigma_x,
# #         sigma_y,
# #         sigma_z,
# #         mass,
# #     )
# #==============================================================================
#     #Since mass is now a NumPy ndarray (thus making writing a Cython inline
#     #function difficult) the call to start_concentration has been replaced by
#     #the line below.  This calculation is performed only once per call to this
#     #function, so the difference in performance ought be trivial.  (Zelenke)
#     initial_concentration = mass / ( sigma_x * sigma_y * sigma_z )
# 
#     for particle_index in range( x.shape[ 0 ] ):
#         particle_x = x[ particle_index ]
#         particle_y = y[ particle_index ]
#         particle_z = z[ particle_index ]
#         particle_mass = mass[ particle_index ]
#         particle_initial_concentration = initial_concentration[ particle_index ]
# 
#         # Skip NaN coordinates.
#         if particle_x != particle_x or particle_y != particle_y or \
#            particle_z != particle_z or particle_mass != particle_mass:
#             continue
# 
#         # Calculate the window of effect of this particle. This defines a
#         # 3D box. Specifically, we're interested in the grid bucket
#         # indices within that box.
#         x_index_min = bucket_index(
#             particle_x - effect_distance_x,
#             grid_x_min,
#             grid_x_max,
#             bucket_x_delta,
#         )
#         x_index_max = bucket_index(
#             particle_x + effect_distance_x,
#             grid_x_min,
#             grid_x_max,
#             bucket_x_delta,
#         )
#         y_index_min = bucket_index(
#             particle_y - effect_distance_y,
#             grid_y_min,
#             grid_y_max,
#             bucket_y_delta,
#         )
#         y_index_max = bucket_index(
#             particle_y + effect_distance_y,
#             grid_y_min,
#             grid_y_max,
#             bucket_y_delta,
#         )
#         z_index_min = bucket_index(
#             particle_z - effect_distance_z,
#             grid_z_min,
#             grid_z_max,
#             bucket_z_delta,
#         )
#         z_index_max = bucket_index(
#             particle_z + effect_distance_z,
#             grid_z_min,
#             grid_z_max,
#             bucket_z_delta,
#         )
# 
#         # For each combination of those indices, calculate a
#         # concentration contribution and add it that grid bucket.
#         bucket_x = bucket_position(
#             x_index_min, grid_x_min, bucket_x_delta,
#         )
#         bucket_y_min = bucket_position(
#             y_index_min, grid_y_min, bucket_y_delta,
#         )
#         bucket_z_min = bucket_position(
#             z_index_min, grid_z_min, bucket_z_delta,
#         )
# 
#         for x_index in range( x_index_min, x_index_max + 1 ):
#             bucket_y = bucket_y_min
# 
#             for y_index in range( y_index_min, y_index_max + 1 ):
#                 bucket_z = bucket_z_min
# 
#                 for z_index in range( z_index_min, z_index_max + 1 ):
#                     grid[ z_index, y_index, x_index ] += concentration( particle_initial_concentration,
#                                                                         particle_x,
#                                                                         particle_y,
#                                                                         particle_z,
#                                                                         bucket_x,
#                                                                         bucket_y,
#                                                                         bucket_z,
#                                                                         sigma_x,
#                                                                         sigma_y,
#                                                                         sigma_z,
#                                                                         ) / conversion_factor
#                     bucket_z += bucket_z_delta
# 
#                 bucket_y += bucket_y_delta
# 
#             bucket_x += bucket_x_delta
# 
#     return None
#==============================================================================



#@cython.boundscheck( False )
@cython.cdivision( True )
def calculate_concentration_timestep_agemassvars(
    np.ndarray[ np.double_t, ndim = 3 ] grid,
    np.ndarray[ np.double_t, ndim = 1 ] x not None,
    np.ndarray[ np.double_t, ndim = 1 ] y not None,
    np.ndarray[ np.double_t, ndim = 1 ] z not None,
    np.ndarray[ np.double_t, ndim = 1 ] t not None,
    np.double_t grid_x_min,
    np.double_t grid_x_max,
    np.double_t grid_y_min,
    np.double_t grid_y_max,
    np.double_t grid_z_min,
    np.double_t grid_z_max,
    np.ndarray[ np.double_t, ndim = 1 ] mass not None,
    np.double_t C_x,
    np.double_t C_y,
    np.double_t C_z,
    np.double_t t_0,
    np.double_t alpha,
    np.double_t conversion_factor,
):
    """
    Calculate a concentration grid based on the given particles for saving to
    file as each time step's grid is calculated. (Zelenke)

    :param grid: concentration grid for a single time step, zeroed
    :type grid: ZxYxZ array of np.double_t
    :param x: x coordinates of particles in meters
    :type x: array of np.double_t
    :param y: y coordinates of particles in meters
    :type y: array of np.double_t
    :param z: z coordinates of particles in meters (positive up)
    :type z: array of np.double_t
    :param t: age of particles in seconds (since their release... not necessarily the same as model start)
    :type t: array of np.double_t
    :param grid_x_min: minimum grid x value in meters
    :type grid_x_min: np.double_t
    :param grid_x_max: maximum grid x value in meters
    :type grid_x_max: np.double_t
    :param grid_y_min: minimum grid y value in meters
    :type grid_y_min: np.double_t
    :param grid_y_max: maximum grid y value in meters
    :type grid_y_max: np.double_t
    :param grid_z_min: minimum grid z value in meters
    :type grid_z_min: np.double_t
    :param grid_z_max: maximum grid z value in meters
    :type grid_z_max: np.double_t
    :param mass: mass of each particle in kg
    :type mass: array of np.double_t
    :param sigma_x_0: initial x dimension of the plume in meters
    :type sigma_x_0: np.double_t
    :param sigma_y_0: initial y dimension of the plume in meters
    :type sigma_y_0: np.double_t
    :param sigma_z_0: initial z dimension of the plume in meters
    :type sigma_z_0: np.double_t
    :param t_0: virtual start time for the original point source in seconds
    :type t_0: np.double_t
    :param alpha: exponent of spreading (unitless)
    :type alpha: np.double_t
    :param conversion_factor: factor to multiply each concentration value by
    :type conversion_factor: np.double_t
    """
#==============================================================================
#     Subtracted 1 from the denominator in the calculation of each of the three
#     delta variables below so that grid dimensions calculated by filter.pyx
#     would match those calculated by grid.py (Zelenke).
#==============================================================================
    cdef np.double_t bucket_x_delta = \
        ( grid_x_max - grid_x_min ) / ( grid.shape[ 2 ] - 1.0 )
    cdef np.double_t bucket_y_delta = \
        ( grid_y_max - grid_y_min ) / ( grid.shape[ 1 ] - 1.0 )
    cdef np.double_t bucket_z_delta = \
        ( grid_z_max - grid_z_min ) / ( grid.shape[ 0 ] - 1.0 )

    # Make one density grid per time step.
    cdef np.int32_t particle_index
    cdef np.double_t sigma_x, sigma_y, sigma_z
    cdef np.double_t effect_distance_x, effect_distance_y, effect_distance_z
    cdef np.double_t initial_concentration
    cdef np.uint32_t x_index_min, x_index_max
    cdef np.uint32_t y_index_min, y_index_max
    cdef np.uint32_t z_index_min, z_index_max
    cdef np.uint32_t x_index, y_index, z_index
    cdef np.double_t particle_x, particle_y, particle_z, particle_t, particle_mass
    cdef np.double_t bucket_x, bucket_y, bucket_z
    cdef np.double_t bucket_y_min, bucket_z_min
    
    

    # Clear the grid at the start.
    grid.fill( 0 )
    
    
    for particle_index in range( x.shape[ 0 ] ):
        particle_x = x[ particle_index ]
        particle_y = y[ particle_index ]
        particle_z = z[ particle_index ]
        particle_t = t[ particle_index ]
        particle_mass = mass[ particle_index ]
        
        
        # Skip NaN coordinates.
        if particle_x != particle_x or particle_y != particle_y or \
           particle_z != particle_z or particle_mass != particle_mass:
            continue
        
        # Calculate how far this particle, given its age, can affect the grid.
        sigma_x = sigma( C_x, particle_t, t_0, alpha )
        sigma_y = sigma( C_y, particle_t, t_0, alpha )
        sigma_z = sigma( C_z, particle_t, t_0, alpha )
        
        effect_distance_x = 3 * sigma_x
        effect_distance_y = 3 * sigma_y
        effect_distance_z = 3 * sigma_z
        
        initial_concentration = start_concentration(
            sigma_x,
            sigma_y,
            sigma_z,
            particle_mass,
        )

        # Calculate the window of effect of this particle. This defines a 3-D
        # box. Specifically, we're interested in the grid bucket indices within
        # that box.
        x_index_min = bucket_index(
            particle_x - effect_distance_x,
            grid_x_min,
            grid_x_max,
            bucket_x_delta,
        )
        x_index_max = bucket_index(
            particle_x + effect_distance_x,
            grid_x_min,
            grid_x_max,
            bucket_x_delta,
        )
        y_index_min = bucket_index(
            particle_y - effect_distance_y,
            grid_y_min,
            grid_y_max,
            bucket_y_delta,
        )
        y_index_max = bucket_index(
            particle_y + effect_distance_y,
            grid_y_min,
            grid_y_max,
            bucket_y_delta,
        )
        z_index_min = bucket_index(
            particle_z - effect_distance_z,
            grid_z_min,
            grid_z_max,
            bucket_z_delta,
        )
        z_index_max = bucket_index(
            particle_z + effect_distance_z,
            grid_z_min,
            grid_z_max,
            bucket_z_delta,
        )

        # For each combination of those indices, calculate a
        # concentration contribution and add it that grid bucket.
        bucket_x = bucket_position(
            x_index_min, grid_x_min, bucket_x_delta,
        )
        bucket_y_min = bucket_position(
            y_index_min, grid_y_min, bucket_y_delta,
        )
        bucket_z_min = bucket_position(
            z_index_min, grid_z_min, bucket_z_delta,
        )
        
        
        for x_index in range( x_index_min, x_index_max + 1 ):
            bucket_y = bucket_y_min

            for y_index in range( y_index_min, y_index_max + 1 ):
                bucket_z = bucket_z_min

                for z_index in range( z_index_min, z_index_max + 1 ):
                    
                    grid[ z_index, y_index, x_index ] += concentration( initial_concentration,
                                                                        particle_x,
                                                                        particle_y,
                                                                        particle_z,
                                                                        bucket_x,
                                                                        bucket_y,
                                                                        bucket_z,
                                                                        sigma_x,
                                                                        sigma_y,
                                                                        sigma_z,
                                                                        ) / conversion_factor
                    bucket_z += bucket_z_delta

                bucket_y += bucket_y_delta

            bucket_x += bucket_x_delta

    return None



@cython.boundscheck( False )
@cython.cdivision( True )
cdef inline np.double_t sigma(
    np.double_t C,
    np.double_t t,
    np.double_t t_0,
    np.double_t alpha
):
    """
    Given some constants and a particular time step in seconds, calculate the
    value of sigma at that time for a particular dimension.
    """
    return ( C * (t + t_0) )**alpha


@cython.boundscheck( False )
@cython.cdivision( True )
cdef inline np.double_t start_concentration(
    np.double_t sigma_x,
    np.double_t sigma_y,
    np.double_t sigma_z,
    np.double_t mass,
):
    """
    Given sigma values for a particular time step, and mass for a particle at
    that time step, calculate the initial concentration.
    """
    return mass / ( sigma_x * sigma_y * sigma_z )



@cython.boundscheck( False )
@cython.cdivision( True )
cdef inline np.double_t concentration(
    np.double_t initial_concentration,
    np.double_t particle_x,
    np.double_t particle_y,
    np.double_t particle_z,
    np.double_t bucket_x,
    np.double_t bucket_y,
    np.double_t bucket_z,
    np.double_t sigma_x,
    np.double_t sigma_y,
    np.double_t sigma_z,
):
    """
    Given several values, calculate the concentration contribution of a given
    particle to a given grid bucket.
    """
    # Denominator literal is (2*pi) ** (3/2)
    return \
        ( initial_concentration / 15.749609945722419 ) * \
        exp(
            -0.5 * ( ( bucket_x - particle_x ) / sigma_x ) ** 2 +
            -0.5 * ( ( bucket_y - particle_y ) / sigma_y ) ** 2 +
            -0.5 * ( ( bucket_z - particle_z ) / sigma_z ) ** 2,
        )



@cython.boundscheck( False )
@cython.cdivision( True )
cdef inline np.uint32_t bucket_index(
    np.double_t position,
    np.double_t bucket_min,
    np.double_t bucket_max,
    np.double_t bucket_delta,
):
    """
    Given a position and data about the grid, return a corresponding bucket
    index.
    """
    # Even if the position is beyond the grid, don't return an index out of
    # bounds.
    if position < bucket_min:
        return 0
    if position >= bucket_max:
        return <np.uint32_t>( ( ( bucket_max - bucket_min ) / bucket_delta ) - 1 )

    return <np.uint32_t>( ( position - bucket_min ) / bucket_delta )




@cython.boundscheck( False )
@cython.cdivision( True )
cdef inline np.double_t bucket_position(
    np.uint32_t bucket_index,
    np.double_t bucket_min,
    np.double_t bucket_delta,
):
    """
    Given a bucket index and data about the grid, return a corresponding
    bucket position.
    """
    return bucket_min + ( bucket_index * bucket_delta )
