==================================================================
Generating a gridded concentration field from Lagrangian particles
==================================================================

Many oil and contaminant transport models rely on a particle tracking approach. In this approach a spill is broken down into many individual particles, each of which represents a parcel of oil or other contaminant. Each parcel can have various properties, such as  mass, density, rise velocity, etc. The movement of each particle is then tracked as it is advected by the currents, winds, etc. This approach has the benefit of simplicity, can capture broad spatial scales, and is not subject to numerical diffusion. It is used in such models as CDOG, GNOME, CBOLT, etc.

However, once the particle tracking model has been run, the end product desired is often not just the location of the contaminant, but the concentration. Thus, there need to be ways to convert a set of particle locations into a gridded field of concentrations. The simplest way to do this is to simply divide your area of interest into a grid, and count the number of particles in each bin multiplied by their mass. This total mass can be divided into the volume of the bin to get the concentration. While simple and easy to compute, the results of this approach are quite sensitive to the bin size chosen, and in regions where the particles are few and far between, can result in patchiness that may not reflect the real distribution of the contaminant.

An alternative approach is to consider each particle not as a point, but as the center of a distribution. The integral of the distribution is the mass of the particle. Thus if you add of the contribution from each individual particle's distribution, you can get the concentration at that point. This is often know as kernel density estimation. This results in a way to calculate the concentration on any desired grid without sensitivity to the grid spacing. The results can, however be sensitive to the width of the distribution used. By having the distribution expand as the particle ages, smooth results can be obtained for a wide range of particle spacings. Physically speaking, the particles can be thought of as diffusing at a spatial scale smaller than the scale of the diffusion in the particle tracking model.

This leaves the question of what kernel to use, and how to compute its width. Given the diffusion as a core natural process, a Gaussian kernel is natural in this application. The code allows the standard deviation of the kernel to grow with time according to a user set power law. This allows the user to have the kernel remain constant with time, (exponent 1) or grow at any appropriate rate. As Fickian diffusion predicts growth of the standard deviation with the square root of time, we chose an exponent of 1/2, and a growth rate in each direction to match the diffusion parameter used in the particle tracking model. The result is that the size of each particle grows at the about the same rate that the particles move away from each-other resulting in a smooth concentration field over time.

Requirements
------------

The following python packages and binary packages are required before building grid_plume packages

* numpy>=1.25
* Cython>=3.0
* netCDF4
* pyproj  

Install
-------

Dependencies
............

All dependencies are available on conda-forge:

`conda install -c conda-forge --file conda_requirements.txt`

Build and Install
.................

::

    pip install -e .

    or

    pip install .

Usage
-----

The gridplume package installs the following scripts::

    gp_grid_plume
    gp_merge_grids
    gp_swept_area


