#!/usr/bin/env python

import numpy as np
from scipy.stats import norm
from scipy import spatial

def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]

def plumeKDE_ellipse(p_posit, p_mass, out_posit, n_neighbors=20, pdf_scale=1.):
    """
    Ellipse version. Using distance to n_neighbors to define the length scales and rotation
    angle of an elliptical kernel.
    """
    # defaults
    # n_neighbors = 20.    # number of neighbors to use to set KDE scale
    mx_dist_scale = 3.  # scale for max distance to look for particles within
    # pdf_scale = 1.         # scale for pdf...scales mean distance to neighbor points

    # test for 2D/3D, assume 2D for now
    p_xy = p_posit

    # call to KDtree
    tree = spatial.KDTree(p_xy)

    # find N nearest neighbors to each particle
    tr_dist, tr_indx = tree.query(p_xy, n_neighbors)

    # mean positions of nearest neighbors (using results from tree call)
    xm = (np.mean(p_xy[tr_indx,0],axis=1))
    ym = (np.mean(p_xy[tr_indx,1],axis=1))

    # arrays of x/y positions of neighbors [np,n_neighbors]
    xtemp = p_xy[tr_indx,0]
    ytemp = p_xy[tr_indx,1]
    #
    for kk in range(len(xm)):
        cov[kk] = np.cov(xtemp[kk,:],ytemp[kk,:])
        vals, vecs = eigsorted(cov)
    
        theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
        w, h = 2 * nstd * np.sqrt(vals)


    # mean dist to these neighbors for pdf scaling
    p_meandist = np.mean(tr_dist, axis=1)
    p_max = np.max(p_meandist)

    # construct out_points tree
    # assume 2D for now
    out_xy = out_posit
    tree_out = spatial.KDTree(out_xy)

    # find output points within a certain distance of each particle position
    p_indx = tree.query_ball_tree(tree_out, mx_dist_scale * p_max)

    # ## loop through particles and sum up KDE influences at each output point
    # (using particle masses to scale KDE)

    # initialize output values
    out_data = np.zeros(np.size(out_posit, axis=0))

    for kk in range(len(p_indx)):
        if p_indx[kk]:
            # find distances to indexed points
            dtemp = spatial.distance_matrix(p_xy[kk:kk + 1, :], out_xy[p_indx[kk], :])

            # scale
            #  use nearest-neighbors mean distances to define lengths of KDE
            #   (scaling by 1/sqrt(2pi) since I'm using a 1D pdf function)
            ntemp = norm.pdf(dtemp[0], loc=0, scale=pdf_scale * p_meandist[kk]) / np.sqrt(2 * np.pi);

            # sum up KDEs at output locations
            out_data[p_indx[kk]] = out_data[p_indx[kk]] + p_mass[kk] * ntemp

    return out_data, p_meandist



def plumeKDE(p_posit, p_mass, out_posit, n_neighbors=20, pdf_scale=1.):
    """
    Original version. Using distance to n_neighbors to define the length scales of
    the kernel. Using a normal distribution for the kernel shape.
    """
    # defaults
    # n_neighbors = 20.    # number of neighbors to use to set KDE scale
    mx_dist_scale = 1000.  # scale for max distance to look for particles within
    # pdf_scale = 1.         # scale for pdf...scales mean distance to neighbor points

    # test for 2D/3D, assume 2D for now
    p_xy = p_posit

    # call to KDtree
    tree = spatial.KDTree(p_xy)

    # find N nearest neighbors to each particle
    tr_dist, tr_indx = tree.query(p_xy, n_neighbors)

    # mean dist to these neighbors for pdf scaling
    p_meandist = np.mean(tr_dist, axis=1)
    p_max = np.max(p_meandist)
    print p_max

    # construct out_points tree
    # assume 2D for now
    out_xy = out_posit
    tree_out = spatial.KDTree(out_xy)

    # find output points within a certain distance of each particle position
    p_indx = tree.query_ball_tree(tree_out, mx_dist_scale * p_max)

    # ## loop through particles and sum up KDE influences at each output point
    # (using particle masses to scale KDE)

    # initialize output values
    out_data = np.zeros(np.size(out_posit, axis=0))

    for kk in range(len(p_indx)):
        if p_indx[kk]:
            # find distances to indexed points
            dtemp = spatial.distance_matrix(p_xy[kk:kk + 1, :], out_xy[p_indx[kk], :])

            # scale
            #  use nearest-neighbors mean distances to define lengths of KDE
            #   (scaling by 1/sqrt(2pi) since I'm using a 1D pdf function)
            ntemp = norm.pdf(dtemp[0], loc=0, scale=pdf_scale * p_meandist[kk]) / np.sqrt(2 * np.pi);

            # sum up KDEs at output locations
            out_data[p_indx[kk]] = out_data[p_indx[kk]] + p_mass[kk] * ntemp

    return out_data, p_meandist


# ---------------------------------------------------------------------------
def plumeKDE_tophat(p_posit, p_mass, out_posit, n_neighbors=5):
    """
    Using distance to n_neighbors to define the length scales of
    the kernel. Using a tophat distribution for the kernel shape.

    """
    # defaults
    # n_neighbors = 5    # number of neighbors to use to set KDE scale

    # test for 2D/3D, assume 2D for now
    p_xy = p_posit

    # call to KDtree
    tree = spatial.KDTree(p_xy)

    # find distances to N nearest neighbors for each particle
    tr_dist, tr_indx = tree.query(p_xy, n_neighbors)

    # max dist to these neighbors for pdf scaling
    p_meandist = np.mean(tr_dist, axis=1)
    p_radii = np.max(tr_dist, axis=1)

    # construct out_points tree
    # assume 2D for now
    out_xy = out_posit
    tree_out = spatial.KDTree(out_xy)

    # ## loop through particles and sum up KDE influences at each output point
    # (using particle masses to scale KDE)

    # initialize output values
    out_data = np.zeros(np.size(out_posit, axis=0))
    out_count = np.zeros(np.size(out_posit, axis=0))

    for kk in range(len(p_radii)):
        # find out_points w/in radius defined by this particle's N-neighbor-radius
        i_tmp = tree_out.query_ball_point(p_xy[kk, :], p_radii[kk])

        out_data[i_tmp] = out_data[i_tmp] + p_mass[kk] / (2. * np.pi * p_radii[kk] ** 2)
        out_count[i_tmp] = out_count[i_tmp] + 1.

    return out_data, p_meandist, out_count


def plumeKDE_adios(p_posit, p_mass, p_density, out_posit, p_age):
    """
    Using Bill's description of oil spreading to define the size of oil patches.

    """
    # input data, change inputs
    n_LE = len(p_posit)  # number of LEs released
    V0 = 0.01;  # initial volume of each LE (m3)
    dbuoy = 0.2  # (rho_h20 - rho_oil)/rho_h20     # relative buoyancy of oil
    # LE_age = 30*60                          # age of LE (sec)

    # defaults (and placeholdrs for input LE data)
    k1, k2 = [1.53, 1.21]  # Fay spreading constants
    nu_h2o = 0.000001;  # water viscosity (m2/sec)
    g = 9.81  # gravity (m2/s)

    ## This approach splits the spreading of the oil into two phases: the initial
    ##    spreading is defined by gravity-inertial forces. This stage scales to
    ##    be relatively quick (compared to typical model time steps), and so we use
    ##    it to give an initial thickness. The next phase is governed by gravity-viscous
    ##    speading and diffusion, and is time-dependent.

    ## NOTE: this assumes all LEs have same initial volume...sensible?


    ## initial area covered by oil release, as defined by Fay intertial spreading
    A0 = np.pi * (k2 ** 4 / k1 ** 2) * (((n_LE * V0) ** 5 * g * dbuoy) / (nu_h2o ** 2)) ** (1. / 6.)

    # # average initial thickness of LE's released during time-step
    # Th_0 = n_LE*V0/A0
    # check for minimum thickness (0.01mm?. based on initial volume of LE's? or
    #    actual current volume)

    # # Fay gravity-viscous and non-diffusion
    dFay = k2 ** 2. / 16. * (g * dbuoy * V0 ** 2 / np.sqrt(nu_h2o * p_age))
    dEddy = 0.033 * p_age ** (4 / 25)

    p_Area = A0 + (dFay + dEddy) * p_age
    p_radii = np.sqrt(p_Area / np.pi)

    # test for 2D/3D, assume 2D for now
    p_xy = p_posit
    # call to KDtree for LE points
    tree = spatial.KDTree(p_xy)

    # construct out_points tree
    # assume 2D for now
    out_xy = out_posit
    tree_out = spatial.KDTree(out_xy)

    # find output points within a certain distance of each particle position
    p_indx = tree.query_ball_tree(tree_out, p_radii)

    # ## loop through particles and sum up KDE influences at each output point
    #      (using particle masses to scale KDE)

    # initialize output values
    out_data = np.zeros(np.size(out_posit, axis=0))

    for kk in range(len(p_posit)):
        if p_indx[kk]:
            # add up at out_points w/in radius defined by this particle's Fay-radius
            out_data[p_indx[kk]] = out_data[p_indx[kk]] + (p_mass[kk] / p_density[kk]) / (2. * np.pi * p_radii ** 2)

    return out_data, p_radii
