#!/usr/bin/env python

"""
makes a GNOME map from the GSSHHS or NOS shorelines
"""
from __future__ import print_function
import os
import json
import numpy as np
import shapely.geometry as sgeo
import fiona


def format_for_bna(polys, levs, dateline=False):
    data = []
    ii = 1
    for lev, poly in zip(levs, polys):
        x = np.array(poly[0])
        y = np.array(poly[1])
        if dateline:
            x = (x < 0).choose(x, x + 360)
        if lev == 0:
            data.append('"Map Bounds", "2", ' + str(len(x)))
        else:
            data.append('"' + str(ii) + '","' + str(int(lev)) + '",' + str(len(x)))
        for lon, lat in zip(x, y):
            data.append("{0:.6f}".format(lon) + ', ' + "{0:.6f}".format(lat))
        ii = ii + 1
    data = '\n'.join(data)
    return data


def make_geojson_feature(poly):
        rings = []
        ring = poly.exterior
        if ring.is_ccw:
            coords = list(ring.coords)
        else:
            coords = list(ring.coords)[::-1]
        rings.append(coords)
        for ring in poly.interiors:
            if not ring.is_ccw:
                coords = list(ring.coords)
            else:
                coords = list(ring.coords)[::-1]
            rings.append(coords)

        feat = {"type": "Feature",
                "geometry": {"type": "Polygon",
                             "coordinates": rings
                             },
                "properties": {}
                }
        return feat


def format_for_geojson(polys, map_bounds=None, spillable_area=None, dateline=False):
    """
    Create geojson of gnome map as a string from a multipolygon object.

    :param polys: A shapely MultiPolygon object with polygons around land
    """
    features = []
    gj = {"type": "FeatureCollection",
          "features": features}
    # add the map_bounds
    if map_bounds is not None:
        feat = make_geojson_feature(map_bounds)
        feat['properties']["poly_type"] = "Map Bounds"
        features.append(feat)
    if spillable_area is not None:
        feat = make_geojson_feature(spillable_area)
        feat['properties']["poly_type"] = "Spillable Area"
        features.append(feat)
    for poly in polys:
        feat = make_geojson_feature(poly)
        feat['properties']["poly_type"] = "land"
        features.append(feat)
    return json.dumps(gj, indent=2)


def load_and_clip_gshhs(data_dir,
                        resolution,
                        wl, el, nl, sl,
                        dateline=False):
    """
    load up the shoreline from gshhs.

    the water levels are clipped from the land levels to result in land-only polygons
    """
    bound_box = (wl, sl, el, nl)
    clip_box = sgeo.box(*bound_box)
    map_bounds = clip_box
    land_polys = sgeo.MultiPolygon()
    for level in [1, 2, 3, 4]:
        filename = os.path.join(data_dir,
                                '{0}/GSHHS_{0}_L{1}.shp'.format(resolution, level))
        print("reading:", filename)
        if not os.path.exists(filename):
            continue
        with fiona.open(filename) as source:
            for polygon in source.filter(bbox=bound_box):
                geom = polygon['geometry']
                if geom.get('type') != 'Polygon':
                    raise ValueError("There is something other than a polygon in the data")
                coords = geom['coordinates']
                if len(coords) != 1:
                    raise ValueError("There is a polygon with multiple rings")
                poly = sgeo.Polygon(coords[0])
                if dateline:  # box depends on the which side of dateline we're worried about
                    poly_bnds = poly.bounds
                    if poly_bnds[0] > -10 and poly_bnds[2] <= 180:
                        clip_box = sgeo.box(wl, sl, 180, nl)
                    else:
                        clip_box = sgeo.box(-180, sl, el, nl)
                clipped_polys = poly.intersection(clip_box)
                clipped_polys = as_multi_polygon(clipped_polys)

                if level % 2:  # odd level, it's land -- add them to the land polys
                    # fixme: there has got to be a better way to merge multipolygons!
                    land_polys = sgeo.MultiPolygon(list(land_polys.geoms) +
                                                   list(clipped_polys.geoms))
                else:  # water -- clip them out of land
                    # clip water out of existing land
                    land_polys = as_multi_polygon(land_polys.difference(clipped_polys))
    return land_polys, map_bounds


def as_multi_polygon(poly):
    if not isinstance(poly, sgeo.MultiPolygon):
        # got a single polygon, make a multiple out of it
        return sgeo.MultiPolygon([poly])
    else:
        return poly





        # if type(clipped_polys) is sgeo.multipolygon.MultiPolygon:
        #     # There are multiple polygons
        #     for clipped_poly in clipped_polys:
        #         try:
        #             bna_polys.append(clipped_poly.exterior.xy)
        #         except AttributeError:
        #             pass
        # elif type(clipped_polys) is sgeo.polygon.Polygon:
        #     try:
        #         bna_polys.append(clipped_polys.exterior.xy)
        #     except AttributeError:
        #         pass



def load_polys_and_clip(shpfile, wl, el, nl, sl, dateline=False, lev=1):
    """
    Loop through polygons in shapefile, if they overlap with user
    specified rectangle, load the polygon and clip
    Note that GSHHS polygons are -180-->180, shorelines that cross the 180 line are
    clipped into west and east polygons, so there is some special casing to handle
    clip rectangles that cross the dateline
    Change lev to 2 if loading lakes
    These *may not* be able to be used to clip trajectories from GAn in Arc?
    """

    bound_box = (wl, sl, el, nl)
    clip_box = sgeo.box(*bound_box)

    if dateline:
        bound_box = (-180, sl, 180, nl)

    bna_polys = []

    with fiona.open(shpfile) as source:
        for s in source.filter(bbox=bound_box):
            geo = s['geometry']['coordinates']
            for c in geo:
                poly = sgeo.Polygon(c)
                if dateline:
                    poly_bnds = poly.bounds
                    if poly_bnds[0] > -10 and poly_bnds[2] <= 180:
                        clip_box = sgeo.box(wl, sl, 180, nl)
                    else:
                        clip_box = sgeo.box(-180, sl, el, nl)

                clipped_polys = poly.intersection(clip_box)
                if type(clipped_polys) is sgeo.multipolygon.MultiPolygon:
                    for clipped_poly in clipped_polys:
                        try:
                            bna_polys.append(clipped_poly.exterior.xy)
                        except AttributeError:
                            pass
                elif type(clipped_polys) is sgeo.polygon.Polygon:
                    try:
                        bna_polys.append(clipped_polys.exterior.xy)
                    except AttributeError:
                        pass

    if len(bna_polys) < 1 and lev == 1:  # all water
        bna_polys.append(clip_box.exterior.xy)
        levels = [2]
    else:  # define levels and add lakes if specified
        levels = [lev for x in range(len(bna_polys))]

    if lev == 1:
        bna_polys.insert(0, ([wl, wl, el, el], [sl, nl, nl, sl]))
        levels.insert(0, 0)

    return bna_polys, levels
