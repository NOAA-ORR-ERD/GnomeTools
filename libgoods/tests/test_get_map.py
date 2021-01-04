#!/usr/bin/env python

"""
testing code for get_map

not proper unit tests, but something...

"""
from __future__ import print_function
import json
import fiona

from shapely.geometry import Polygon, MultiPolygon

from libgoods.get_map import load_and_clip_gshhs, format_for_geojson

GSHHS_dir = "/Users/chris.barker/Temp/GSHHG/gshhg-shp-2.3.7/GSHHS_shp"
# A few test cases:

# Piece of great lakes with bit of shoreline and island:
# wl, el, nl, sl,
gl1 = (-85.273, -84.894, 45.724, 45.62)

# Piece of great lakes with multiple shorelines and islands:
# wl, el, nl, sl,
gl2 = (-86.0, -84.814, 46.132, 45.527)


# simple straight piece of long island
# wl, el, nl, sl,
li = (-72.84, -72.616, 40.808, 40.695)

# wl, el, nl, sl,
wa_coast = (-124.478, -124.294, 47.648, 47.567)

# island doesn't appear to be in the gshhs shoreline
wa_coast_with_island = (-124.548, -124.295, 47.70, 47.629)

# Puget Sound
puget_sound = (-124.0, -122.0, 49.0, 48.0)


# Resolution options:
# [f]ull
# [h]igh
# [i]termediate
# [l]ow
# [c]ourse


def test_west_lake_erie():
    # wl, el, nl, sl,
    bounds = (-83.859, -81.343, 42.255, 41.189)
    land, map_bounds = load_and_clip_gshhs(GSHHS_dir, 'f', *bounds)
    gj = format_for_geojson(land,
                            # map_bounds=map_bounds
                            )
    open("lake_erie.geojson", 'w').write(gj)


def test_lake_st_clair():
    # wl, el, nl, sl,
    bounds = (-83.166, -82.255, 42.849, 42.226)
    land, map_bounds = load_and_clip_gshhs(GSHHS_dir, 'f', *bounds)
    gj = format_for_geojson(land,
                            # map_bounds=map_bounds
                            )
    open("lake_st_clair.geojson", 'w').write(gj)


def test_great_lakes():
    land, map_bounds = load_and_clip_gshhs(GSHHS_dir, 'f', *gl1)
    gj = format_for_geojson(land,
                            map_bounds=map_bounds
                            )
    open("great_lakes_1.geojson", 'w').write(gj)

def test_great_lakes2():
    land, map_bounds = load_and_clip_gshhs(GSHHS_dir, 'f', *gl2)
    gj = format_for_geojson(land,
                            # map_bounds=map_bounds
                            )
    open("great_lakes_2.geojson", 'w').write(gj)


def test_wa_coast():
    """ simple straight shoreline """

    land, map_bounds = load_and_clip_gshhs(GSHHS_dir, 'i', *wa_coast)
    gj = format_for_geojson(land, map_bounds)

    open("wa_shore.geojson", 'w').write(gj)


def test_puget_sound():
    """ lots of islands """

    land, map_bounds = load_and_clip_gshhs(GSHHS_dir, 'f', *puget_sound)
    gj = format_for_geojson(land, map_bounds)

    open("puget_sound.geojson", 'w').write(gj)


def test_south_florida():
    """ Lake Okeechobee is in there! """
    # wl, el, nl, sl,
    south_florida = (-83.0, -79.5, 27.5, 25.0)
    land, map_bounds = load_and_clip_gshhs(GSHHS_dir, 'h', *south_florida)
    gj = format_for_geojson(land, map_bounds)

    open("south_florida.geojson", 'w').write(gj)


# def test_wa_coast_with_island():
#     """ simple straight shoreline with islands offshore
#           disabled, as the island doesn't appear tobe in the data
#     """

#     land, map_bounds = load_and_clip_gshhs(GSHHS_dir, 'f', *wa_coast_with_island)
#     gj = format_for_geojson(land, map_bounds)

#     open("wa_shore_island.geojson", 'w').write(gj)


def test_long_island():
    """ a barrier island -- 3 polys """

    land, map_bounds = load_and_clip_gshhs(GSHHS_dir, 'i', *li)
    gj = format_for_geojson(land)

    open("long_island_shore.geojson", 'w').write(gj)


def test_format_for_geojson():
    poly1 = Polygon([(-78, 41), (-75.5, 41), (-76, 44), (-78, 45)])
    poly2 = Polygon([(-76, 46), (-75, 46), (-75.5, 48)])
    mpoly = MultiPolygon([poly1, poly2])
    bounds = Polygon([(-80, 40), (-80, 50), (-75, 50), (-76, 40)])
    spillable = Polygon([(-78, 40), (-76, 49), (-76, 49), (-76, 42)])
    gj = format_for_geojson(mpoly,
                            map_bounds=bounds,
                            spillable_area=spillable)

    open("simple_test.geojson", 'w').write(gj)

    gjdict = json.loads(gj)

    print(gj)

    assert gjdict['type'] == "FeatureCollection"
    assert len(gjdict['features']) == 4

