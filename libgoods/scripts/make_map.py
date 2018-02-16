#!/usr/bin/env python

"""
script to make a GNOME map (in BNA format) from GSHHS or NOS shoreline data
"""

# from libgoods import get_map
import os
import sys


# data_dir = 'C:\\Users\\amy.macfadyen\\PyProjects\\GOODS\\trunk\\static\\coastline_data\\gshhs_shp'
# data_dir = 'C:\\Users\\amy.macfadyen\\PyProjects\\GOODS\\trunk\\static\\coastline_data\\nos_shoreline'


def main(coast_data='GSHHS',  # GSHHS or NOS
         data_dir='',
         bounds=((-130, 40), (-120, 50)),
         resolution='h',  # can be f/h/i/l/c (from highest res to lowest)
         ):
    lon_min = bounds[0][0]
    lat_min = bounds[0][1]
    lon_max = bounds[1][0]
    lat_max = bounds[1][1]


    if coast_data == 'GSHHS':
        bna_polys = []
        for level in [1, 2, 3, 4]:
            if level == 4 and resolution == 'c':
                continue  # course res has no level 4 data
            shpfile = os.path.join(data_dir,
                                   '{0}/GSHHS_{0}_L{1}.shp'.format(resolution, level))
            print("loading:", shpfile)
            assert os.path.isfile(shpfile)
            # bna_polys2, levs = get_map.load_polys_and_clip(shpfile,
            #                                               lon_min,
            #                                               lon_max,
            #                                               lat_max,
            #                                               lat_min,
            #                                               lev=level)
            bna_polys.extend(bna_polys2)
    #         levs.extend(levs2)
    # elif coast_data == 'NOS':
    #     shpfile = os.path.join(data_dir,'nosshore.shp')
    #     bna_polys, levs = get_map.load_polys_and_clip(shpfile,lon_min,lon_max,nl,lat_min)

    # data = get_map.format_for_bna(bna_polys, levs)

    # f = open(coast_data + '_coast.bna','w')
    # f.write(data)
    # f.close()


HELP = """
make_map.py

Makes a GNOME BNA from GSHHS or NOS Shoreline data

make_map.py data_dir bounds [res]

data_dir: is the path to the directory the data are in

bounds: are the bounds you want the map for, in the form:

        "min_lon,min_lat,max_lon,max_lat"

res: the desired resolution (GSHHS data). One of:
     'f', 'h', 'i', 'l', 'c'
     fine, high, intermediate, low, course
"""


if __name__ == "__main__":

    if len(sys.argv) == 1:
        print HELP
        sys.exit(1)
    try:
        data_dir = sys.argv[1]
    except IndexError:
        print HELP
        sys.exit(2)
    try:
        bounds = sys.argv[2]
        min_lon, min_lat, max_lon, max_lat = [float(s) for s in bounds.split(',')]
    except IndexError:
        print HELP
        sys.exit(3)
    except ValueError:
        print "bounds must be a single string with 4 comma separated values"
        sys.exit(3)


    main(coast_data='GSHHS',
         data_dir="/Users/chris.barker/Temp/GSHHG/gshhg-shp-2.3.7/GSHHS_shp"
         )

    # # make sure all files are there
    # for res in ['f', 'h', 'i', 'l', 'c']:
    #         main(coast_data='GSHHS',
    #              data_dir="/Users/chris.barker/Temp/GSHHG/gshhg-shp-2.3.7/GSHHS_shp",
    #              resolution=res
    #              )
