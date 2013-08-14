"""
TAP Setup.py

Master set up script for a TAP run

All the data required to set up and build TAP cubes + site.txt file should be in here

"""

import os, datetime

RootDir =  os.path.split(__file__)[0]
print "Loading TAP Setup data from:", RootDir

## Cube Builder Data
ReceptorType = "Grid" # should be either "Grid" or "Polygons" (only grid is supported at the moment)
CubeType = "Volume" # should be either "Volume" or "Cumulative"

## CubeDataType options are: 'float32', 'uint16', or 'uint8'
##   float32 gives better precision for lots of LEs
##   uint8 saves disk space -- and is OK for about 1000 LEs
##   uint16 is a mid-point -- probably good to 10,000 LEs or so
CubeDataType = 'float32' 

## Batch GNOME Data:
## You can use multiple machines in parallel to do many runs
NumTrajMachines = 12

# Files with time series records in them used by GNOME
# These are used to compute the possible time files. The format is:
# It is a list of one or more time files. each file is desribed with a tuple:
#  (file name, allowed_gap_length, type)
#    file_name is a string
#    allowed_gap_length is in hours. It indicates how long a gap in the time
#         series records you will allow GNOME to interpolate over.
#    type is a string describing the type of the time series file. Options
#         are: "Wind", "Hyd" for Wind or hydrology type files
# if set to None, model start an end times will be used
#TimeSeries = [("WindData.OSM", datetime.timedelta(hours = 6), "Wind" ),]
TimeSeries = None

# 2003 data
DataStartEnd = (datetime.datetime(1992, 11, 1),
                datetime.datetime(1999, 11, 20)
                 )
# 2009 data:
##DataStartEnd = (datetime.datetime(2000, 1, 1),
##               datetime.datetime(2008, 5, 06)
##               )

DataGaps = ( )


# specification for how you want seasons to be defined:
#  a list of lists:
#  [name, (months) ]
#    name is a string for the season name  
#    months is a tuple of integers indicating which months are in that season

Seasons = [["90 days", range(1,13) ],
               ]

StartTimeFiles = [(os.path.join(RootDir, s[0]+'Starts.txt'), s[0]) for s in Seasons]

# number of start times you want in each season:
#NumStarts = 500
NumStarts = 200

# Length of release in hours
ReleaseLength = 24 * 90 #24 hrs * XX days

ModelTimeStep = 3 * 60 # minutes
OutputTimeStep = 3 * 60# minutes

# name of the GNOME SAV file you want to use
# note: GNOME locks it (for a very brief time when loading) 
# which means you need a separate copy for each
# instance of GNOME you want to run (OR just don't start multiple GNOMES too quickly)
SaveFileName = os.path.join('Q:\EastCoast_TAP\SaveFiles\CompoundMovers_1999.sav')

# location of desired GNOME to use -- if this path is wrong or this is not specified
# then the auto starting of all the trajectories won't work
# (and at this point it will not exit gracefully)
GNOME = 'Q:\EastCoast_TAP\GNOME\Gnome_NetCDF.exe'                            

# number of Lagrangian elements you want in the GNOME run
NumLEs = 36000

LE_Windage = (0.0001, 0.03) # windage as a fraction (0.0 not allowed?)
                             # set to None for the default
                            
# we only have "MediumCrude"  in the data for now (see OilWeathering.py)
# OilWeatheringType = 'MediumCrude'
OilWeatheringType = 'FL_Straits_MediumCrude'  # use None for no weathering -- weathering can be
                          # post-processed by the TAP viewer for instantaneous
                          # releases

#If ReceptorType is Grid, you need these, it defines the GRID
## this one to include the Bahamas
class Grid:
	pass
Grid.min_lat = 22.0 # decimal degrees
Grid.max_lat = 45.0
Grid.min_long = -98.0
Grid.max_long = -65.0
Grid.num_lat = 138
Grid.num_long = 165


TrajectoriesPath = "Trajectories" # relative to RootDir
#TrajectoriesRootname = "FlStr_Traj"

CubesPath = "Cubes"
CubesRootNames = ["FlSt" for i in StartTimeFiles] # built to match the start time files

CubeStartSitesFilename = os.path.join(RootDir, "FlStrStartSites.txt")

CubeStartSites = [x.split("#", 1)[0].strip() for x in open(CubeStartSitesFilename).readlines()]
CubeStartSites = [x for x in CubeStartSites if x]


## TAP Viewer Data (for SITE.TXT file)
##
TAPViewerSource = "Q:\EastCoast_TAP\TAP_Viewer" # where the TAP view, etc lives.

MapName = "East Coast"
MapFileName, MapFileType = ("EastCoast.bna", "BNA")
days = [1, 3, 5, 7, 10, 15, 20, 30, 50, 70, 90, 120]
                            
OutputTimes = [24*i for i in days] # output times in hours(calculated from days
#OutputUserStrings = ["5 hours",
#                     "6 hours",
#                     "12 hours",
#                     "1 day",
#                     "2 days",
#                     "3 days",
#                     ]
OutputUserStrings = ["%i days"%i for i in days]

TrajectoryRunLength = OutputTimes[-1]

PresetLOCS = ["500 barrels", "1000 barrels", "2000 barrels", "5000 barrels", "10000 barrels"]
PresetSpillAmounts = ["75000 barrels", "1000 barrels", "10000 barrels", "50000 barrels"]

## setup for the Viewer
TAPViewerPath = "TAP_Viewer"
