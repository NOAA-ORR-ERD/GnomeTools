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
NumTrajMachines = 2

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

# time span of your data set
DataStartEnd = (datetime.datetime(1994, 8, 1, 1),
                datetime.datetime(1994, 8, 31, 23)
                )

DataGaps = ( )


# specification for how you want seasons to be defined:
#  a list of lists:
#  [name, (months) ]
#    name is a string for the season name  
#    months is a tuple of integers indicating which months are in that season

# could do 
Seasons = [["All_year", range(1,13) ],
               ]
# # example for specifying season
# Seasons = [["Winter", [11, 12, 1, 2, 3]],
#            ["Summer",[4, 5, 6, 7, 8, 9, 10]],
#           ]

# You don't need to do anything with this
StartTimeFiles = [(os.path.join(RootDir, s[0]+'Starts.txt'), s[0]) for s in Seasons]

# number of start times you want in each season:
#NumStarts = 500
NumStarts = 10

# # Length of release in hours
# ReleaseLength = 24 * 90 #24 hrs * XX days

# Instantanious
ReleaseLength = 0

# name of the GNOME SAV file you want to use
# note: GNOME locks it (for a very brief time when loading) 
# which means you need a separate copy for each
# instance of GNOME you want to run (OR just don't start multiple GNOMES too quickly)
PyGnome_script = "script_ucla"

# number of Lagrangian elements you want in the GNOME run
NumLEs = 100
                            
# we only have "MediumCrude"  in the data for now (see OilWeathering.py)
OilWeatheringType = None
# OilWeatheringType = 'FL_Straits_MediumCrude'  # use None for no weathering -- weathering can be
#                           # post-processed by the TAP viewer for instantaneous
#                           # releases

#If ReceptorType is Grid, you need these, it defines the GRID
## this one to include the Bahamas
class Grid:
	pass
Grid.min_lat = 32.0 # decimal degrees
Grid.max_lat = 36.0
Grid.min_long = -123.0
Grid.max_long = -117.0
Grid.num_lat = 20
Grid.num_long = 30


TrajectoriesPath = "Trajectories" # relative to RootDir
#TrajectoriesRootname = "FlStr_Traj"

CubesPath = "Cubes"
CubesRootNames = ["SB_" for i in StartTimeFiles] # built to match the start time files

CubeStartSitesFilename = os.path.join(RootDir, "SB_sites.txt")

# this code reads the file
#CubeStartSites = [x.split("#", 1)[0].strip() for x in open(CubeStartSitesFilename).readlines()]
#CubeStartSites = [x for x in CubeStartSites if x]


## TAP Viewer Data (for SITE.TXT file)
##
TAPViewerSource = "TAP_Viewer" # where the TAP view, etc lives.

MapName = "SC Bight"
MapFileName, MapFileType = ("coast.bna", "BNA")

days = [1, 3, 5, 7, 10, 15, 20, 30, 50, 70, 90, 120]
OutputTimes = [24*i for i in days] # output times in hours(calculated from days

# output time in hours
OutputTimes = [6, 12, 24, 48, 72]

OutputUserStrings = ["6 hours",
                     "12 hours",
                     "1 day",
                     "2 days",
                     "3 days",
                     ]

# this is calculated from the OutputTimes
TrajectoryRunLength = OutputTimes[-1]

PresetLOCS = ["5 barrels", "10 barrels", "20 barrels"]
PresetSpillAmounts = ["200 barrels", "10 barrels"]

## setup for the Viewer
TAPViewerPath = "TAP_Viewer"
