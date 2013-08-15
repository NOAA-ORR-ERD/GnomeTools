"""
TAP Setup.py

Master set up script for a TAP run

All the data required to set up and build TAP cubes + site.txt file should be in here

This is an example, using a San Diego Save file

"""

import sys, os, datetime

RootDir =  os.path.split(__file__)[0]
print "Loading TAP Setup data from:", RootDir

## Full Path to GNOME executable
GNOME = os.path.join(RootDir, "gnome")
if sys.platform == 'darwin':
    GNOME += ".app"
else: # assume it's Windows
    GNOME += ".exe"

## Cube Builder Data
ReceptorType = "Grid" # should be either "Grid" or "Polygons" (only grid is supported at the moment)

## CubeType can be "Volume" or "Cumulative"
##   "Volume" counts the maximum amount of oil in a given receptor since
##            the beginning of the spill
##   "Cumulative" counts all the oil that has passed through a given receptor
##                since the beginning of the spill
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
# These are used to comute the possibel time files. The format is:
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

# And example time range:
DataStartEnd = (datetime.datetime(2004, 1, 1),
                datetime.datetime(2007, 12, 31)
               )

##
##DataGaps = (( datetime.datetime(1999, 1, 3), datetime.datetime(1999, 2, 2)  ),
##            ( datetime.datetime(2001, 4, 5, 12), datetime.datetime(2001, 4, 20) ),
##            )

DataGaps = ( ) # none in this case

# specification for how you want seasons to be defined:
#  a list of lists:
#  [name, (months) ]
#    name is a string for the season name  
#    months is a tuple of integers indicating which months are in that season

Seasons = [["Winter", (12, 1, 2, 3, 4) ],
           ["Summer", (5, 6, 7, 8, 9, 10, 11) ],
           ]

StartTimeFiles = [(os.path.join(RootDir, s[0]+'Starts.txt'), s[0]) for s in Seasons]

# number of start times you want in each season:
NumStarts = 10 # small number for the test 

# Length of release in hours
ReleaseLength = 12 # hours

ModelTimeStep  =  30 # minutes

# name of the GNOME SAV file you want to use
SaveFileName = os.path.join(RootDir, "SanDiego.SAV") 

# number of Lagrangian elements you want in the GNOME run
NumLEs = 500 # should be over 1000, in general, fewer for the sample

LE_Windage = (0.01, 0.04) # windage as a fraction (0.0 not allowed?)
                            # set to None for the default

# we only have "MediumCrude"  in the data for now (see OilWeathering.py)
# OilWeatheringType = 'MediumCrude'
OilWeatheringType = None  # use None for no weathering -- weathering can be
                          # post-processed by the TAP viewer for instantaneous
                          # releases

#If ReceptorType is Grid, you need these, it defines the GRID
class Grid:
	pass
Grid.min_lat = 32.533333 # decimal degrees
Grid.max_lat = 32.738889
Grid.min_long = -117.3
Grid.max_long = -117.092
Grid.num_lat = 40
Grid.num_long = 30

TrajectoriesPath = "Trajectories" # relative to RootDir
TrajectoriesRootname = "SDT_Traj"

CubesPath = "Cubes"
CubesRootNames = ["SDT_" for i in StartTimeFiles] # built to match the start time files
                                                  # cubes Root name can only be 4 characters

CubeStartSitesFilename = os.path.join(RootDir, "SD_StartSites.txt")

CubeStartSites = [x.split("#", 1)[0].strip() for x in open(CubeStartSitesFilename).readlines()]
CubeStartSites = [x for x in CubeStartSites if x]


## TAP Viewer Data (for SITE.TXT file)
##
TAPViewerSource = os.path.join(RootDir, "TAP_ViewerSource") 


MapName = "San Diego"
MapFileName, MapFileType = ("SanDiegoMap.bna", "BNA")
#days = [ 30, 60, 90, 120]
outputdays = range(1, 5)
OutputTimes = [24*i for i in outputdays] # output times in hours(calculated from days
#OutputUserStrings = ["5 hours",
#                     "6 hours",
#                     "12 hours",
#                     "1 day",
#                     "2 days",
#                     "3 days",
#                     ]

# this to auto-generate
OutputUserStrings = ["%i days"%i for i in outputdays]

TrajectoryRunLength = OutputTimes[-1]

PresetLOCS = ["1 barrels", "10 barrels", "100 barrels", "1000 barrels"]
PresetSpillAmounts = ["100 barrels", "1000 barrels", "10000 barrels"]

## setup for the Viewer
TAPViewerPath = "SD_Test_TAP" # relative to RootDir

