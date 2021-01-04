"""
TAP Setup.py
Master set up script for a TAP run
All the parameters required to set up and build TAP cubes + site.txt file should be in here

"""
import os
from netCDF4 import Dataset, MFDataset, num2date
import numpy as np
import datetime

RootDir =  os.path.split(__file__)[0]
print "Loading TAP Setup data from:", RootDir

#******************************************************************************
# Set some paths for output data and TAP viewer
#******************************************************************************
TrajectoriesPath = "st_simons_output\Trajectories" # relative to RootDir
ImagesPath = "st_simons_output\Images"
CubesPath = "st_simons_output\Cubes_n" 
TAPViewerSource = RootDir # where the TAP viewer, etc lives.
TAPViewerPath = os.path.join('st_simons_output',"TapView")

#******************************************************************************
# Forcing data and Map
#******************************************************************************
Data_Dir = os.path.join(RootDir,'st_simon_input')

#If not None a GridCurrentMover is added to model
curr_fn = os.path.join(Data_Dir,'TidalCurrent.cur')
tide_fn = os.path.join(Data_Dir,'Entrancesouth_C_2019-cst.txt')
curr_topo = None

#If not None a GridWindMover or WindMover is added to model
wind_fn = None
wind_topo = None
wind_type = 'Point' #Point or Grid
wind_data = ('10','45') #mag, direction


#If not None a random mover is added to the model
diff_coef = 10000

# time span of your data set - based on currents (use MFdataset for filelist)
tstart = datetime.datetime(2019,10,1,0,0)
tend = datetime.datetime(2019,10,31,0,0)
DataStartEnd = (tstart,tend)

DataGaps = ( ) #not sure this is well supported?

MapName = "StSimonSound"
MapFileDir = os.path.join('.',Data_Dir)
MapFileName, MapFileType = ('StSimonsSound.bna', "BNA")

#******************************************************************************
# Params for model start times, spill durations, output time steps
#******************************************************************************
TimeSeries = None #see older Setup_TAP for docs -- do we still need?

#If you dont't want random starts use this to set evenly spaced start times
#on this hourly interval from the model start time 
EvenStarts = 1

# specification for how you want seasons to be defined:
#  a list of lists:
#  [name, (months) ] eg. ['Winter', [12,1,2]]
Seasons = [
          ["AllYear", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]]
          # ["Winter", [12, 1, 2, 3, 4, 5]],
          # ["Summer",  [6, 7, 8, 9, 10, 11 ]], 
          ]      
StartTimeFiles = [(os.path.join(RootDir, s[0]+'Starts.txt'), s[0]) for s in Seasons]
          
# Number of model runs you want in each season:
NumStarts = 24

# section for defining already run Trajectory files for BuildCubes
RunFiles = [] #???

#Output time in Days
#days = [1, 2, 3]
#Units = ' days'
#OutputTimes = [24*i for i in days] # output times in hours(calculated from days
#OutputUserStrings = [str(i) + Units for i in days]

#Output time in Hours
OutputTimes = range(3,73,12)
Units = ' hours'
OutputUserStrings = [str(i) + Units for i in OutputTimes]

GnomeOutputTimestepHours = .25 

TrajectoryRunLength = OutputTimes[-1]

#******************************************************************************
# Params for spill
#******************************************************************************
ReleaseLength = 0  # in hours (0 for instantaneous)
NumLEs = 1000 # number of Lagrangian elements you want in the GNOME run
                      
# we only have "MediumCrude"  in the data for now (see OilWeathering.py)
# TODO: add option to run GNOME with weathering                      
OilWeatheringType = None 

#Spill locations
CubeStartSitesFilename = os.path.join(RootDir, "amy_sites.txt")
CubeStartFilter = []   # January ?? no idea what this is for

PresetLOCS = ["5 barrels", "10 barrels", "20 barrels"]
PresetSpillAmounts = ["1000 barrels", "100 barrels"]

#******************************************************************************
# Params for setting up grid and building cubes
#******************************************************************************
## Cube Builder Data
ReceptorType = "Grid" # should be either "Grid" or "Polygons" (only grid is supported)
CubeType = "Volume" # should be either "Volume" or "Cumulative"
## CubeDataType options are: 'float32', 'uint16', or 'uint8'
##   float32 gives better precision for lots of LEs
##   uint8 saves disk space -- and is OK for about 1000 LEs
##   uint16 is a mid-point -- probably good to 10,000 LEs or so
CubeDataType = 'float32' 

class Grid:
	pass
Grid.min_lat = 30.5 # decimal degrees
Grid.max_lat = 31.25
Grid.dlat = 0.01       #  

Grid.min_long = -81.6
Grid.max_long = -81.2
Grid.dlong = 0.01       # 17km wide cells at 70N, 15 at 75N, 23 at 65N

Grid.num_lat = int(np.ceil(np.abs(Grid.max_lat - Grid.min_lat)/Grid.dlat) + 1)
Grid.num_long = int(np.ceil(np.abs(Grid.max_long - Grid.min_long)/Grid.dlong) + 1)

#******************************************************************************
# Calculated params - probably don't need to mess with
#******************************************************************************
CubesRootNames = ["Arc_" for i in StartTimeFiles] # built to match the start time files


# this code reads the file
CubeStartSites = [x.split("#", 1)[0].strip() for x in open(CubeStartSitesFilename).readlines()]
CubeStartSites = [x for x in CubeStartSites if x]
