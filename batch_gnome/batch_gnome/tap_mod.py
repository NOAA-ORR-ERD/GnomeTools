#!/usr/bin env python

"""
module TAP_mod

Contains assorted functions, etc for computing TAP stuff

Functions include:
read_binLE(filename):


ReadTrajectory(filename,CheckFlag = 1,HeaderData = {})

    Returns a Numpy array of LE positions, and a dictionary of the Header data.

transform(data, n=1.5):

    converts a NumPy array of Floats to 8 bit integer for TAP II binary cubes                
    This is reversed by: (Cube**n)/(255**n) * num_LE;

reverse_transform(data, n=1.5):
    convert from 8 bit integer NumPy array to Float array from TAP II binary cubes

CompTAPIICube(FileList,OutputTimes):
    Computes a TAP II cube from the list of Files in FileList
    at the output times listed in OutputTimes ( a sequence of hours from 0)


Classes include:

Receptor_set:

Constructor:
    Receptor_set(sites):
        sites is a list of NumPy arrays of Floats


Attributes:
    sites: a list of (N,2) NumPy arrays of Floats, each is a receptor site
    bb:    a list of (4,) NumPy arrays of Floats, each is the bounding box
           for the corasponding site. : array((max_x,min_x,max_y,min_y),Float)

Methods:
    comp_hits(LEs,Hit_Table):
       A function that computes the intersections of the LEs

    LEs is a M X N X 2 NumPy array (of Floats).
       M is the number of time steps
       N is the number of LEs

    Hit_Table is a NumPy array of np.uint8 of size (Num_LEs, Num_sites),
    it hold the values of 0 or 1, depending on whether a given LE has hit a given sight.
    ***Hit_Table is ALTERED by this function!!!***


"""

##from Tap_ext.check_receptors import hit_test, Grid_hit_test
#import check_receptors as CR
##hit_test = check_receptors.hit_test
##Grid_hit_test = check_receptors.Grid_hit_test

#from numpy import *
#from Hazmat.TAP.TAP_ext import NumericExtras
#byteswap = NumericExtras.byteswap
#changetype = NumericExtras.changetype
import sys
import time
#import Geometry

import string # depreciated, but it's there

import numpy as np

class TAPError(Exception):
    pass


def read_binLE(filename): # A faster version (hopefully!)

    """
    Function that reads a binary LE file, and returns a (NX2) NumPy array
    of LE coordinates, in Lat-Long format, and a list of Beach codes
    (-50 = beached)
    
    NOTE: the binary LE files are assumed to be Big-endian (Mac format)
    """
    
    
    import struct
    
    # test for endian format:
    if sys.byteorder == 'big':
        endian = 'B'
    elif sys.byteorder == 'little':
        endian = 'L'
        
        
    file = open(filename,'rb')
    
    header_format = '>10shhhhhffl'
    header =  file.read(struct.calcsize(header_format))
    
    (version1,day,month,year,hour,minute,current_time,version2,num_LE) = struct.unpack(header_format,header)
    
    # I'm not using struct anymore, but these data are handy to know:
    # le_format = '>ffffflll'
    #        (Lat,Long,Release_time,AgeWhenReleased,BeachHeight,nMap,pollutant,WindKey)=struct.unpack(le_format,record)
    
    data  = np.fromfile(file, dtype=np.uint8) # get all the data as bytes
    data.shape = (num_LE, 8*4) # reshape into records
    
    Beached = data[:,4*4:5*4].copy() # extract the Beach code.
    Beached = Beached.view(dtype=np.float32).reshape((num_LE,)) # reconfigure as C float

    data = data[:,:4*2].copy() # extract the lat and long
    data = data.view(dtype=np.float32).reshape( (num_LE,2) ) # reconfigure as C float

    # Do we need to byteswap?
    if sys.byteorder == 'little':
        Beached.byteswap(True)
        data.byteswap(True)
        
    # Switch Lat and Long and Assume positive is the western hemisphere
    LEs = np.zeros_like(data)
    LEs[:,0] = -data[:,1]
    LEs[:,1] = data[:,0]
    
    return LEs, Beached

    
def ReadTrajectory_Old(filename, CheckFlag = 1,HeaderData = {}): # Read a binary trajectory file
    """
    Reads the Batch-GNOME trajectory file stored in filename.
    
    Call:
    (Trajectory,(NumTimesteps,NumLEs),HeaderData,flags) = ReadTrajectory(filename,CheckFlag,HeaderData)
    
    Arguments:
    *filename : the name (full path is a good idea) of the trajectroy file
    *Checkflags: and integer flag:
      -if Checkflag = 0, don't read the header
      -if Checkflag = 1, read eveything
      -if Checkflag = 2, read only the Header
      - default is 1
    *HeaderData: a dictionary of the header data. Only needed in Checkflag == 0
    
    Returns:
    *the LE positions (Trajectory) as a (Num_times X Num_LEs X 2)
    NumPy array of Python Floats (C doubles)
    
    *(NumTimesteps,NumLEs) as a tuple of integers
    
    *HeaderData as a dictionary of the data in the header
    
    *Flags is a (NumTimesteps X NumLEs) NumPy array of Unsigned single byte integers
    
    NOTE: the text header uses DOS style line endings (\r\n)
    """
    file = open(filename,'rb') # open as a binary file
    
    if CheckFlag:
        # read the header:
        # Note that I'm only bothering to parse out the info I really need.
        HeaderData = {}
        line = string.strip(file.readline())
        NumLines =  1
        while line[:13] != "[BINARY DATA]": # do until the end of the header
            if string.strip(line)[0] == "[":
                key = line[1:string.find(line,']')]
                if not ( key == "FILE INFORMATION"
                         or key == "RUN PARAMETERS"
                         or key == "BINARY FORMAT"):
                    raise TAPError("%s doesn't look like a TAPRUN file to me!!"%
                                   filename)
            else:
                type,data = string.split(line,':',1)
                HeaderData[type] = string.strip(data)
            line = string.strip(file.readline())
            NumLines = NumLines + 1
        HeaderData["NumLines"] = NumLines
        
        # Check the header info
        if HeaderData["File type"] != "TAPRUN":
            raise("File Type is wrong!!")
        if string.strip(string.split(HeaderData["Run duration"])[1]) != "hours":
            raise("Run Duration should be in hours!! in \n%s"%filename)
        if not (HeaderData['Endian'] == "big" or HeaderData['Endian'] == "little"):
            raise TAPError('Endian needs to be either "big" or "little" in \n%s'%filename)
            
        # Change some Header Data to numbers
        HeaderData["Number of LEs"] = int(HeaderData["Number of LEs"])
        HeaderData['Number of output steps'] = int(HeaderData['Number of output steps'])
        HeaderData["Output time step in seconds"] = float(HeaderData["Output time step in seconds"])
        
        if not ((HeaderData["Bit flags"] == "CHAR notReleased beached offMaps evaporated notOnSurface")
                and (HeaderData["Longitude"] == "LONG" and HeaderData["Latitude"] == "LONG")):
            raise  TAPError("This is not a binary format I know how to deal with")
    else:
        for i in range(HeaderData['NumLines']):
            file.readline()
            ## now read the binary part:
            
    if CheckFlag < 2:
        NumTimesteps = HeaderData['Number of output steps']
        NumLEs = HeaderData['Number of LEs']
        
        NumBytes = 9 # number of bytes in a record
        
        # Read the file
        data = fromstring(file.read(NumBytes*NumTimesteps*NumLEs),uint8)
        file.close()
        
        try:
            data.shape = (NumTimesteps,NumLEs,NumBytes)
        except ValueError:
            raise TAPError(" The size of the data in the Trajectory file:\n%s\n Does not match the header"%
                           filename )

        # extract LE data:
        LEs = data[:,:,:8].copy()
        # extract the flags
        flags = data[:,:,8].copy()
        del data

        #LEs = fromstring(LEs.tostring(),np.int32)
        changetype(LEs,np.int32)
        
        if HeaderData['Endian'] != sys.byteorder:
            byteswap(LEs)
            
        LEs.shape = (NumTimesteps,NumLEs,2)
        
        # convert to numpy float32 (C float)
        Trajectory = LEs.astype(np.float32)
        
        # convert to degrees (in place)
        divide(Trajectory, array((1e6), np.float32), Trajectory)
        
        
        flags.shape = (NumTimesteps, NumLEs)
        return (Trajectory, (NumTimesteps,NumLEs), HeaderData,flags)
    else:
        return (None, (None, None), HeaderData, None)


def ReadTrajectory(filename, CheckFlag = 1, HeaderData = {}): # Read a binary trajectory file
    ## fixme -- this reads the 3-d version, but tosses the z coord. 
    """
    Reads the Batch-GNOME trajectory file stored in filename.
    
    Call:
    (Trajectory,(NumTimesteps,NumLEs),HeaderData,flags) = ReadTrajectory(filename,CheckFlag,HeaderData)
    
    Arguments:
    *filename : the name (full path is a good idea) of the trajectroy file
    *Checkflags: and integer flag:
      -if Checkflag = 0, don't read the header
      -if Checkflag = 1, read eveything
      -if Checkflag = 2, read only the Header
      - default is 1
    *HeaderData: a dictionary of the header data. Only needed in Checkflag == 0
    
    Returns:
    *the LE positions (Trajectory) as a (Num_times X Num_LEs X 2)
    NumPy array of numpy float32 (C floats)
    
    *(NumTimesteps,NumLEs) as a tuple of integers
    
    *HeaderData as a dictionary of the data in the header
    
    *Flags is a (NumTimesteps X NumLEs) NumPy array of Unsigned single byte integers
    
    NOTE: the text header uses DOS style line endings (\r\n)
    """
    file = open(filename,'rb') # open as a binary file
    
    if CheckFlag:
        # read the header:
        # Note that I'm only bothering to parse out the info I really need.
        HeaderData = {}
        line = file.readline().strip()
        NumLines =  1
        while line[:13] != "[BINARY DATA]": # do until the end of the header
            if line.strip()[0] == "[":
                key = line[1:string.find(line,']')]
                if not ( key == "FILE INFORMATION"
                         or key == "RUN PARAMETERS"
                         or key == "BINARY FORMAT"):
                    raise TAPError("%s doesn't look like a TAPRUN file to me!!"%
                                   filename)
            else:
                type,data = string.split(line,':',1)
                HeaderData[type] = string.strip(data)
            line = string.strip(file.readline())
            NumLines = NumLines + 1
        HeaderData["NumLines"] = NumLines
        
        # Check the header info
        if HeaderData["File type"] != "TAPRUN":
            raise("File Type is wrong!!")
        if string.strip(string.split(HeaderData["Run duration"])[1]) != "hours":
            raise("Run Duration should be in hours!! in \n%s"%filename)
        if not (HeaderData['Endian'] == "big" or HeaderData['Endian'] == "little"):
            raise TAPError('Endian needs to be either "big" or "little" in \n%s'%filename)
            
        # Change some Header Data to numbers
        HeaderData["Number of LEs"] = int(HeaderData["Number of LEs"])
        HeaderData['Number of output steps'] = int(HeaderData['Number of output steps'])
        HeaderData["Output time step in seconds"] = float(HeaderData["Output time step in seconds"])
        
        ##fixme: this should really check the binary format better!
        if not ( (HeaderData["Bit flags"] == "CHAR notReleased beached offMaps evaporated notOnSurface")
                  and (HeaderData["Longitude"] == "LONG")
                  and (HeaderData["Latitude"] == "LONG")
                  and (HeaderData["z"] == "DOUBLE") ):
            raise  TAPError("This is not a binary format I know how to deal with")
    else:
        for i in range(HeaderData['NumLines']):
            file.readline()
            ## now read the binary part:
            
    if CheckFlag < 2:
        NumTimesteps = HeaderData['Number of output steps']
        NumLEs = HeaderData['Number of LEs']
        
        DataType = np.dtype([("long","i4"), ("lat", "i4"), ("z", "f8"), ("flag","u1")])
       
        # Read the file
        #data = fromstring(file.read(NumBytes*NumTimesteps*NumLEs),np.uint8)
        data = np.fromfile(file, DataType)
        file.close()
        try:
            data.shape = (NumTimesteps, NumLEs)
        except ValueError:
            raise TAPError(" The size of the data in the Trajectory file:\n%s\n Does not match the header"%
                           filename )

        # extract LE data:
#        LELat = data['lat']
#        LELong = data['long']
#        LELat.shape = LELong.shape = (NumTimesteps,NumLEs,1)
        
        LEs = np.empty((NumTimesteps, NumLEs, 2), dtype=np.int32)

        LEs[:,:,0] = data['long']
        LEs[:,:,1] = data['lat']

#        LEs = np.c_[(LELong, LELat)]

        # extract the flags
        flags = data['flag'].copy()
        del data

        #LEs = fromstring(LEs.tostring(),np.int32)
        #changetype(LEs,np.int32)
        
        if HeaderData['Endian'] != sys.byteorder:
            LEs.byteswap(True)
            
        #LEs.shape = (NumTimesteps,NumLEs,2)
        
        # convert to python Floats (C doubles)
        Trajectory = LEs.astype(np.float32)
        
        # convert to degrees (in place)
        np.divide(Trajectory, np.array((1e6), np.float32),Trajectory)
        
        flags.shape = (NumTimesteps, NumLEs)
        return (Trajectory, (NumTimesteps, NumLEs), HeaderData, flags)
    else:
        return (None, (None, None), HeaderData, None)
        
        
class Receptor_set:
    """
    Constructor:
    Receptor_set(sites):
        sites is a list of NumPy arrays of Floats
    
    
    Attributes:
    *sites: a list of (N,2) NumPy arrays of Floats, each is a receptor site
    *bb:    a list of (4,) NumPy arrays of Floats, each is the bounding box
           for the corasponding site. : array((max_x,min_x,max_y,min_y),Float)
    
    Methods:
    
    *comp_hits(LEs,Hit_Table):
       A function that computes the intersections of the LEs
    
    *LEs is a M X N X 2 NumPy array (of Floats).
       M is the number of time steps
       N is the number of LEs
    
    Hit_Table is a NumPy array of np.uint8 of size (Num_LEs, Num_sites),
    it hold the values of 0 or 1, depending on whether a given LE has hit a given sight.
    ***Hit_Table is ALTERED by this function!!!***
    
    """
    
    def __init__(self, polyset):
        # polyset is a set of polygons, as returned by read_bna
        self.sites = []
        self.bb = []
        for type,site in polyset:
            self.sites.append(site)
            self.bb.append(array((max(site[:,0]), min(site[:,0]), max(site[:,1]), min(site[:,1])),np.float))
    def NumSites(self):
        return len(self.sites)
    def __len__(self):
        return len(self.sites)
        
    def comp_hits(self,LEs,Hit_Table,StartTimestep = 0):
        import check_receptors as CR
        CR.hit_test(LEs,self.sites,self.bb,Hit_Table,StartTimestep)
        return None
        
## conversion functions for binary TAP II cubes

def transform(data, num_LE, n = 1.5, dtype=np.uint8):
    """
    Convert number of LEs to 8 bit integer for storage in Cubes               
    This is reversed by: (Cube**n)/(255**n) * num_LE;
    input data is a Numpy array

    """
    # make sure the input is float32
    data = np.asarray(data, dtype=np.float32)

    if dtype is np.float32:
        data = data / num_LE 
    elif dtype is np.uint8:
        data = np.ceil((data / num_LE)**(1.0/n)*255)  
    elif dtype is np.uint16:
        data = np.ceil((data / num_LE)**(1.0/n)*65535)  

    return data.astype(dtype)
    
def reverse_transform(data, num_LE = 1000, n = 1.5, dtype=np.uint8):
    """
    NOTE: this isn't used, but is handy for testing
    also, it isn't exact, because of the use of "ceil" in both cases 
    """
    # convert from integer NumPy array to float32 array
    if dtype is np.uint8:
        max_val = 255
    elif dtype is np.uint16:
        max_val = 65535
    
    data = data.astype(np.float)
    data = np.ceil((data/max_val)**n * num_LE)

    return data

    
def CompTAPIICube(FileList, OutputTimes, Receptors):
    # read the header of the first one:
    (junk,junk,HeaderData,junk) = ReadTrajectory(FileList[0],2)
    
    NumTimesteps = HeaderData['Number of output steps']
    NumLEs = HeaderData['Number of LEs']
    TimeStep = HeaderData["Output time step in seconds"]

    OutputSteps = np.array(([0] + OutputTimes),np.int)*60*60/int(TimeStep)
    OutputSteps2 = np.array(OutputTimes,np.int)*60*60/int(TimeStep)
    
    # Allocate Cube
    NumSpills,NumSites,NumTimes = len(FileList),len(Receptors),len(OutputTimes)

##    ##Old version    
##    Cube = np.zeros((NumTimes,NumSites,NumSpills),np.uint8)
##    start = time.time()
##    for SpillNum in range(NumSpills):
##        print "computing spill number %i"%(SpillNum,)
##        HitTable = np.zeros((NumLEs,NumSites),np.int) - 1
##        (Trajectory,(NumTimesteps,NumLEs),HeaderData,flags) = ReadTrajectory(FileList[SpillNum],1,HeaderData)
##        for StepNum in range(len(OutputTimes)):
##            Start,Finish = OutputSteps[StepNum],OutputSteps[StepNum+1]+1
##            Receptors.comp_hits(Trajectory[Start:Finish,:,:],HitTable,Start)
##            Cube[StepNum,:,SpillNum] = transform(sum((HitTable >= 0),0),NumLEs)
##    print "Version 1 took %f seconds"%(time.time()-start)

    Cube = np.zeros((NumTimes,NumSites,NumSpills),np.uint8)
    #print "cube is : %i MB"%(Cube.size * Cube.itemsize / 1024 / 1024)
    #raw_input("about to start building cubes: hit enter to continue")
    start = time.time()
    # pre-create the Hit Table
    HitTable = np.empty((NumLEs,NumSites),np.int)
    for SpillNum in range(NumSpills):
        #print "computing spill number %i"%(SpillNum,)
        #print "NumLEs, NumSites:", NumLEs, NumSites
        HitTable[:] = -1 # set HitTable to "nothing"
        #raw_input("created the hit table:  hit enter to continue")
        (Trajectory,(NumTimesteps,NumLEs),HeaderData,flags) = ReadTrajectory(FileList[SpillNum],1,HeaderData)
        print "read the trajectory"
        #raw_input("About to compute hits: hit enter to continue")
        Receptors.comp_hits(Trajectory, HitTable)
        del Trajectory # don't need it anymore -- do need the memory
        #raw_input("compiling the hit table: hit enter to continue")
        for StepNum in range(len(OutputTimes)):
            #print "Step Num: %i of %i"%(StepNum, len(OutputTimes) )
            hitSites = ((HitTable <= OutputSteps2[StepNum]) & (HitTable >= 0))
            Cube[StepNum,:,SpillNum] = transform(sum(hitSites), NumLEs)
        #raw_input("About to delete extra arrays: hit enter to continue")
        del hitSites
        #raw_input("deleted extra arrays: hit enter to continue")
##    ## test to see if they are the same:
##    print "Version 2 took %f seconds"%(time.time()-start)
##    if alltrue(Cube.flat == Cube2.flat):
##        print "They are all equal"
##    else:
##        print "They are NOT all equal"
##    return Cube, Cube2
    return Cube

def CompThicknessCubeOld(FileList, OutputTimes, Grid, Weather = None, data_type='byte'):
    """

    CompThicknessCube computes the average thickness of the oil over
    each receptor site. It only works for grid receptors.

    If Weather is not None it must be a tuple of A,B values, where A and
    B are the coefficient of a weathering curve:

    Volume = InitialVolume * (A*log( TimeSinceRelease) ) + B)

    (I think, I didn't write this down when I wrote the code!)

    If Weather is None, then there is no change in volume.
        
    """
    # G = Grid
    # min_long
    # max_long, min_lat, max_lat,num_lat,num_long = Grid
    # dlat = (G.max_lat-G.min_lat) / G.num_lat
    # dlong = (G.max_long-G.min_long) / G.num_long

    min_long = Grid.min_long
    min_lat  = Grid.min_lat
    num_lat  = Grid.num_lat
    num_long = Grid.num_long
    dlat     = Grid.dlat
    dlong    = Grid.dlong
    NumSites = Grid.num_cells

    if Weather:
        A, B = Weather
    ## read the header of the first trajectory file: It is assumed that
    ## the others will match..no check is made to assure this, but it
    ## will crash if anything is very wrong.
    else:
        pass
        #print "no weathering coeff"
        
    (junk,junk,HeaderData,junk) = ReadTrajectory(FileList[0],2)

    NumTimesteps = HeaderData['Number of output steps']
    NumLEs = HeaderData['Number of LEs']
    TimeStep = HeaderData["Output time step in seconds"]
    TimeStepHours = TimeStep / 3600
    OutputSteps = np.array( [0] + OutputTimes, np.int32 )*60*60/int(TimeStep) # in units of time step
    # Allocate the Cube
    NumSpills, NumTimes = len(FileList), len(OutputTimes)
    Cube = np.zeros((NumTimes,NumSites,NumSpills), np.float32)

    start = time.time() # just for timing how long it takes to run

    ## Loop through each individual trajectory
    for SpillNum in range(NumSpills):
        #print "computing spill number %i"%(SpillNum,)
        ## Load the trajectory data (see function in this file)
        ## Trajectory is a large (NumTimesteps X NumLEs X 2) Numeric array with the positions of all the lEs at all time steps
        ## flags is an array (NumTimesteps X NumLEs) of the one-byte flags associated with the LEs
        (Trajectory,(NumTimesteps,NumLEs),HeaderData,flags) = ReadTrajectory(FileList[SpillNum],1,HeaderData)
        ## Hits = np.zeros((NumSites), np.int32) # note: could this be an int16?
        VolTable = np.zeros((NumSites), np.float64) # this will store the Maximum volume in each grid box.
        released = np.logical_not(flags & 1) # make a "released" flag.
        #del flags

        ## Step through the Cube output time steps
        for step in xrange( len(OutputSteps) -1 ):
            ## step through the Trajectory time steps between each Cube Timestep
            for t in xrange(OutputSteps[step], OutputSteps[step+1]):
                Vol = np.zeros((NumSites), np.int32) # this stores the volume (number of LEs)
                LEs = Trajectory[t,:,:] # the slice out of the Trajectory array for this time step
                LE_flags = flags[t,:]
                if Weather:
                    putmask(ReleasedSteps, (released[t,:] & (ReleasedSteps == -1) ), t ) # set release time step for the ones just released
                    WeatherTime = (t - ReleasedSteps) * TimeStepHours
                    ValidLEs = np.nonzero(WeatherTime >= 1)
                    Weathered_Vol_1 = A*log( np.take(WeatherTime,ValidLEs) ) + B # compute just the ones over one hour
                    Weathered_Vol =  np.zeros((NumLEs,), np.float)
                    put( Weathered_Vol, ValidLEs, Weathered_Vol_1)
                    putmask( Weathered_Vol, WeatherTime < 1, 1 )
                    # just to make sure
                    if len(np.nonzero(Weathered_Vol < 0)) > 0 or len(np.nonzero(Weathered_Vol > 1)) > 0:
                        raise TAPError("Weathered Volume out of range")
                else:
                    Weathered_Vol = np.ones((NumLEs,), np.float)
                    #remove the "not released LEs"
                    Weathered_Vol[(LE_flags & 1).astype(bool)] = 0.0
                ### What boxes are the LEs in?
                i = ( (LEs[:,0] - min_long) / dlong ).astype(np.int32)
                j = ( (LEs[:,1] - min_lat ) / dlat ).astype(np.int32)
                # check for whether the LEs are in the bounds of the grid at all
                InGrid =  np.nonzero( (i < num_long) &
                                      (i >= 0) &
                                      (j < num_lat) &
                                      (j >= 0) )[0] # want the indices in the first dimension
                ## the indexes of the LEs that are within the grid
                ## translated from (i,j) to a one-d array
                site_ind = (np.take(i, InGrid) * num_lat ) + np.take(j, InGrid)
                # The Volumes of the LEs that are in the grid
                LE_Vol = np.take(Weathered_Vol, InGrid)
                ## add the volume of each LE to the total volume in it's grid box
                for k in xrange(len(site_ind)):
                    Vol[site_ind[k]] += LE_Vol[k]
                ## put the Volume in the grid box only if it's larger than the volume there
                VolTable = np.maximum(Vol, VolTable)
            ## put the max volume in the Cube at this Cube time step
            #Cube[step,:,SpillNum] = transform(VolTable, NumLEs)
            Cube[step,:,SpillNum] = VolTable
    #print "cube took %s seconds to generate"%(time.time() - start)

    return Cube

def CompThicknessCubeTimestepOld(Grid, LE_positions, LE_mass=None, flags=None, flag_bitmask_to_ignore = 1+4+8+16 ):
    """
    compute how much is in each cell for only one timestep
    """
    min_long = Grid.min_long
    min_lat  = Grid.min_lat
    num_lat  = Grid.num_lat
    num_long = Grid.num_long
    dlat     = Grid.dlat
    dlong    = Grid.dlong
    NumSites = Grid.num_cells

    LEs = LE_positions
    
    Vol = np.zeros((NumSites), np.float32) # this stores the volume (number of LEs)
    NumLEs = LEs.shape[0]

    Weathered_Vol = np.ones((NumLEs,), np.float)

    ### What boxes are the LEs in?
    i = ( (LEs[:,0] - min_long) / dlong ).astype(np.int32)
    j = ( (LEs[:,1] - min_lat ) / dlat ).astype(np.int32)
    # check for whether the LEs are in the bounds of the grid at all
    InGrid =  np.nonzero( (i < num_long) &
                          (i >= 0) &
                          (j < num_lat) &
                          (j >= 0) )[0] # want the indices in the first dimension
    ## the indexes of the LEs that are within the grid
    ## translated from (i,j) to a one-d array
    site_ind = (np.take(i, InGrid) * num_lat ) + np.take(j, InGrid)
    # The Volumes of the LEs that are in the grid
    LE_Vol = np.take(Weathered_Vol, InGrid)
    ## add the volume of each LE to the total volume in it's grid box
    for k in xrange(len(site_ind)):
        Vol[site_ind[k]] += LE_Vol[k]
    
    return Vol.reshape(num_long, num_lat) 

from tap_comp_volume import comp_volume
#from cy_tap_comp_volume import comp_volume

#DDR import nc_particles
from post_gnome import nc_particles

def CompThicknessCube(FileList, OutputTimes, Grid, Weather=None):

    """
    CompThicknessCube computes the average thickness of the oil over
    each receptor site. It only works for grid receptors.
    
    Filelist is a list of netcdf file names: one for each trajectory
    
    OutputTimes is a sequence of output times, in hours, from the beginning of the run.

    Grid is a Grid object, specifying the grid parameters
    
    If Weather is not None it must be a tap_comp_volume.weather_curve object

    If Weather is None, then there is no change in volume.
        
    """

        
    ## read the header of the first trajectory file: It is assumed that
    ## the others will match..no check is made to assure this, but it
    ## will crash if anything is very wrong.
    #(junk, junk, HeaderData, junk) = ReadTrajectory(FileList[0],2)

    # read the trajectory data from the first netcdf file
    #print "**************"
    #print "getting header info from file:", FileList[0]
    print "nc_particles module:", nc_particles.__file__
    #DDR traj_file = nc_particles.nc_particle_file(FileList[0])
    traj_file = nc_particles.Reader(FileList[0])

    if traj_file.get_units('age') != 'seconds':
        raise ValueError("particle age units in netcdf file must be in seconds")
    
    #DDR NumTimesteps = traj_file.num_times
    NumTimesteps = len(traj_file.times)
    MaxNumLEs = traj_file.particle_count[:].max() 
    
    TimeStep = traj_file.times[1] - traj_file.times[0] # assume constant timestep!
    #DDR traj_file.nc.close()
    traj_file.nc.close()
    TimeStepHours = TimeStep.total_seconds() / 3600.00
    ## OutputTimes should already be in hours
    OutputSteps = (np.array([0] + OutputTimes) / TimeStepHours).astype(np.int32) # in integer units of time step
    # Allocate the Cube
    NumSpills, NumSites, NumTimes = len(FileList), Grid.num_cells, len(OutputTimes)
    ## fixme: need to make this a float cube!
    Cube = np.zeros((NumTimes,NumSites,NumSpills), np.float32)

    start = time.time() # just for timing how long it takes to run
    print OutputSteps

    ## Loop through each individual trajectory
    for SpillNum in range(NumSpills):
        #print "computing spill number %i"%(SpillNum,)
        # read new trajectory file:
        #print "working with file:", FileList[SpillNum]
        #DDR traj_file = nc_particles.nc_particle_file(FileList[SpillNum])
        traj_file = nc_particles.Reader(FileList[SpillNum])

        VolTable = np.zeros((NumSites), np.float32) # this will store the Maximum volume in each grid box.

        ## Step through the Cube output time steps
        for step in xrange( len(OutputSteps) - 1 ):
            ## step through the Trajectory time steps between each Cube Timestep
            for t in xrange(OutputSteps[step], OutputSteps[step+1]):

                LE_lat = traj_file.get_timestep_single_var(t, 'latitude')
                LE_long = traj_file.get_timestep_single_var(t, 'longitude')
                LE_positions = np.column_stack((LE_long, LE_lat))
                NumLEs = LE_positions.shape[0]
                LE_age = traj_file.get_timestep_single_var(t, 'age').astype(np.float32) / 3600.00 # age needs to be in hours
                #print "age:", LE_age
                #NOTE: for TAP -- we assume that are the particles have unit mass at the start
                #      so we don't read it from the file
                LE_mass = np.ones((NumLEs,), dtype = np.float32)
                #print "before"
                #print LE_mass
                if Weather:
                    #print "weathering the LEs"
                    LE_mass = Weather.weather(LE_mass, LE_age)
                #DDR flags = traj_file.get_timestep_single_var(t, 'flag').astype(np.uint8)
                flags = traj_file.get_timestep_single_var(t, 'status_codes').astype(np.uint8)
                Vol = comp_volume(LE_positions, LE_mass, flags, Grid)
                # keep the largest volume computed between output timesteps
                VolTable = np.maximum(Vol.flat, VolTable)
            ## put the max volume in the Cube at this Cube time step
            #Cube[step,:,SpillNum] = transform(VolTable, MaxNumLEs)
            Cube[step,:,SpillNum] = VolTable
        traj_file.close()
    #print "cube took %s seconds to generate"%(time.time() - start)
    return Cube


class Grid:
    def __init__(self, min_long, max_long, min_lat, max_lat,num_lat,num_long):
        self.min_long = min_long
        self.max_long = max_long
        self.min_lat  = min_lat 
        self.max_lat  = max_lat 
        self.num_lat  = num_lat    
        self.num_long = num_long   
        self.dlat = float(max_lat-min_lat) / num_lat
        self.dlong = float(max_long-min_long) / num_long

    def _get_num_cells(self):
        return self.num_lat * self.num_long
    num_cells = property(_get_num_cells)
    
    def WriteBNA(self, file):
        print "Writing the BNA data"
        self.dlat =  float(self.max_lat- self.min_lat ) / self.num_lat
        self.dlong = float(self.max_long-self.min_long) / self.num_long
        polynum = 1
        for i in range(self.num_long):
             for j in range(self.num_lat):
                 file.write('"#%i","1",5\n'%polynum)
                 file.write('%11.6f, %11.6f\n'%(self.min_long + self.dlong*i ,     self.min_lat + self.dlat*j))     
                 file.write('%11.6f, %11.6f\n'%(self.min_long + self.dlong*(i+1) , self.min_lat + self.dlat*j))
                 file.write('%11.6f, %11.6f\n'%(self.min_long + self.dlong*(i+1) , self.min_lat + self.dlat*(j+1)))
                 file.write('%11.6f, %11.6f\n'%(self.min_long + self.dlong*i ,     self.min_lat + self.dlat*(j+1)))
                 file.write('%11.6f, %11.6f\n'%(self.min_long + self.dlong*i ,     self.min_lat + self.dlat*j))    
                 polynum += 1
        return True

    
    def __str__(self):
        msg = []
        msg.append("Grid object:")
        msg.append("min_lat: %s"%self.min_lat)
        msg.append("max_lat: %s"%self.max_lat)
        msg.append("min_long: %s"%self.min_long)
        msg.append("max_long: %s"%self.max_long)
        msg.append("num_lat: %s"%self.num_lat)
        msg.append("num_long: %s"%self.num_long)
        msg.append("delta-lat: %s"%self.dlat)
        msg.append("delta-long: %s"%self.dlong)
        
        return "\n".join(msg)

##def Grid_hit_test(LEs, grid, Hit_Table,Start_step):
##    N_times, N_LEs = LEs.shape[:2]

##    min_long  = grid.min_long
##    max_long = grid.max_long
##    min_lat  = grid.min_lat 
##    max_lat  = grid.max_lat 
##    num_lat  = grid.num_lat    
##    num_long = grid.num_long   
##    dlat =     grid.dlat
##    dlong =    grid.dlong

##    for T_ind in range(1,N_times): # loop over timesteps
##        for LE_ind in range(N_LEs): # loop over LEs
##            LE_line = (tuple(LEs[T_ind-1,LE_ind,:]),tuple(LEs[T_ind,LE_ind,:])) # LE-movement segment
##            ((px1,py1),(px2,py2)) = LE_line
##            ## did the LE move?
##            if (LE_line[0] != LE_line[1]):
##                ## compute bounding box
##                ## what box is the first point in?
##                i0 = int(floor((LE_line[0][0] - min_long) / dlong))
##                j0 = int(floor((LE_line[0][1] - min_lat ) / dlat ))
##                ## what box is second point in?
##                i1 = int(floor((LE_line[1][0] - min_long) / dlong))
##                j1 = int(floor((LE_line[1][1] - min_lat ) / dlat ))

##                ## sort them
##                i0,i1 = min(i0,i1), max(i0,i1)
##                j0,j1 = min(j0,j1), max(j0,j1)
##                #if (LE_ind == 4 ): print "i0,i1 = %i,%i,  j0,j1 = %i,%i\n"%(i1,i1,j0,j1);

##                ## check the boxes:
##                for i in range(i0, i1+1):
##                    for j in range(j0, j1+1):
##                        p1 = (min_long + i*dlong, min_lat + j*dlat)
##                        p2 = (min_long + (i+1)*dlong, min_lat + j*dlat)
##                        p3 = (min_long + (i+1)*dlong, min_lat + (j+1)*dlat)
##                        p4 = (min_long + i*dlong, min_lat + (j+1)*dlat)

##                        s1 = SideOfLineCheck(px1,py1,px2,py2,p1[0],p1[1])
##                        #if (LE_ind == 4 ): print "i = %i, j=%i, s1 = %g"%(i,j,s1)
##                        if s1 > 0:
##                            side = 1
##                        else:
##                            side = 0
##                        for p in (p2,p3,p4):
##                            s2 = SideOfLineCheck(px1,py1,px2,py2,p[0],p[1])
##                            if (s2 > 0) <> side:
##                                site_ind = i*num_lat + j
##                                if not Hit_Table[LE_ind,site_ind]:
##                                    Hit_Table[LE_ind,site_ind] = T_ind + Start_step
##                                break
##    return None


## Utility Functions ( should these be in geometry.py? )

def SideOfLineCheck(x1,y1,x2,y2,Px,Py):
    """ Given a line segment x1,y1 to x2,y2
    it checks to see if point Px,Py is to the right
    or to the left of the line segment looking from
    point x1,y1 to point x2,y2.
    If D is positive, then the point Px,Py is to the LEFT of the
    line segment.  If D is negative, P is to the right of segment.
    If D is zero then, P is on the segment
    If D =0 then that means that the point P is on the line
    defined by the two points...they may not be on the segment

    The check is done by taking the
    cross product of the vectors x1,y1 to x2,y2

    and x1,y1 to Px,Py
    """

    dx = x2 - x1
    dy = y2 - y1
    dxp = Px - x1
    dyp = Py - y1

    return CrossProduct(dx,dy,dxp,dyp)


def CrossProduct(x1,y1,x2,y2):
                # Given vectors x1,y1 and x2,y2
                # this routine returns the cross product
                # which is also the determinant

    return x1*y2 - y1*x2


def test1():    
    ## Now some test code
    ## I think this is obsolete

    min_lat = 100.
    max_lat = 120.
    min_long = 0.
    max_long = 20.

    num_lat = 10
    num_long = 10

    num_LE = 2

    LE1 = np.zeros((num_LE,2),np.float)
    LE2 = np.zeros((num_LE,2),np.float)

    LE1[:,:] = ( ((3.1, 104.9),(20, 100)) )
    LE2[:,:] = ( ((6.8, 113.2),(0,  120)) )

    Hit_Table = np.zeros((num_LE,num_lat*num_long),np.int)


    grid = Grid(min_long, max_long, min_lat, max_lat,num_lat,num_long)


    CR.Grid_hit_test(LE1, LE2, grid, Hit_Table)
    Hit_Table.shape = (2,num_long,num_lat)
    print "after:\n", Hit_Table[0,:,:]
    print "after:\n", Hit_Table[1,:,:]


def test_bin_LE():
    test_file_name = "../Tests/BINARY_LE_TEST.FORCST"
    
    LEs, Beached = read_binLE(test_file_name)
    print LEs
    print Beached
    

if __name__ == "__main__":
    test_bin_LE()    


    
    

    







