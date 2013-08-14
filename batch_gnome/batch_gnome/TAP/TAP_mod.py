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

    Hit_Table is a NumPy array of UnsignedInt8 of size (Num_LEs, Num_sites),
    it hold the values of 0 or 1, depending on whether a given LE has hit a given sight.
    ***Hit_Table is ALTERED by this function!!!***


"""

##from Tap_ext.check_receptors import hit_test, Grid_hit_test
#from Hazmat.TAP.TAP_ext import check_receptors as CR
##hit_test = check_receptors.hit_test
##Grid_hit_test = check_receptors.Grid_hit_test

#from Numeric import *
#from Hazmat.TAP.TAP_ext import NumericExtras
#byteswap = NumericExtras.byteswap
#changetype = NumericExtras.changetype
import sys
import time
#import Geometry

import numpy as N

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
    
    data  = fromstring(file.read(8*4*num_LE),'1') # get all the data as bytes
    data.shape = (num_LE,8*4) # reshape into records
    
    Beached = data[:,4*4:5*4].copy() # extract the Beach code.
    #Beached = fromstring(Beached.tostring(),'f') # reconfigure as C float
    changetype(Beached,'f')
    data = data[:,:4*2].copy() # extract the lat and long
    #data = fromstring(data.tostring(),'f') # reconfigure as C float
    changetype(data,'f')

    # Do we need to byteswap?
    if sys.byteorder == 'little':
        byteswap(Beached)
        byteswap(data)
    Beached = Beached.tolist()
        
    data.shape = (num_LE,2)
    # Switch Lat and Long and Assume western hemisphere
    LEs = zeros([num_LE,2],Float)
    LEs[:,0] = -data[:,1]
    LEs[:,1] = data[:,0]
    
    return LEs,Beached
    
def ReadTrajectory_Old(filename,CheckFlag = 1,HeaderData = {}): # Read a binary trajectory file
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
            print flags
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
        data = fromstring(file.read(NumBytes*NumTimesteps*NumLEs),UnsignedInt8)
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

        #LEs = fromstring(LEs.tostring(),Int32)
        changetype(LEs,Int32)
        
        if HeaderData['Endian'] != sys.byteorder:
            byteswap(LEs)
            
        LEs.shape = (NumTimesteps,NumLEs,2)
        
        # convert to python Floats (C doubles)
        Trajectory = LEs.astype(Float)
        
        # convert to degrees (in place)
        divide(Trajectory,array((1e6),Float),Trajectory)
        
        
        flags.shape = (NumTimesteps,NumLEs)
        return (Trajectory,(NumTimesteps,NumLEs),HeaderData,flags)
    else:
        return (None,(None,None),HeaderData,None)


def ReadTrajectory(filename,CheckFlag = 1,HeaderData = {}): # Read a binary trajectory file
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

        DataType = N.dtype([("long","i4"), ("lat", "i4"), ("flag","b1")])
       
        # Read the file
        #data = fromstring(file.read(NumBytes*NumTimesteps*NumLEs),UnsignedInt8)
        data = N.fromfile(file, DataType)
        file.close()

        try:
            data.shape = (NumTimesteps,NumLEs)
        except ValueError:
            raise TAPError(" The size of the data in the Trajectory file:\n%s\n Does not match the header"%
                           filename )

        # extract LE data:
#        LELat = data['lat']
#        LELong = data['long']
#        LELat.shape = LELong.shape = (NumTimesteps,NumLEs,1)
        
        LEs = N.empty((NumTimesteps, NumLEs, 2), dtype=N.int32)

        print LEs.shape
        LEs[:,:,0] = data['long']
        LEs[:,:,1] = data['lat']

#        LEs = N.c_[(LELong,LELat)]

        # extract the flags
        flags = data['flag'].copy()
        del data

        #LEs = fromstring(LEs.tostring(),Int32)
        #changetype(LEs,Int32)
        
        if HeaderData['Endian'] != sys.byteorder:
            LEs.byteswap(True)
            
        #LEs.shape = (NumTimesteps,NumLEs,2)
        
        # convert to python Floats (C doubles)
        Trajectory = LEs.astype(Float)
        
        # convert to degrees (in place)
        N.divide(Trajectory,N.array((1e6),Float),Trajectory)
        
        
        flags.shape = (NumTimesteps,NumLEs)
        return (Trajectory,(NumTimesteps,NumLEs),HeaderData,flags)
    else:
        return (None,(None,None),HeaderData,None)
        
        
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
    
    Hit_Table is a NumPy array of UnsignedInt8 of size (Num_LEs, Num_sites),
    it hold the values of 0 or 1, depending on whether a given LE has hit a given sight.
    ***Hit_Table is ALTERED by this function!!!***
    
    """
    
    def __init__(self, polyset):
        # polyset is a set of polygons, as returned by read_bna
        self.sites = []
        self.bb = []
        for type,site in polyset:
            self.sites.append(site)
            self.bb.append(array((max(site[:,0]), min(site[:,0]), max(site[:,1]), min(site[:,1])),Float))
    def NumSites(self):
        return len(self.sites)
    def __len__(self):
        return len(self.sites)
        
    def comp_hits(self,LEs,Hit_Table,StartTimestep = 0):
        CR.hit_test(LEs,self.sites,self.bb,Hit_Table,StartTimestep)
        return None
        
        # conversion functions for binary TAP II cubes

        
def transform(data, num_LE = 1000, n = 1.5):
    """
    Convert number of LEs to to 8 bit integer for storage in Cubes               
    This is reversed by: (Cube**n)/(255**n) * num_LE;
    input data is a Numpy array

    """
    data = ceil((data.astype(Float)/num_LE)**(1.0/n)*255)  
    return data.astype(UnsignedInt8)

def VolToCube(data,TotalVolume, n = 1.5):
    """
    Convert actual Volume to 8 bit integer for storage in Cubes               
    This is reversed by: (Cube**n)/(255**n) * num_LE;
    input data is a Numpy array
    """
#    data = ceil((data.astype(Float)/TotalVolume)**(1.0/n)*255)  
    data = ceil((data/TotalVolume)**(1.0/n)*255)  
    return data.astype(UnsignedInt8)
    
def reverse_transform(data, num_LE = 1000, n = 1.5):
    """
    NOTE: this isn't used, but is handy for testing
    also, it isn't exact, because of the use of "ceil" in both cases 
    """
    # convert from 8 bit integer NumPy array to Float array
    data = data.astype(Float)
    data = ceil((data/255)**n*num_LE)
    return data
    
def CompTAPIICube(FileList,OutputTimes,Receptors):
    # read the header of the first one:
    (junk,junk,HeaderData,junk) = ReadTrajectory(FileList[0],2)
    
    NumTimesteps = HeaderData['Number of output steps']
    NumLEs = HeaderData['Number of LEs']
    TimeStep = HeaderData["Output time step in seconds"]

    OutputSteps = array(([0] + OutputTimes),Int)*60*60/int(TimeStep)
    OutputSteps2 = array(OutputTimes,Int)*60*60/int(TimeStep)
    
    # Allocate Cube
    NumSpills,NumSites,NumTimes = len(FileList),len(Receptors),len(OutputTimes)

##    ##Old version    
##    Cube = zeros((NumTimes,NumSites,NumSpills),UnsignedInt8)
##    start = time.time()
##    for SpillNum in range(NumSpills):
##        print "computing spill number %i"%(SpillNum,)
##        HitTable = zeros((NumLEs,NumSites),Int) - 1
##        (Trajectory,(NumTimesteps,NumLEs),HeaderData,flags) = ReadTrajectory(FileList[SpillNum],1,HeaderData)
##        for StepNum in range(len(OutputTimes)):
##            Start,Finish = OutputSteps[StepNum],OutputSteps[StepNum+1]+1
##            Receptors.comp_hits(Trajectory[Start:Finish,:,:],HitTable,Start)
##            Cube[StepNum,:,SpillNum] = transform(sum((HitTable >= 0),0),NumLEs)
##    print "Version 1 took %f seconds"%(time.time()-start)

    Cube = zeros((NumTimes,NumSites,NumSpills),UnsignedInt8)
    start = time.time()
    for SpillNum in range(NumSpills):
        print "computing spill number %i"%(SpillNum,)
        HitTable = zeros((NumLEs,NumSites),Int) - 1
        (Trajectory,(NumTimesteps,NumLEs),HeaderData,flags) = ReadTrajectory(FileList[SpillNum],1,HeaderData)
        Receptors.comp_hits(Trajectory,HitTable)
        for StepNum in range(len(OutputTimes)):
            hitSites = ((HitTable <= OutputSteps2[StepNum]) & (HitTable >= 0))
            Cube[StepNum,:,SpillNum] = transform(sum(hitSites),NumLEs)
##    ## test to see if they are the same:
##    print "Version 2 took %f seconds"%(time.time()-start)
##    if alltrue(Cube.flat == Cube2.flat):
##        print "They are all equal"
##    else:
##        print "They are NOT all equal"
##    return Cube, Cube2
    return Cube

def CompThicknessCube(FileList,OutputTimes,Receptors,Weather = None):
    """

    CompThicknessCube computes the average thickness of the oil over
    each receptor site. It only works for grid receptors.

    If Weather is not None it must be a tuple of A,B values, where A and
    B are the coefficient of a weathering curve:

    Volume = InitialVolume * (A*log( TimeSinceRelease) ) + B)

    (I think, I didn't write this down when I wrote the code!)

    If Weather is None, then there is no change in volume.
        
    """

    min_long, max_long, min_lat, max_lat,num_lat,num_long = Receptors.grid
    dlat = (max_lat-min_lat) / num_lat
    dlong = (max_long-min_long) / num_long

    if Weather:
        A, B = Weather
    ## read the header of the first trajectory file: It is assumed that
    ## the others will match..no check is made to assure this, but it
    ## will crash if anything is very wrong.

    (junk,junk,HeaderData,junk) = ReadTrajectory(FileList[0],2)

    NumTimesteps = HeaderData['Number of output steps']
    NumLEs = HeaderData['Number of LEs']
    TimeStep = HeaderData["Output time step in seconds"]
    TimeStepHours = TimeStep / 3600

    OutputSteps = array([0] + OutputTimes,Int)*60*60/int(TimeStep)
    
    # Allocate the Cube
    NumSpills,NumSites,NumTimes = len(FileList),len(Receptors),len(OutputTimes)
    Cube = zeros((NumTimes,NumSites,NumSpills),UnsignedInt8)

    start = time.time() # just for timing how long it takes to run

    ## Loop through each individual trajectory
    for SpillNum in range(NumSpills):
        print "computing spill number %i"%(SpillNum,)

        ## Load the trajectory data (see function in this file)
        ## Trajectory is a large (NumTimesteps X NumLEs X 2) Numeric array with the positions of all the lEs at all time steps
        ## flags is an array (NumTimesteps X NumLEs) of the one-byte flags associated with the LEs
        (Trajectory,(NumTimesteps,NumLEs),HeaderData,flags) = ReadTrajectory(FileList[SpillNum],1,HeaderData)
        Hits = zeros((NumSites),Int)
        VolTable = zeros((NumSites),Float) # this will store the Maximum volume in each grid box.
        released = logical_not(flags & 1) # make a "released" flag.
        del flags
        ReleasedSteps = zeros((NumLEs,),Int) - 1 #Array to store timestep that LEs were released (-1 if not released)

        ## Step through the Cube output time steps
        for step in xrange( len(OutputSteps) -1 ):
            ## step through the Trajectory time steps between each Cube Timestep
            for t in xrange(OutputSteps[step],OutputSteps[step+1]):
                Vol = zeros((NumSites),Int) # this stores the volume (number of LEs)
                LEs = Trajectory[t,:,:] # the slice out of the Trajectory array for this time step
                if Weather:
                    putmask(ReleasedSteps, (released[t,:] & (ReleasedSteps == -1) ), t ) # set release time step for the ones just released
                    WeatherTime = (t - ReleasedSteps) * TimeStepHours
                    ValidLEs = nonzero(WeatherTime >= 1)
                    Weathered_Vol_1 = A*log( take(WeatherTime,ValidLEs) ) + B # compute just the ones over one hour
                    Weathered_Vol =  zeros((NumLEs,), Float)
                    put( Weathered_Vol, ValidLEs, Weathered_Vol_1)
                    putmask( Weathered_Vol, WeatherTime < 1, 1 )
                    # just to make sure
                    if len(nonzero(Weathered_Vol < 0)) > 0 or len(nonzero(Weathered_Vol > 1)) > 0:
                        raise TAPError("Weathered Volume out of range")
                else:
                    Weathered_Vol = ones((NumLEs,), Float)
                ### What boxes are the LEs in?
                i = ( (LEs[:,0] - min_long) / dlong ).astype(Int)
                j = ( (LEs[:,1] - min_lat ) / dlat ).astype(Int)
                # check for whether the LEs are in the bounds of the grid at all
                InGrid =  nonzero( (i < num_long) &
                                    (i >= 0) &
                                    (j < num_lat) &
                                    (j >= 0) )
                ## the indexes of the LEs that are within the grid
                ## translated form (i,j) to a one-d array
                site_ind = (take(i,InGrid) * num_lat ) + take(j,InGrid)
                # The Volumes of the LEs that are in the grid
                LE_Vol = take(Weathered_Vol, InGrid)
                ## add the volume of each LE to the total volume in it's grid box
                ## this is where you would change the volume of the LE according to your weathering curve 
                for k in xrange(len(site_ind)):
                	Vol[site_ind[k]] += LE_Vol[k]
                ## put the Volume in the grid box only if it's larger than the volume there
                VolTable = maximum(Vol, VolTable)
            ## put the max volume in the Cube at this Cube time step
            Cube[step,:,SpillNum] = transform(VolTable,NumLEs)
    return Cube

class GridReceptor_set:
    """
    Constructor:
    Receptor_set(min_long, max_long, min_lat, max_lat,num_lat,num_long):
    
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
    
    Hit_Table is a NumPy array of type Int of size (Num_LEs, Num_sites),
    it hold the values of 0 or 1, depending on whether a given LE has hit a given sight.
    ***Hit_Table is ALTERED by this function!!!***
    
    """
    
    def __init__(self, min_long, max_long, min_lat, max_lat,num_lat,num_long):
        self.grid = (min_long, max_long, min_lat, max_lat,num_lat,num_long)

##    def BuildPolygonSet(self):
##        #import Geometry
##        (min_long, max_long, min_lat, max_lat,num_lat,num_long) = self.grid
##        dlat = float(max_lat-min_lat) / num_lat
##        dlong = float(max_long-min_long) / num_long
        
##        points = zeros((5*self.num_lat*self.num_long, 2),Float)
##        indices = zeros((num_lat*num_long, 2),Int)
        
##        polynum = 0
##        for i in range(self.num_lat):
##            for j in range(self.num_long):
##                indices[polynum,:] = (polynum*5, polynum*5+5)
##                points[polynum*5,:] = [min_long + dlong*i , min_lat + dlat*j]
##                points[polynum*5+1,:] = [min_long + dlong*(i+1) , min_lat + dlat*j]
##                points[polynum*5+2,:] = [min_long + dlong*(i+1) , min_lat + dlat*(j+1)]
##                points[polynum*5+3,:] = [min_long + dlong*i , min_lat + dlat*(j+1)]
##                points[polynum*5+4,:] = [min_long + dlong*i , min_lat + dlat*j]
##                polynum += 1
##        return Geometry.PolygonSet((points,indices))                

    def WriteBNA(self,file):
        print "Writing the BNA data"
        (min_long, max_long, min_lat, max_lat,num_lat,num_long) = self.grid
        dlat = float(max_lat-min_lat) / num_lat
        dlong = float(max_long-min_long) / num_long
        print num_lat, num_long
        polynum = 1
        for i in range(num_long):
             for j in range(num_lat):
                 print "writing polygon # %i"%polynum
                 file.write('"#%i","1",5\n'%polynum)
                 file.write('%11.6f, %11.6f\n'%(min_long + dlong*i , min_lat + dlat*j))     
                 file.write('%11.6f, %11.6f\n'%(min_long + dlong*(i+1) , min_lat + dlat*j))
                 file.write('%11.6f, %11.6f\n'%(min_long + dlong*(i+1) , min_lat + dlat*(j+1)))
                 file.write('%11.6f, %11.6f\n'%(min_long + dlong*i , min_lat + dlat*(j+1)))
                 file.write('%11.6f, %11.6f\n'%(min_long + dlong*i , min_lat + dlat*j))    
                 polynum += 1
        return True


    def NumSites(self):
        return self.grid[4] * self.grid[5]
    def __len__(self):
        return self.NumSites()
        
    def comp_hits(self,LEs,Hit_Table,StartTimestep = 0):
        CR.Grid_hit_test(LEs, self.grid, Hit_Table, StartTimestep)
        return None


##class Grid:
##    def __init__(self, min_long, max_long, min_lat, max_lat,num_lat,num_long):
##        self.min_long = min_long
##        self.max_long = max_long
##        self.min_lat  = min_lat 
##        self.max_lat  = max_lat 
##        self.num_lat  = num_lat    
##        self.num_long = num_long   
##        self.dlat = float(max_lat-min_lat) / num_lat
##        self.dlong = float(max_long-min_long) / num_long

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


if __name__ == "__main__":
    
    ## Now some test code
    ## I think this is obsolete

    min_lat = 100.
    max_lat = 120.
    min_long = 0.
    max_long = 20.

    num_lat = 10
    num_long = 10

    num_LE = 2

    LE1 = zeros((num_LE,2),Float)
    LE2 = zeros((num_LE,2),Float)

    LE1[:,:] = ( ((3.1, 104.9),(20, 100)) )
    LE2[:,:] = ( ((6.8, 113.2),(0,  120)) )

    Hit_Table = zeros((num_LE,num_lat*num_long),Int)


    grid = Grid(min_long, max_long, min_lat, max_lat,num_lat,num_long)


    CR.Grid_hit_test(LE1, LE2, grid, Hit_Table)
    Hit_Table.shape = (2,num_long,num_lat)
    print "after:\n", Hit_Table[0,:,:]
    print "after:\n", Hit_Table[1,:,:]

    
    

    







