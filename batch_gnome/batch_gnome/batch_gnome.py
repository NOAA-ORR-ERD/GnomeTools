# this is a module that includes stuff useful for batch GNOME processing

import os
from datetime import datetime, timedelta


class SpillFromFile:
    def __init__(self,
                 Name,
                 FileName,
                 StartTime=None,
                 Windage=None,
                 ):
        

        self.Name = Name
        self.FileName = FileName
        self.StartTime = StartTime
        self.Windage = Windage
    
    def CommandString(self):
        msg = []
        msg.append('MESSAGE createSpill;TO model')
        msg.append('LeFilePath %s '%self.FileName)
        msg.append('NAME %s'%self.Name)
        if self.Windage is not None:
            msg.append('windageA %f; windageB %f'%self.Windage)
        if self.StartTime is not None:
            msg.append('startRelTime %s'%DT2str(self.StartTime))
        return "; ".join(msg) + ";\n"

    
class PointSpill:
    def __init__(self,
                 Name, 
                 NumLEs,
                 StartTime,
                 StartPosition, # expecting signed floats: (long, lat)
                 EndTime = None,
                 EndPosition = None,
                 Windage=None,
                 ):
        self.Name = Name
        self.NumLEs = NumLEs
        self.StartTime = StartTime
        self.EndTime = EndTime
        self.StartPosition  = StartPosition 
        self.EndPosition  = EndPosition 
        self.Windage = Windage

    def CommandString(self):
        msg = []
        msg.append('MESSAGE createSpill;TO model')
        msg.append('NAME %s'%self.name)
        msg.append('numLEs %i'%self.NumLEs)
        if self.Windage is not None:
            msg.append('windageA %f; windageB %f'%self.Windage)
        msg.append('startRelTime %s'%DT2str(self.StartTime))
        if self.EndTime is not None:
            msg.append('endRelTime %s'%DT2str(self.EndTime))
        msg.append('startRelPos %f, %f'%self.StartPosition)
        if self.EndPosition is not None:
            msg.append('endRelPos %f, %f'%self.EndPosition)
        
        return "; ".join(msg) + ";\n"

class ChangeSpill():
    """
    Holds the things you want to change about an existing spill
    """
    def __init__(self,
                 Name, 
                 NumLEs = None,
                 StartTime = None,
                 StartPosition = None, # expecting signed floats: (long, lat)
                 EndTime = None,
                 EndPosition = None,
                 Windage=None,
                 ):
        self.Name = Name
        self.NumLEs = NumLEs
        self.StartTime = StartTime
        self.EndTime = EndTime
        self.StartPosition  = StartPosition 
        self.EndPosition  = EndPosition 
        self.Windage = Windage

    def CommandString(self):
        msg = []
        msg.append('MESSAGE setField;TO %s'%self.Name)
        if self.NumLEs is not None:
            msg.append('numLEs %i'%self.NumLEs)
        if self.Windage is not None:
            msg.append('windageA %f; windageB %f'%self.Windage)
        if self.StartTime is not None:
            msg.append('startRelTime %s'%DT2str(self.StartTime))
        if self.EndTime is not None:
            msg.append('endRelTime %s'%DT2str(self.EndTime))
        if self.StartPosition is not None:
            msg.append('startRelPos %f, %f'%self.StartPosition)
        if self.EndPosition is not None:
            msg.append('endRelPos %f, %f'%self.EndPosition)
        
        return "; ".join(msg) + ";\n"

class JustRun:
    def __init__(self,
                 Duration,
                 StartTime,
                 OutputPath,
                 ):
        self.StartTime = StartTime     
        self.Duration =  Duration   
        self.OutputPath = OutputPath    

    def CommandString(self):
        msg = []
        msg.append('MESSAGE run; TO model')
        msg.append('runDurationInHrs %i'%round(TD2hours(self.Duration)))
        msg.append('startTime %s'%DT2str(self.StartTime))
        msg.append('outputPath %s'%self.OutputPath)

        return "; ".join(msg) + ";\n"
      
        
        
class TapRun:
    def __init__(self,
                 StartTime,
                 EndTime,
                 StartPosition,
                 OutputFileName,
                 OutputFilePath = "",
                 Windage=None,
                 Note = None):
        self.StartTime = StartTime     
        self.EndTime = EndTime
        self.StartPosition  = StartPosition 
        self.OutputFileName = OutputFileName
        self.OutputFilePath = OutputFilePath
        self.Windage = Windage
        self.Note       = Note
        


    def __repr__(self):
        return ("TapRun Object: Properties are: \n"+
                "\tStartTime     = " + `self.StartTime      `+"\n"+
                "\tEndTime       = " + `self.EndTime      `+"\n"+
                "\tStartPosition = " + `self.StartPosition  `+"\n"+
                "\tOutputFileName= " + `self.OutputFileName `+"\n"+
                "\tOutputFilePath= " + `self.OutputFilePath `+"\n"+
                "\tNote      = " + `self.Note       `+"\n")
        

class CommandFile:
    def __init__(self,
                 CommandFileName = None,
                 SaveFilePath    = "",
                 SaveFileName    = None,
                 ModelTimeStep   = None,
                 OutputTimeStep  = None,
                 RunLength       = None,
                 NumLEs      = 1000,
                 SetTimer    = False,
                 ):

        self.CommandFileName = CommandFileName
        self.SaveFilePath    = SaveFilePath   
        self.SaveFileName    = SaveFileName   
        self.ModelTimeStep   = ModelTimeStep 
        self.OutputTimeStep  = OutputTimeStep
        self.RunLength       = RunLength     
        self.NumLEs          = NumLEs    
        self.SetTimer        = SetTimer
        self.CloseWhenDone   = True # the default is to close GNOME when done
        
        self.Runs = []

    def __repr__(self):
        return ("CommandFile Object: Properties are: \n"+
                "CommandFileName= " + `self.CommandFileName `+"\n"+
                "SaveFilePath   = " + `self.SaveFilePath    `+"\n"+
                "SaveFileName   = " + `self.SaveFileName    `+"\n"+
                "ModelTimeStep  = " + `self.ModelTimeStep   `+"\n"+
                "OutputTimeStep = " + `self.OutputTimeStep  `+"\n"+
                "RunLength      = " + `self.RunLength       `+"\n"+
                "NumLEs     = " + `self.NumLEs      `+"\n"+
                "Runs     	= "	+ `self.Runs       `+"\n")

    def write(self):
        if not (self.CommandFileName):
            ("There must be a Command File name and path")
        if self.ModelTimeStep and self.OutputTimeStep:
            if  self.OutputTimeStep % self.ModelTimeStep:
                raise("OutputTimeStep should be an integer multiple of ModelTimeStep\n"
                       "OutputTimeStep: %g, ModelTimeStep: %g"% (self.OutputTimeStep, self.ModelTimeStep))
        file = open(self.CommandFileName,'wt')
        file.write("[GNOME COMMAND FILE]\n")
        file.write("-- This Command File was written by the BatchGnome.py module --  \n")
        file.write("-- It is set up to use a GNOME save file that had been pre set up\n")
        file.write("--\n")
        file.write("-- this would clear all maps etc in case we are not launching \n")
        file.write("MESSAGE close; TO model; \n")
        file.write("--\n--load the command file\n")
        if (self.SaveFileName):
            file.write("Message open; TO model; PATH "+self.SaveFileName+"\n")
        else:
            file.close()
            raise(" There must be a Save File Specified")
        file.write("--\n")
        if self.SetTimer:
            file.write("MESSAGE STARTTIMER; TO model\n\n")
        file.write("-- Now the actual runs \n")
        if self.Runs:
            for Run in self.Runs:
                file.write("MESSAGE runSpill;TO model;")
                file.write("runDurationInHrs %i; "%(self.RunLength))
                file.write("numLEs %i; "%(self.NumLEs,))
                if Run.Windage is not None:
                    file.write("windageA %5.4f; "%Run.Windage[0])
                    file.write("windageB %5.4f; "%Run.Windage[1])
                if self.ModelTimeStep:
                    file.write("timeStepInMinutes %i; "%(self.ModelTimeStep))
                if self.OutputTimeStep:
                    file.write("outputStepInMinutes %i; "%(self.OutputTimeStep))
                file.write("startRelTime %s; "%Run.StartTime)
                if Run.EndTime:
                    file.write("endRelTime %s; "%Run.EndTime)
                file.write("startRelPos %s; "%Run.StartPosition)
                file.write("netcdfPath %s;\n"%os.path.join(Run.OutputFilePath, Run.OutputFileName) )
        else:
            file.close()
            raise("There must be at least one run included")
        if self.SetTimer:
            file.write("MESSAGE STOPTIMER; TO model\n\n")
        if self.CloseWhenDone:
            file.write("\nMESSAGE quit; TO model;\n")

        file.close()
#        try:
#            import MacOS
#            print "Setting Type and Creator for Mac"
#            MacOS.SetCreatorAndType(self.CommandFileName,'COSM','.SAV')
#        except ImportError:
#            # must not be on a mac
#            pass

    def AddRun(self,run):
        self.Runs.append(run)
        
    def AddRuns(self,runs):
        self.Runs.extend(runs)

    def ClearRuns(self):
        self.Runs = []
            
def ConvertToNW(location):
    """
    takes a location as a signed long, lat string and converts to NSEW string
    """
    long, lat = map(float, location.split(",") )
    if not -180 <= long <= 180:
            raise "Longitude out of range"
    if not -90 <= lat <= 90:
            raise "latitude out of range"
    if long < 0:
            EW = "W"
    else:
            EW = "E"
    if lat > 0:
            NS = "N"
    else:
            NS = "S"
    long = abs(long)
    lat = abs(lat)
    return "%f %s %f %s"%(long, EW, lat, NS)


def str2DT(dt_string):
    """
    convert a datetime expressed in the GNOME format to a datetime object
    """
    day, month, year, hour, minute = [int(i) for i in dt_string.split(',')]
    dt = datetime(year, month, day, hour, minute)
    return dt

def DT2str(dt):
    """
    convert a datetime to a string in the GNOME format.
    """
    dt_string = "%3i, %3i, %5i, %3i, %3i"%(dt.day, dt.month, dt.year, dt.hour, dt.minute)
    return dt_string

def TD2hours(td):
    """
    convert a timedelta to integer hours
    """
    return td.days*24 + td.seconds/3600.0 + td.microseconds/3600000.0













