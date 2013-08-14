#!/usr/bin/env python

"""
Code to read a set of OSSM time files, and generate a set of valid start times.
"""

from datetime import datetime, timedelta
import os, random

random.seed(1) # so that all runs get the same sequence

import BatchGnome

class GapSet:
    def __init__(self, Times, MaxGap):
        """
        Times is a list of datetimes
        It is assumed that they are in order
        MaxGap is the maximum size of gaps in hours
        """
        MaxGap = timedelta(hours=MaxGap)
        self.Start = Times[0]
        self.End = Times[-1]
        Gaps = []
        PrevTime = Times[0]
        #Gapstart = PrevTime
        #GapEnd = PrevTime
        for time in Times[1:]:
            if (time - PrevTime) > MaxGap:
                Gaps.append( [PrevTime, time] )
            PrevTime = time

        # merge contiguous gaps:
        NewGaps = []
        if len(Gaps) > 1:
            while Gaps:
                while len(Gaps) > 1 and Gaps[0][1] == Gaps[1][0]:
                    Gaps[0] = ( (Gaps[0][0],Gaps[1][1]) )
                    del Gaps[1]
                NewGaps.append( Gaps[0] )
                del Gaps[0]

        self.Gaps = NewGaps

    def __getitem__(self,i):
        try:
            return self.Gaps[i]
        except IndexError:
            raise
        

    def MergeGaps(self, Gaps2):
        """
        MergeGaps() merges the gaps from this record with the one input.

        """
        Gaps1  = self.Gaps

        if not Gaps1:
            NewGaps = Gaps2.Gaps
        else:
            NewGaps = []
            for gap2 in Gaps2.Gaps:
                while Gaps1:
                    gap1 = Gaps1[0]
                    #print gap1, gap2
                    if gap2[1] < gap1[0]: # before the gap
                        #print "gap2 before gap1"
                        NewGaps.append(gap2)
                        break
                    elif gap2[0] > gap1[1]: # after the gap
                        #print "gap2 after gap1"
                        NewGaps.append(gap1)
                        del Gaps1[0]
                    else: # They intersect
                        #print "They intersect"
                        NewGaps.append( min(gap1[0],gap2[0]), max(gap1[1],gap2[1]) )
                        del Gaps1[0]
        self.Gaps = NewGaps

    def TimeInGap(self, time, span = timedelta(0)):
        start = time
        end = time + span
        print "start:", start
        print "end:", end
        for gap in self.Gaps:
            if end < gap[0] or  start > gap[1]:
                continue
            else:
                return True
        return False

def GetTimes(filename, type):
    print "loading data from :",filename
    data = open(filename).readlines()
    Times = []
    if type == "Hyd":
        data = data[3:]
    elif type == "Wind":
        data = data[5:]
    else:
        raise "Unsupported File type:"+type
    for line in data:
        day,month,year,hour,min = map(int,line.split(",")[:5])
        if year < 100: # assumes all dates after 2000 are 4 digit
            year += 1900
        Times.append(datetime(year,month,day,hour,min))

    return Times

def FindGaps():
    Gaps = []
    for file in files:
        Times = GetTimes(file[0],file[2])
        Gaps.append(GapSet(Times,file[1]))
    
    # find the overlap of the records:
    # find latest start time:
    Start, End = Gaps[0].Start, Gaps[0].End
    
    for gaps in Gaps[1:]:
        if gaps.Start > Start:
            Start = gaps.Start
        if gaps.End < End:
            End = gaps.End
            
    print " the overlap in the records is:"
    print "Start:", Start
    print "End:  ", End
    
    # Merge the Gaps
    
    AllGaps = Gaps[0]
    for gaps in Gaps[1:]:
        AllGaps.MergeGaps(gaps)
    return Start, End, AllGaps

def FindStarts(Start, End, AllGaps, RunTime):
    # find the start times:
    TotalHours = Timedelta2Hours(End - Start - RunTime)
    StartHours = []
    StartTimes = [[] for i in range(len(setup.Seasons))]
    Done = 0
    while not Done:
        StartHour = random.randint(0,TotalHours/6)*6
        if not StartHour in StartHours:
            StartHours.append(StartHour)
            if len(StartHours) >= TotalHours:
                print "I can't find %i valid starts in the records!"%NumStarts
                raise ValueError()
            StartTime = Start + timedelta(hours = StartHour)
            if not AllGaps.TimeInGap(StartTime,RunTime):
                for i in range(len(setup.Seasons)):
                    if StartTime.month in setup.Seasons[i][1]:
                        StartTimes[i].append(StartTime)
                        if len(StartTimes[i]) >= setup.NumStarts:
                            print "done with", setup.Seasons[i][0]
                            setup.Seasons[i][1] = []
                        StillToDo = 0
                        for season in setup.Seasons:
                            if season[1]:
                                StillToDo = 1
                        if not StillToDo:
                            Done = 1
            else:
                print StartTime, "is in a gap"
    for i in range(len(setup.Seasons)):
        stats = {}
        outfilename = os.path.join(setup.RootDir, setup.Seasons[i][0]+"Starts.txt")
        outfile = open(outfilename, 'w')
        print "Writing:", outfilename
        #print setup.Seasons[i][0]
        for time in StartTimes[i]:
            stats[time.year] = stats.setdefault(time.year, 1) + 1
            outfile.write(BatchGnome.DT2str(time)+"\n")
        outfile.close()
        for year, num in stats.items():
            print year, num

def Timedelta2Hours(delta):
    return int(round(delta.days*24 + delta.seconds/3600.))

class EmptyGapSet:
    def TimeInGap(self, time, span = timedelta(0)):
        return False

class SimpleGapSet(GapSet):
    # same thing, but with a different __init__
    def __init__(self, Gaps):
        self.Gaps = Gaps
        

if __name__ == "__main__":
    
    from TAP_Setup import setup

    RunTime = timedelta(hours = setup.TrajectoryRunLength)
    
    if setup.TimeSeries is None:
        Start, End = setup.DataStartEnd
        if setup.DataGaps:
            FindStarts(Start, End, SimpleGapSet(setup.DataGaps), RunTime)
        else:  
            ## fixme: isthe required, or would an empty list be fine for gaps            
            FindStarts(Start, End, EmptyGapSet(), RunTime)
    else:
        files = setup.TimeSeries
        Start, End, AllGaps = FindGaps()
        FindStarts(Start, End, AllGaps, RunTime)










