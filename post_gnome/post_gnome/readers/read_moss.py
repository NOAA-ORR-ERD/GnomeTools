#!/usr/bin/env python

import os, sys
from __future__ import print_function

class MossPolygon():
    def __init__(self):
 
        self.outerBoundary = [] # a coordinate list
        self.innerBoundaries = [] # a list of coordinate lists
        
    def Print(self):
        print("outerBoundary:",self.outerBoundary)
        print("innerBoundaries:",self.innerBoundaries)



def ReadMossPolygons(mossBaseFileName, printDiagnostic = False):
    '''
    .ms1 file is fixed format
    
    Header lines 56 characters
    char 1-5 : Item Number(NEGATIVE IF THE COORDINATES ARE LON/LAT)
    char 16-45: Attribute Name
    char 51-55: Number of coord. pairs
    
    X,Y, Coordinate pairs 23 characters
    01-11: x coordinate
    12-22: y coordinate
    
    Long,Lat pairs
    char 01-10: LONGITUDE
    char 11-20: LATITUDE
    char 21-22: FLAG 
        0-NORMAL 
        1-INDICATES FIRST POINT OF ISLAND POLYGON
    
    '''
    # import shapely here so it won't be imported if not needed
    from shapely.geometry import Polygon, MultiPolygon

    # what about MAPBOUND, EXTENDEDOUTLOOKTHREAT
    attributesToRead = ["MAPLAND","FORECASTHEAVY","FORECASTMEDIUM","FORECASTLIGHT","FORECASTUNCERTAINTY"]
    
    landBoundariesList = [] 
    heavyBoundariesList = [] 
    mediumBoundariesList = [] 
    lightBoundariesList = []
    uncertaintyBoundariesList = []
     
    extension = ".ms1"
    if os.path.exists(mossBaseFileName + extension):
        inFile = file(mossBaseFileName + extension, 'rU')
        alreadyReadNextLine = False
        while True:  # read the header lines     
            if not alreadyReadNextLine:
                line = inFile.readline()
            alreadyReadNextLine = False
            
            if not line: break # end of this file
            if line.strip() == "" : continue; # blank line
            
            itemNum = int(line[0:5])
            assert itemNum < 0 # we expect long/lat values, so the itemNum should be negative
            
            attributeName = line[15:45].strip()

            numCoordinates = int(line[50:55])
            
            
            # note: for some versions of GNOME analyst
            # the number of coordinates in MAPLAND overflows and is reported as a negative number or incorrect number.
            # To support those files, we will ignore the number and just read lines based on the length of the lines
              
            if printDiagnostic: 
                print(itemNum,attributeName,numCoordinates)
            
            readingOuterBoundary = True
            coordList = []
            
            outerBoundary = None
            innerBoundaries = []
             
            if numCoordinates <= 0:
                print("*** ignoring bad %s header line value: numCoordinates: %d ***" % (attributeName,numCoordinates))
            
            # read the points for the polygon for this header           
            while True: #for i in range(0,numCoordinates):
                # since we don't want to rely on numCoordinates
                # we need to look to see if this is the end of this block of coordinates
                ##################
                line = inFile.readline()
                
                # The lines are fixed format,
                # read until we find a header line
                # header lines lines are longer than coordinate lines 
                if (not line) or len(line.strip()) > 50 : # must be end of file or a header line
                    alreadyReadNextLine = True 
                    break # out of this while loop
                
                if attributeName in attributesToRead:
                    #process this line
                    longitudeStr = line[0:10].strip()
                    latitudeStr = line[10:20].strip()
                    flag = line[20:22].strip()
                    
                    if flag == "1" :
                        #then we are starting a new "inner hole"
                        
                        # enforce having the last point of the polygon equal the first point
                        if len(coordList) > 0 and coordList[0] != coordList[-1]:
                            coordList.append(coordList[0])

                        # save the previous coordList
                        if len(coordList) >= 4: # less than 4 would be a degenerate case                          
                            if readingOuterBoundary:
                                outerBoundary = coordList
                            else:
                                innerBoundaries.append(coordList)
                                
                        # reset the coordinate list 
                        coordList = []                                         
                        readingOuterBoundary = False
                    
                    coordList.append((float(longitudeStr),float(latitudeStr)))
                ############# end of while loop
                
            # finished reading the header
            # record the lists we have filled in
            
            if not attributeName in attributesToRead:
                continue # on to the next header line 
            
            # enforce having the last point of the polygon equal the first point
            if len(coordList) > 0 and coordList[0] != coordList[-1]:
                coordList.append(coordList[0])
                
            #filter out degenerate cases
            if len(coordList) < 4: # less than 4 would be a degenerate case
                print("*** ignoring degenerate polygon ***")
                continue # on to the next header line 
                
            # save the coordinate list                             
            if readingOuterBoundary:
                outerBoundary = coordList
            else:
                innerBoundaries.append(coordList)
                
            
            #save thisPolygon
            if len (outerBoundary) > 0 :
                # outerBoundary,innerBoundaries = VerifyAndFixGnomeAnalystPolygon(attributeName,outerBoundary,innerBoundaries)
                if len (outerBoundary) > 0 :
                    thisPolygon = (outerBoundary,innerBoundaries)
                    if attributeName == "MAPLAND": landBoundariesList.append(thisPolygon)
                    elif attributeName == "FORECASTHEAVY": heavyBoundariesList.append(thisPolygon)
                    elif attributeName == "FORECASTMEDIUM": mediumBoundariesList.append(thisPolygon)
                    elif attributeName == "FORECASTLIGHT": lightBoundariesList.append(thisPolygon)
                    elif attributeName == "FORECASTUNCERTAINTY": uncertaintyBoundariesList.append(thisPolygon) 
                
                
        inFile.close()
        
    # convert the lists of MossPolygons to shapely MultiPolygons
    if len(landBoundariesList) == 0: landPolygons = None
    else: 
        landPolygons = MultiPolygon(landBoundariesList)
        if not landPolygons.is_valid:
            # try analysing and fixing the problem
            landPolygons = DiagnoseAndFixMultiPolygon("MAPLAND",landBoundariesList)
    
    if len(heavyBoundariesList) == 0: heavyPolygons = None
    else: 
        heavyPolygons = MultiPolygon(heavyBoundariesList)
        if not heavyPolygons.is_valid:
            heavyPolygons =DiagnoseAndFixMultiPolygon("FORECASTHEAVY",heavyBoundariesList)
    
    if len(mediumBoundariesList) == 0: mediumPolygons = None
    else: 
        mediumPolygons = MultiPolygon(mediumBoundariesList)
        if not mediumPolygons.is_valid:
            mediumPolygons = DiagnoseAndFixMultiPolygon("FORECASTMEDIUM",mediumBoundariesList)
    
    if len(lightBoundariesList) == 0: lightPolygons = None
    else: 
        lightPolygons = MultiPolygon(lightBoundariesList)
        if not lightPolygons.is_valid:
            lightPolygons = DiagnoseAndFixMultiPolygon("FORECASTLIGHT",lightBoundariesList)
                              
    if len(uncertaintyBoundariesList) == 0: uncertaintyPolygons = None
    else: 
        uncertaintyPolygons = MultiPolygon(uncertaintyBoundariesList)
        if not uncertaintyPolygons.is_valid:
            uncertaintyPolygons = DiagnoseAndFixMultiPolygon("FORECASTUNCERTAINTY",uncertaintyBoundariesList)
        
   
    # clip the oil contours to the shoreline
    # note: we need to check that the polygons are valid before trying to clip to prevent shapely from crashing
    if landPolygons != None :
        if landPolygons.is_valid == False:
            print("*** landPolygons is not valid. We will not clip to the shoreline. ***")
        else :
            if heavyPolygons != None: 
                if heavyPolygons.is_valid == False:
                    print("*** heavyPolygons is not valid. It will not be clipped to the shoreline. ***")
                elif heavyPolygons.intersects(landPolygons):
                    print("clipping heavyPolygons to shoreline")
                    heavyPolygons = heavyPolygons.difference(landPolygons)
                    
            if mediumPolygons != None: 
                if mediumPolygons.is_valid == False:
                    print("*** mediumPolygons is not valid. It will not be clipped to the shoreline. ***")
                elif mediumPolygons.intersects(landPolygons):
                    print("clipping mediumPolygons to shoreline")
                    mediumPolygons = mediumPolygons.difference(landPolygons)
                    
            if lightPolygons != None: 
                if lightPolygons.is_valid == False:
                    print("*** lightPolygons is not valid. It will not be clipped to the shoreline. ***")
                elif lightPolygons.intersects(landPolygons):
                    print("clipping lightPolygons to shoreline")
                    lightPolygons = lightPolygons.difference(landPolygons)
                    
            # note: JerryM wonders we should clip the uncertainty to the shoreline.  That it looks better as simple polygons going over the land.
            if uncertaintyPolygons != None: 
                if uncertaintyPolygons.is_valid == False:
                    print("*** uncertaintyPolygons is not valid. It will not be clipped to the shoreline or oil polygons ***")
                else:
                    print("clipping uncertaintyPolygons to shoreline and oil polygons")
                    
                    for polygons,nameOfPolygons in [(lightPolygons,"lightPolygons"),(mediumPolygons,"mediumPolygons"),(heavyPolygons,"heavyPolygons"),(landPolygons,"landPolygons")]:
                        if polygons != None: 
                            print("taking difference with",nameOfPolygons)
                            newUncertaintyPolygons = uncertaintyPolygons.difference(polygons)
                            print("finished taking difference")
                            if not newUncertaintyPolygons.is_valid:
                                #print "Uncertainty Polygon is no longer valid after taking difference with",polygons,nameOfPolygons
                                s = "*** uncertaintyPolygons have not been clipped to %s ***"%(nameOfPolygons)
                                print(s)
                            else:
                                uncertaintyPolygons = newUncertaintyPolygons
                                
                    
   
           
    return   (landPolygons,heavyPolygons,mediumPolygons,lightPolygons,uncertaintyPolygons)


def VerifyAndFixGnomeAnalystPolygon(attributeName,outerBoundary,innerBoundaries):
    '''
    GNOME Analyst sometimes generated invalid polygons.
    
    One case was with a MAPLAND polygon and seemed to be caused by a clipping routine (perhaps in DOGS) where the outer boundary crossed itself.
    There is nothing we can automatically due for suce a case, so we omit such polygons and inform the user.
    
    Another case had a spurious hole inside another hole n contours output by GNOME Analyst.
    . Rather than figure it out in GNOME Analyst,
    A will instead just leave out any bad holes. If we miss having a hole in the polygon, it is probably not a big deal.
    '''
    # first see if the polygon is valid
    shapelyPolygon = Polygon(outerBoundary,innerBoundaries)
    if shapelyPolygon.is_valid:
        return (outerBoundary,innerBoundaries) # this is a good polygon
    # alright, check to see if the outerBoundary is valid
    shapelyPolygon = Polygon(outerBoundary)
    if not shapelyPolygon.is_valid:
        # the problem is one we can't fix
        print("*** ERROR: bad %s outerBoundary ***" % (attributeName))
        print("Can't fix this problem, so omitting this polygon.")
        if attributeName ==  "MAPLAND": print("** Minor error: The result is that the oil will not be clipped to this island.")
        else: print("*** IMPORTANT ERROR:  an oiled area is being omitted.")
        print(outerBoundary)
        return ([],[]) # can't fix it , so return empty lists, so that this part gets skipped
    # if we get here the outer boundary is OK
    # try eliminating any bad holes
    goodHoles = []
    innerBoundariesToTry = []
    for hole in innerBoundaries:
        shapelyPolygon = Polygon(outerBoundary,goodHoles + [hole])
        if shapelyPolygon.is_valid:
            goodHoles.append(hole)
        else:
            print("*** omitting bad hole in %s ***"%(attributeName))
            print(hole)
                  
    return (outerBoundary,goodHoles)
 
 
def DiagnoseAndFixMultiPolygon(attributeName,boundariesList, printDiagnostic = False):
    print("---")
    print("*** %s is invalid. Running diagnostic... ***"%(attributeName))
    # Note: we have already tested each polygon, so the problem is when we add them to a multipolygon
    # the likely cause is overlapping polygons, so we will automatically union those polygons 
    #printDiagnostic = True
    goodPolygons = []
    badPolygons = []
    useSlowerCode = False
    numPolygons = len(boundariesList)
    if useSlowerCode:
        for i in range(0,numPolygons) :
            if i % 20 == 0:
                print("processing polygon", i)
            thisPolygon = boundariesList[i]
            shapelyMultiPoly = MultiPolygon(goodPolygons + [thisPolygon])
            if shapelyMultiPoly.is_valid:
                goodPolygons.append(thisPolygon)
            else:
                print("*** problem when adding polygon %d ***"%(i))
                print(thisPolygon)
                badPolygons.append(thisPolygon)
    else:
        ''' it can be very slow trying these one at a time, 
        so lets try jumping when we can and dropping back to one at a time when we have a problem
        '''
        numPolygons = len(boundariesList)
        if printDiagnostic: print("numPolygons",numPolygons)
        
        i = 0
        lineReported = 0
        while True:
            if i >= numPolygons:
                if printDiagnostic: print("reached end of list")
                break # end of list
            
            lineToReport = (i/100)*100
            if lineToReport != lineReported:
                lineReported = lineToReport
                if i > 0:
                    print("processing polygon %d of %d"%(i,numPolygons))
                
            continueToTop = False
            for numToJump in [100,50,20,10]:
                if continueToTop: continue
                if i + numToJump > numPolygons:
                    numToJump = numPolygons - i;
                    if printDiagnostic: print("adjusting numToJump to",numToJump)
                #try the jump
                if printDiagnostic: print("trying jump")
                polygonsToAdd = boundariesList[i:i+numToJump]
                shapelyMultiPoly = MultiPolygon(goodPolygons + polygonsToAdd)
                if shapelyMultiPoly.is_valid:
                    goodPolygons = goodPolygons + polygonsToAdd
                    i = i + numToJump
                    if printDiagnostic: print("jumped to %i"%(i))
                    continueToTop = True
                
            if continueToTop: continue
                            
            #resort to one at a time
            if printDiagnostic: print("resorting to one at a time. i = %d"%(i))
            for j in range(0,numToJump):              
                thisPolygon = boundariesList[i]
                shapelyMultiPoly = MultiPolygon(goodPolygons + [thisPolygon])
                if shapelyMultiPoly.is_valid:
                    goodPolygons.append(thisPolygon)
                else:
                    if True: print("polygon %d is bad"%(i))
                    badPolygons.append(thisPolygon)
                i = i+1

    # now automatically generate a "fixed" multipolygon by unioning the bad polygons
    shapelyMultiPoly = MultiPolygon(goodPolygons)
    if len(badPolygons) > 0: 
        print("---")
        print("*** Handling bad polygons ***")
        for poly in badPolygons:
            print("Bad polygon:", poly)
            islandPoly = MultiPolygon([poly])
            if islandPoly.is_valid:
                print("Fixed: The bad polygon was a valid polygon, it has been unioned to the whole.")
                shapelyMultiPoly = shapelyMultiPoly.union(islandPoly)
            else:
                # try to fix this polygon
                '''
                Outer boundary problems:
                Sometimes with clipping of shoreline, there are cases of the outer boundary of the
                shoreline crossing itself.  There is not really anything we can do about such cases.
                '''
                # check to see if the outer boundary is valid
                outerBoundary,innerBoundaries = poly
                
                islandPoly = MultiPolygon([(outerBoundary,[])])
                if not islandPoly.is_valid:
                    print("Outer boundary is not valid.")
                    print("Unable to fix this polygon.")
                    if attributeName == "MAPLAND":
                        print("This is a minor error, it just means the oil contours will not be clipped to this part of the land.")
                    else:
                        print("This is a serious error.  Part of the %s area will be missing."%(attributeName))
                    print("---")
                    continue # we will omit this polygon
                
                print("The outer boundary of the polygon is valid.")
                
                
                ##############################
                '''
                Hole problems:
                There are two kinds of problems that can occur with GNOME Analyst contours.
                
                Sometimes the holes stick slightly out of the outerboundary.
                In such a case we will rely on the fact that the holes were just areas to be removed.
                
                Sometimes there seem to be holes within holes. 
                I'm not sure why GNOME Analyst is doing that, but we will assume that 
                the holes were just areas to be removed.
                '''
                numHoles = len(innerBoundaries)
                if numHoles > 0: print("Examing the %d holes..."%(numHoles))
                assert numHoles > 0 # the only way to get to this part of the code is for a hole to be causing the problem
                numOmittedHoles = 0
                for innerBoundary in innerBoundaries:
                    hole = Polygon(innerBoundary)
                    if not hole.is_valid:
                        numOmittedHoles = numOmittedHoles + 1
                        print("Hole is not valid:", hole)
                        print("Omitting this hole. This is a minor error.")
                    else:
                        # subtract this hole from the islandPolygon
                        islandPoly = islandPoly.difference(hole)
                    
                if numOmittedHoles > 0: 
                    print("Partially fixed: %d invalid holes were not subtracted from this polygon, but this polygon has been unioned to the whole."%(numOmittedHoles))
                elif numHoles > 0: 
                    print("Fixed: all holes successfully subtracted from this polygon and the polygon unioned to the whole.")
                else : 
                    print("Fixed: polygon unioned to the whole.")
                
                shapelyMultiPoly = shapelyMultiPoly.union(islandPoly)
                
            print("---")
        
    return shapelyMultiPoly
     
        
def ReadSpillInfo(mossBaseFileName):
    '''
    0, SPILLID: 
    0, FROM: 
    0, CONTACT: 
    0, ISSUED: 13:20, 7/24/09
    0, VALIDFOR: 06:00, 7/17/09
    0, ADDLEDATA: 
    0, CAVEAT: This trajectory was produced by GNOME (General NOAA Operational Modeling Environment),
    0, CAVEAT: and should be used for educational and planning purposes only--not for a real response. In the
    0, CAVEAT: event of an oil or chemical spill in U.S. waters, contact the U.S. Coast Guard National Response
    0, CAVEAT: Center at 1-800-424-8802.
    0, CAVEAT: 
    0, FROMLOGO: gnome.bmp
    0, LIGHTHEAVYRANGE: 3,5    
    '''
    spillInfo = {} 
    extension = ".ms3"
    if os.path.exists(mossBaseFileName + extension):
        "reading spill info from:", mossBaseFileName + extension
        inFile = open(mossBaseFileName + extension, 'rU')
        while True:       
            line = inFile.readline()
            if not line: break # end of this file
            #print line
            if line.strip() == "" : continue; # blank line
            if line[0:2] != "0," :# lines that start with "0," are the header lines
                break; # we are past the header lines, we don't need to read anymore of this file
            else:
                line = line[2:].strip()
                key = line.split()[0]
                value = line[len(key):].strip()
                # remove the : char from the keyword (should be the last char
                assert key[-1] == ':'
                key = key[:-1] # remove the last char
                if key in spillInfo: 
                    if spillInfo[key] == "":
                        spillInfo[key] = value
                    else: # concatenate the strings
                        spillInfo[key] = spillInfo[key] + " " + value
                else:
                    spillInfo[key] = value
                    
                #print key,value
                
        inFile.close()
        
        # CAVEAT and FROMLOGO were not in the original NOAA standard, so add them if they are not there
        for key in ["CAVEAT","FROMLOGO","LIGHTHEAVYRANGE"]:
            spillInfo.setdefault(key, "") # setdefault will set the key to the default only if it is not already set
                
    else:
        print("File does not exist:", mossBaseFileName + extension)
    return spillInfo

    
'''
0, SPILLID: 
0, FROM: 
0, CONTACT: 
0, ISSUED: 13:20, 7/24/09
0, VALIDFOR: 06:00, 7/17/09
0, ADDLEDATA: 
0, CAVEAT: This trajectory was produced by GNOME (General NOAA Operational Modeling Environment),
0, CAVEAT: and should be used for educational and planning purposes only--not for a real response. In the
0, CAVEAT: event of an oil or chemical spill in U.S. waters, contact the U.S. Coast Guard National Response
0, CAVEAT: Center at 1-800-424-8802.
0, CAVEAT: 
0, FROMLOGO: noaa.bmp
'''

class LE():
    def __init__(self):
 
        self.itemNum = None
        self.latitudeStr = None
        self.longitudeStr = None
        #
        self.type = None
        self.pollutant = None
        self.depthInMetersStr = None
        self.massInKilogramsStr = None
        self.densityInGramsPerCubicCentimeterStr = None
        self.ageInHoursSinceReleaseStr = None # tech doc said age in seconds, but it has always been age in hours
        self.status = None
        
    def Print(self):
        print(self.itemNum,self.latitudeStr,self.longitudeStr,self.type,self.status)
    def __str__(self): 
        return str( (self.itemNum, self.latitudeStr, self.longitudeStr, self.type,self.status), )


def ReadMossLeFiles(mossBaseFileName):
    """
    Reads an LE moss file, and returns a pair of  lists of LE objects , (blackLEs,redLEs)     
    """

    blackLEs = []
    redLEs = []
    leGroups = [(blackLEs,".ms4",".ms5"),(redLEs,".ms6",".ms7")]
            
    for LEs,coordinateExtension,propertyExtension in leGroups:
        if not os.path.exists(mossBaseFileName + coordinateExtension):
            continue # only read the files if we have the coordinates file
        
        coordinatesFile = file(mossBaseFileName + coordinateExtension, 'rU')
        
        if(os.path.exists(mossBaseFileName + propertyExtension)):
            propertiesFile = file(mossBaseFileName + propertyExtension, 'rU')
        else :
            propertiesFile = None
            
        while True:       
            # read the LE coordinates
            # first line holds the item number (negative) 
            line = coordinatesFile.readline()
            if not line: break # end of this file
            if line.strip() == "" : break; # blank line, end of this file
            assert line[15:23] == "LE POINT" 
            thisLE = LE()        
            thisLE.itemNum = int(line[0:5])
            # the next line holds the coordinates
            line = coordinatesFile.readline()
            thisLE.longitudeStr,thisLE.latitudeStr = line.split()[:2]
             
            # read the LE properties
            if(propertiesFile):
                # Note: we will rely on the LEs being in the same order, but we verify that
                line = propertiesFile.readline()
                properties = line.split(",")
                #print properties
                i = 0
                      
                assert thisLE.itemNum == int(properties[i]) # verify the LEs (itemNums) are in the same order
                i = i + 1
                 
                thisLE.type,i = properties[i].strip(),i+1
                assert thisLE.type in ["ABSOLUTEMASS", "RELATIVEMASS", "ABSOLUTEPROBABILITY", "RELATIVEPROBABILITY"]
               
                thisLE.pollutant,i = properties[i].strip(),i+1
                assert thisLE.pollutant in ["GAS", "JP4", "JP5", "DIESEL",
                    "IFO", "BUNKER", "LIGHTCRUDE", "MEDIUMCRUDE", "HEAVYCRUDE", "LAPIO",
                    "CONSERVATIVE"]
                
                thisLE.depthInMetersStr,i = properties[i].strip(),i+1
                thisLE.massInKilogramsStr,i = properties[i].strip(),i+1
                thisLE.densityInGramsPerCubicCentimeterStr,i = properties[i].strip(),i+1
                thisLE.ageInHoursSinceReleaseStr,i = properties[i].strip(),i+1
                thisLE.status,i = properties[i].strip(),i+1
                assert thisLE.status in ["INWATER", "ONBEACH", "OFFMAP"]
                             
            LEs.append(thisLE)
            #thisLE.Print()
                     
        coordinatesFile.close()
        if propertiesFile != None:
            propertiesFile.close()
    
    return (blackLEs,redLEs)


def ReadMossFiles(inputPath):
    # read each of the files that exist
    # verify that we have a file of the right extension
    if not os.path.exists(inputPath):
        print("File does not exist:",inputPath)
        return
    allowedExtentions = [".MS1",".MS2",".MS3",".MS4",".MS5",".MS6",".MS7"]
    extention = inputPath[-4:].upper() #
    if not extention.upper() in allowedExtentions:
        print("File name must end in ", allowedExtentions)
        return
    
 
def FilterLEsByStatus(LEs, statusVal):
    matchingLEs = []
    
    for le in LEs:
        if le.status == statusVal:
            matchingLEs.append(le)
              
    if len(matchingLEs) == 0:
        matchingLEs = None
    return matchingLEs

#############################
#############################           
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        filePath = sys.argv[1]
        ext = filePath.split('.')[-1]
        if len(ext) == 3 and ext[:2].lower() == 'ms':
            # this looks like it has a moss extension, so strip it.
            filePath = filePath[:-4]
            
        printDiagnostic = True
        ReadMossPolygons(filePath,printDiagnostic)
    else:
        ''' 
        If called with no argument, run the sample files.
        '''
        print("***************")
        print("No argument provided.  Processing the sample files...")
        print("***************")
        os.chdir(os.path.dirname(sys.argv[0]))

        inputBasePath = "samplemossfiles/test1/tests"
        #spillInfo = ReadSpillInfo(inputBasePath)
        #print spillInfo
        #blackLEs,redLEs = ReadMossLeFiles(inputBasePath)
        printDiagnostic = False
        
        #inputBasePath = "samplemossfiles/mcity/Mcity.14.michael"
        inputBasePath = "samplemossfiles/boston test/boston test"
        inputBasePath = "samplemossfiles/bad map/3-23-0600"
    
        landPolygons,heavyPolygons,mediumPolygons,lightPolygons,uncertaintyPolygons = ReadMossPolygons(inputBasePath,printDiagnostic)
        #print heavyPolygons,mediumPolygons,lightPolygons
        print(landPolygons)
        print(heavyPolygons)
        print(mediumPolygons)
        print(lightPolygons)
     
    print("done")