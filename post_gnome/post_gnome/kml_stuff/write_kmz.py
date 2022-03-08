#!/usr/bin/env python
"""
write_kmz.py

module that supports writng kml / kmz for trajectory files
"""
from __future__ import print_function
from mako.template import Template
import zipfile, base64
import os,sys
import kml_templates

def write_kmz(kmz_filename, times, best_guess_positions):
    """
    write a kmz file

    at this point, only "best guess" particles are supported

    :param kmz_filename: filename of resulting kmz file

    :param times: list of datetimes of the time for each time step

    :param best_guess_positions: lat-lon of particle positions.

    """
    basefilename = kmz_filename.rstrip(".kmz")
    kmz_filename = basefilename + ".kmz"

    mytemplate = Template(kml_templates.LE_Template)

    print("rendering: ", kmz_filename)

    for i, dt in enumerate(times):
        bestGuessLEs = best_guess_positions[i]
        description = "A description as a placeholder"
        spillInfo = "Sample Spill info"
        
        fromLogo = "" # turn this off for now

        #LEs_Floating = ReadMoss.FilterLEsByStatus(bestGuessLEs,"INWATER")
        #LEs_Beached = ReadMoss.FilterLEsByStatus(bestGuessLEs,"ONBEACH")
        
        #uncertaintyLEs_Floating = ReadMoss.FilterLEsByStatus(uncertaintyLEs,"INWATER")
        #uncertaintyLEs_Beached = ReadMoss.FilterLEsByStatus(uncertaintyLEs,"ONBEACH")
        
   
        mytemplate = Template(mkl_templates.LE_Template)
    
        # ( Z = zulu, or can have +-hh:mm for conversion to UTC
        #mossTime = spillInfo["VALIDFOR"] #  example: 15:00, 08/17/09
        #gmtOffset = 0 # code goes here
        beginTime = dt.isoformat()
        endTime = None
        #print "mossTime %s\n beginTime %s" %(mossTime,beginTime)
        
    
        altitude = 0.1 # raise the best guess LEs up above the water so that they are always on top of the uncertainty LEs
        LE_Sets = [
                   (LEs_Floating,"%s - %s" %(mossTime,"Floating Splots (Best Guess)"),"#YellowDotIcon",altitude),
                   (LEs_Beached,"%s - %s" %(mossTime,"Beached Splots (Best Guess)"),"#YellowXIcon",altitude)
                   ]
        
        # altitude = 0 # on the water surface
        # UncertaintyLE_Sets = [
        #            (uncertaintyLEs_Floating,"%s - %s" %(mossTime,"Uncertainty Floating Splots"),"#RedDotIcon",altitude),
        #            (uncertaintyLEs_Beached,"%s - %s" %(mossTime,"Uncertainty Beached Splots"),"#RedXIcon",altitude)
        #            ]
        
        # Polygon_Sets = [
        #                 (forecastHeavyPolygons,"%s - %s" %(mossTime,"Forecast Heavy"),"#HeavyStyle",altitude),
        #                (forecastMediumPolygons,"%s - %s" %(mossTime,"Forecast Medium"),"#MediumStyle",altitude),
        #                (forecastLightPolygons,"%s - %s" %(mossTime,"Forecast Light"),"#LightStyle",altitude),
        #                (uncertaintyPolygons,"%s - %s" %(mossTime,"Uncertainty"),"#UncertaintyStyle",altitude)                  
        #               ]
        
        seriesItem = MossSeriesItem(beginTime,endTime,LE_Sets)
        BestGuess_Series.append(seriesItem)
        
        seriesItem = MossSeriesItem(beginTime,endTime,UncertaintyLE_Sets)
        Uncertainty_Series.append(seriesItem)
       
        #seriesItem = MossSeriesItem(beginTime,endTime,Polygon_Sets)
        #Polygon_Series.append(seriesItem)
    
    # set the end times
    for timeSeries in [BestGuess_Series,Uncertainty_Series,Polygon_Series]:
        for i in range(0,len(timeSeries)-1):
            timeSeries[i].endTime = timeSeries[i+1].beginTime
        
    kmlstring = mytemplate.render(BestGuess_Series =BestGuess_Series, 
                                  Uncertainty_Series =Uncertainty_Series, 
                                  Polygon_Series =Polygon_Series, 
                                  name=Name,
                                  description = description,
                                  fromLogo = fromLogo)
                                   
    kmzfile = zipfile.ZipFile(outfilename, 'w', compression=zipfile.ZIP_DEFLATED)
    
    kmzfile.writestr('dot.png', base64.b64decode(kml_templates.DotData))
    kmzfile.writestr('x.png', base64.b64decode(kml_templates.XData))
    kmzfile.writestr('scalebar.png', base64.b64decode(kml_templates.ScaleBarData))
    
    kmzfile.writestr(Name+".kml", kmlstring)
    kmzfile.close() 
        
    if printDiagnostic: print(kmlstring)
    
    print("done")


