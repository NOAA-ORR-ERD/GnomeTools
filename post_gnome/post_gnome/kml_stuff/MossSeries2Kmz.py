#!/usr/bin/env python

"""
This builds a kml/kmz file for Google earth from
a series of GNOME moss files

In this version the template and pngs are embedded.

"""

from mako.template import Template
import zipfile, base64
import os,sys
import ReadMoss
import Moss2Kmz
import MossKmlTemplates


def MossTime2GoogleTime(mossTime, gmtOffset):
        # mossTime = spillInfo["VALIDFOR"] # 15:00, 08/17/09
        #print mossTime
        timeStr,dateStr = mossTime.split(",")
        dateStr = dateStr.strip()
        monthStr,dayStr,yearStr = dateStr.split("/")
        year = int(yearStr)
        if gmtOffset == 0 :
            offsetStr = "Z"
        elif gmtOffset  > 0 :
            offsetStr = "+02d:00"%(gmtOffset)
        else :
            offsetStr = "-02d:00"%(gmtOffset)

        if year < 50 :
            fourDigitYearStr = "20"+yearStr
        elif year <= 99 :
            fourDigitYearStr = "19"+yearStr
        else:  # already a 4 digit year
            fourDigitYearStr = yearStr

        googleStr = "%s-%s-%sT%s:00%s" % (fourDigitYearStr,monthStr,dayStr,timeStr,offsetStr)
        #print googleStr
        return googleStr


class MossSeriesItem():
    def __init__(self,beginTime,endTime,setsOfObjects):

        self.beginTime = beginTime
        self.endTime = endTime
        self.setsOfObjects = setsOfObjects




def buildkmzSeries(mossSeriesFilename,printDiagnostic = False):
    Name = os.path.basename(mossSeriesFilename) # mossfilename may include the path, so get just the file name

    assert Name.endswith("000") #require that the last three chars of the name be "000"
    baseName = Name[:-3] # chop off the last three characters
    basePathAndName = mossSeriesFilename[:-3]

    seriesDescription = GetDescription(mossSeriesFilename)

    outfilename = basePathAndName+".kmz"
    print "rendering: %s"%outfilename

    BestGuess_Series = []
    Uncertainty_Series = []
    Polygon_Series = []

    for i in range(0,999) :
        mossfilename =  "%s%03d"%(basePathAndName,i)
        print "Processing",mossfilename

        extension = ".ms3"
        if not os.path.exists(mossfilename + extension):
            break # no more moss files in the series

        bestGuessLEs,uncertaintyLEs = ReadMoss.ReadMossLeFiles(mossfilename)
        description = Moss2Kmz.GetDescription(mossfilename)
        spillInfo = ReadMoss.ReadSpillInfo(mossfilename)

        fromLogo = spillInfo["FROMLOGO"] # usually either a gnome or noaa logo
        #if os.path.exists(os.path.join(os.path.dirname(mossSeriesFilename) + "scalebar.png")):
        #fromLogo = "scalebar.png" # try this
        fromLogo = "" # turn this off for now

        LEs_Floating = ReadMoss.FilterLEsByStatus(bestGuessLEs,"INWATER")
        LEs_Beached = ReadMoss.FilterLEsByStatus(bestGuessLEs,"ONBEACH")

        uncertaintyLEs_Floating = ReadMoss.FilterLEsByStatus(uncertaintyLEs,"INWATER")
        uncertaintyLEs_Beached = ReadMoss.FilterLEsByStatus(uncertaintyLEs,"ONBEACH")

        landPolygons,forecastHeavyPolygons,forecastMediumPolygons,forecastLightPolygons,uncertaintyPolygons = ReadMoss.ReadMossPolygons(mossfilename)

        #print description

        mytemplate = Template(LE_Template)

        # ( Z = zulu, or can have +-hh:mm for conversion to UTC
        mossTime = spillInfo["VALIDFOR"] #  example: 15:00, 08/17/09
        gmtOffset = 0 # code goes here
        beginTime = MossTime2GoogleTime(mossTime, gmtOffset)
        endTime = None
        #print "mossTime %s\n beginTime %s" %(mossTime,beginTime)


        altitude = 0.1 # raise the best guess LEs up above the water so that they are always on top of the uncertainty LEs
        LE_Sets = [
                   (LEs_Floating,"%s - %s" %(mossTime,"Floating Splots (Best Guess)"),"#YellowDotIcon",altitude),
                   (LEs_Beached,"%s - %s" %(mossTime,"Beached Splots (Best Guess)"),"#YellowXIcon",altitude)
                   ]

        altitude = 0 # on the water surface
        UncertaintyLE_Sets = [
                   (uncertaintyLEs_Floating,"%s - %s" %(mossTime,"Uncertainty Floating Splots"),"#RedDotIcon",altitude),
                   (uncertaintyLEs_Beached,"%s - %s" %(mossTime,"Uncertainty Beached Splots"),"#RedXIcon",altitude)
                   ]

        Polygon_Sets = [
                        (forecastHeavyPolygons,"%s - %s" %(mossTime,"Forecast Heavy"),"#HeavyStyle",altitude),
                       (forecastMediumPolygons,"%s - %s" %(mossTime,"Forecast Medium"),"#MediumStyle",altitude),
                       (forecastLightPolygons,"%s - %s" %(mossTime,"Forecast Light"),"#LightStyle",altitude),
                       (uncertaintyPolygons,"%s - %s" %(mossTime,"Uncertainty"),"#UncertaintyStyle",altitude)
                      ]

        seriesItem = MossSeriesItem(beginTime,endTime,LE_Sets)
        BestGuess_Series.append(seriesItem)

        seriesItem = MossSeriesItem(beginTime,endTime,UncertaintyLE_Sets)
        Uncertainty_Series.append(seriesItem)

        seriesItem = MossSeriesItem(beginTime,endTime,Polygon_Sets)
        Polygon_Series.append(seriesItem)




    if (
        len(BestGuess_Series) == 0
        and len(Uncertainty_Series) == 0
        and len(Polygon_Series) == 0
        ):
        print "no files processed"
        return

    # set the end times
    for timeSeries in [BestGuess_Series,Uncertainty_Series,Polygon_Series]:
        for i in range(0,len(timeSeries)-1):
            timeSeries[i].endTime = timeSeries[i+1].beginTime


    kmlstring = mytemplate.render(
                                  BestGuess_Series =BestGuess_Series,
                                  Uncertainty_Series =Uncertainty_Series,
                                  Polygon_Series =Polygon_Series,
                                  name=Name,
                                  description = description,
                                  fromLogo = fromLogo)

    kmzfile = zipfile.ZipFile(outfilename, 'w', compression=zipfile.ZIP_DEFLATED)

    kmzfile.writestr('dot.png', base64.b64decode(MossKmlTemplates.DotData))
    kmzfile.writestr('x.png', base64.b64decode(MossKmlTemplates.XData))
    kmzfile.writestr('scalebar.png', base64.b64decode(MossKmlTemplates.ScaleBarData))

    kmzfile.writestr(Name+".kml", kmlstring)
    kmzfile.close()

    if printDiagnostic: print kmlstring

    print "done"


def GetDescription(mossfilename):
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

    spillInfo = ReadMoss.ReadSpillInfo(mossfilename)
    if spillInfo == None:  return None


    str = ""

    for label,key in [("Valid for","VALIDFOR"),
                      ("From","FROM"),
                      ("Contact","CONTACT"),
                       ("Issued","ISSUED")]:
        if spillInfo[key]: str = str+"<b>%s:</b> %s<br>\n"%(label,spillInfo[key])

    if spillInfo["CAVEAT"]:
        str = str + "<br>\n"+ spillInfo["CAVEAT"]

    return str


#############################################
### the templates:



################################
################################

LE_Template = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>${name}</name>
    <open>1</open>
      % if description :
      <description><![CDATA[${description}]]></description>
      % endif

""" + (
   MossKmlTemplates.LeStyles +
   MossKmlTemplates.ContourStyles
)  +  """

% for seriesType, series in [("Best Guess",BestGuess_Series),("Uncertainty",Uncertainty_Series)]:
% if series and len(series) > 0 :
 <Folder>
 <name>${seriesType}</name>
% for seriesItem in series:
% for LEs,LEsname,LEsstyle,altitude in seriesItem.setsOfObjects:
   % if LEs :
     <Placemark>
      <name>${LEsname}</name>
      <styleUrl>${LEsstyle}</styleUrl>
      %if seriesItem.beginTime or seriesItem.endTime:
      <TimeSpan id="ID">
          % if seriesItem.beginTime :
          <begin>${seriesItem.beginTime}</begin>     <!-- kml:dateTime -->
          % endif
          % if seriesItem.endTime :
          <end>${seriesItem.endTime}</end>         <!-- kml:dateTime -->
          % endif
      </TimeSpan>
      % endif
      <MultiGeometry>
        % for le in LEs:
             <Point>
                 % if altitude == 0 :
                     <coordinates>${le.longitudeStr},${le.latitudeStr}</coordinates>
                 % else :
                     <altitudeMode>relativeToGround</altitudeMode>
                     <coordinates>${le.longitudeStr},${le.latitudeStr},${altitude}</coordinates>
                 % endif
             </Point>
        % endfor

      </MultiGeometry>
     </Placemark>
    %endif
% endfor
% endfor
 </Folder>
% endif
% endfor

""" + (
   MossKmlTemplates.LogoOverlay
)  +  """

  </Document>
</kml>
"""



if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        filePath = sys.argv[1]
        ext = filePath.split('.')[-1]
        if len(ext) == 3 and ext[:2].lower() == 'ms':
            # this looks like it has a moss extension, so strip it.
            filePath = filePath[:-4]
        buildkmzSeries(filePath)
    else:
        '''
        If called with no argument, run the sample files.
        '''
        print "***************"
        print "No argument provided.  Processing the sample files..."
        print "***************"
        os.chdir(os.path.abspath(os.path.dirname(sys.argv[0])))
        for filePath in ["samplemossfiles/jerryseries/jerryseries000"]:
            buildkmzSeries(filePath)


