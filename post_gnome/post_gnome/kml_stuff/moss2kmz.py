#!/usr/bin/env python

"""
A simple prototype of building a kml/kmz file for Google earth from
GNOME moss files

In this version the template and pngs are embedded.

This needs additions like drawing LEs on the beach differently, not drawing
LEs that are off teh mpa, or evaporated, or ...

"""
from __future__ import print_function
from mako.template import Template
import zipfile, base64


def ReadLEMossFile(filename):
    """
    Reads an LE moss file,a nd returns the coordinates of the LEs
    
    NOTE: this may not be very robust!
    """
    mossfile = file(filename, 'rU')
    LEs = []
    while True:
        mossfile.readline()
        line = mossfile.readline()
        if not line: break
        LEs.append(line.split()[:2])
    return LEs

def buildkmz(mossfilename):
    BestGuess =  ReadLEMossFile(mossfilename+".ms4")
    Uncertainty = ReadLEMossFile(mossfilename+".ms6")
    Name = mossfilename

    mytemplate = Template(LE_Template)

    outfilename = mossfilename+".kmz"
    print("rendering: ", outfilename)

    kmlstring = mytemplate.render(LEs=BestGuess, LEs_U=Uncertainty, IconName="DotIcon", name=Name)

#    kmzfile = zipfile.ZipFile(outfilename, 'w')
    kmzfile = zipfile.ZipFile(outfilename, 'w', compression=zipfile.ZIP_DEFLATED)
    kmzfile.writestr('RedDot.png', base64.b64decode(RedDotData))
    kmzfile.writestr('BlackDot.png', base64.b64decode(BlackDotData))
    kmzfile.writestr('YellowDot.png', base64.b64decode(BlackDotData))
    kmzfile.writestr(Name+".kml", kmlstring)
    kmzfile.close()       


### the template:
LE_Template = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>${name}</name>

    <Style id="RedDotIcon">
      <IconStyle>
         <scale>1</scale>
         <Icon>
            <href>RedDot.png</href>
         </Icon>
         <hotSpot x="0.5"  y="0.5" xunits="fraction" yunits="fraction"/>  
      </IconStyle>
      <LabelStyle>
         <color>00000000</color>
      </LabelStyle>
    </Style>

    <Style id="BlackDotIcon">
      <IconStyle>
         <scale>1</scale>
         <Icon>
            <href>BlackDot.png</href>
         </Icon>
         <hotSpot x="0.5"  y="0.5" xunits="fraction" yunits="fraction"/>  
      </IconStyle>
      <LabelStyle>
         <color>00000000</color>
      </LabelStyle>
    </Style>

    <Style id="YellowDotIcon">
      <IconStyle>
         <scale>1</scale>
         <Icon>
            <href>YellowDot.png</href>
         </Icon>
         <hotSpot x="0.5"  y="0.5" xunits="fraction" yunits="fraction"/>  
      </IconStyle>
      <LabelStyle>
         <color>00000000</color>
      </LabelStyle>
    </Style>

    <StyleMap id="BlackDotPair">
        <Pair>
            <key>normal</key>
                <styleUrl>#BlackDotIcon</styleUrl>
        </Pair>
        <Pair>
            <key>highlight</key>
                <styleUrl>#YellowDotIcon</styleUrl>
        </Pair>
    </StyleMap>

    <StyleMap id="RedDotPair">
        <Pair>
            <key>normal</key>
            <styleUrl>#RedDotIcon</styleUrl>
        </Pair>
        <Pair>
            <key>highlight</key>
            <styleUrl>#RedDotIcon</styleUrl>
        </Pair>
    </StyleMap>

    <Placemark>
      <name>Uncertainty</name>
      <styleUrl>#RedDotPair</styleUrl>
      <MultiGeometry>
        % for x, y in LEs_U:
                 <Point>
                   <coordinates>${x},${y}</coordinates>
                 </Point>
        % endfor
      </MultiGeometry>
    </Placemark>

    <Placemark>
      <name>Best Guess</name>
      <styleUrl>#BlackDotPair</styleUrl>
      <MultiGeometry>
        % for x, y in LEs:
                 <Point>
                   <coordinates>${x},${y}</coordinates>
                 </Point>
        % endfor
      </MultiGeometry>
    </Placemark>

 
  </Document>
</kml>
"""

# These icons (these are base64 encoded 3-pixel sized dots in a 32x32 transparent PNG)    
BlackDotData  = "iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAABYgAAAWIBXyfQUwAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAAABaSURBVFiF7dWxDQAhDEPR79NNd+x765mGCSykNI6UgoLwhIvINpP1jL5eQAEFFFBAAcCbXpQk4DvH3+latR01sACfXumc8QiU/tytCGLArRqPoIACCiiggAI2js1pdbjIwXwAAAAASUVORK5CYII="
RedDotData    = "iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAABZAAAAWQB3ySUyAAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAAABaSURBVFiF7daxCcAwDETRr5Dp4n2z3rnxBMKg5gtUuLD84Aq5kjBZz+jrAgQIECBAAPC2b1YV8J3TT3etJuk1rEBOr+6c8Qiq/SG5FEEfcKnGIxAgQIAAAQI2p2Rrcf2y0MYAAAAASUVORK5CYII="
YellowDotData = "iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAABZAAAAWQB3ySUyAAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAAACESURBVFiF7ZUxCsMwEATnDPpO8in9KC87v8fSpkiauLKFiRHswjVXrAY0cCGJO7Pc+roBDGAAAxhgYoAotKj0SHokLSpEGaqSdH42qhr6mY060hVD57hHAo/ddmXR82zVpA6I16Hdoa4RB0T5epBq5Of/Kf9z4MJM6oABDGAAAxjgwrwBd2LRMaa/WmgAAAAASUVORK5CYII="

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        name = sys.argv[1]
        ext = name.split('.')[-1]
        if len(ext) == 3 and ext[:2] == 'ms':
            # this looks like it has a moss extesion, so strip it.
            name = name[:-4]
    else:
        name = "tests"
        
    buildkmz(name)
    

