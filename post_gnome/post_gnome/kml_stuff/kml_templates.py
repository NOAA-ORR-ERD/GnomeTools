"""
kml_templates

various templates for writting GNOME KML files
"""


LeStyles = """
    <Style id="RedDotIcon">
      <IconStyle>
         <scale>0.2</scale>
         <color>ff0000ff</color>
         <Icon>
            <href>dot.png</href>
         </Icon>
          <hotSpot x="0.5"  y="0.5" xunits="fraction" yunits="fraction"/>  
      </IconStyle>
      <LabelStyle>
         <color>00000000</color>
      </LabelStyle>
    </Style>

    <Style id="BlackDotIcon">
      <IconStyle>
         <scale>0.2</scale>
         <Icon>
            <href>dot.png</href>
         </Icon>
         <color>ff000000</color>
         <hotSpot x="0.5"  y="0.5" xunits="fraction" yunits="fraction"/>  
      </IconStyle>
      <LabelStyle>
         <color>00000000</color>
      </LabelStyle>
    </Style>
    
    <Style id="YellowDotIcon">
      <IconStyle>
         <scale>0.2</scale>
         <Icon>
            <href>dot.png</href>
         </Icon>
         <color>ff00efef</color>
         <hotSpot x="0.5"  y="0.5" xunits="fraction" yunits="fraction"/>  
      </IconStyle>
      <LabelStyle>
         <color>00000000</color>
      </LabelStyle>
    </Style>
 
    <Style id="RedXIcon">
      <IconStyle>
         <scale>0.2</scale>
         <color>ff0000ff</color>
         <Icon>
            <href>x.png</href>
         </Icon>
          <hotSpot x="0.5"  y="0.5" xunits="fraction" yunits="fraction"/>  
      </IconStyle>
      <LabelStyle>
         <color>00000000</color>
      </LabelStyle>
    </Style>

    <Style id="BlackXIcon">
      <IconStyle>
         <scale>0.2</scale>
         <Icon>
            <href>x.png</href>
         </Icon>
         <color>ff000000</color>
         <hotSpot x="0.5"  y="0.5" xunits="fraction" yunits="fraction"/>  
      </IconStyle>
      <LabelStyle>
         <color>00000000</color>
      </LabelStyle>
    </Style>
    
    <Style id="YellowXIcon">
      <IconStyle>
         <scale>0.2</scale>
         <Icon>
            <href>x.png</href>
         </Icon>
         <color>ff00efef</color>
         <hotSpot x="0.5"  y="0.5" xunits="fraction" yunits="fraction"/>  
      </IconStyle>
      <LabelStyle>
         <color>00000000</color>
      </LabelStyle>
    </Style>

"""


ContourStyles = """
    
    <Style id="HeavyStyle">
        <LineStyle><color>ff9E0000</color><width>1</width></LineStyle>  
        <PolyStyle><fill>1</fill><color>e09E0000</color></PolyStyle>
    </Style>
    
    <Style id="MediumStyle">
        <LineStyle><color>ffff7326</color><width>1</width></LineStyle>  
        <PolyStyle><fill>1</fill><color>e0ff7326</color></PolyStyle>
    </Style>
    
    <Style id="LightStyle">
        <LineStyle><color>ffffff00</color><width>1</width></LineStyle>  
        <PolyStyle><fill>1</fill><color>e0ffff00</color></PolyStyle>
    </Style>
    
    <Style id="UncertaintyStyle">
        <LineStyle><color>e0dcdcdc</color><width>2</width></LineStyle>  
        <PolyStyle><fill>1</fill><color>e0dcdcdc</color></PolyStyle>
    </Style>
"""

SaveOrigUncertainty = """
    <Style id="UncertaintyStyle">
        <LineStyle><color>ffff00ff</color><width>2</width></LineStyle>  
        <PolyStyle><fill>0</fill><color>e0ff00ff</color></PolyStyle>
    </Style>

"""

LogoOverlay = """
    % if fromLogo != "" :
    <Folder><name>Screen Overlays</name>
    <ScreenOverlay> 
    <name>Logo</name> 
    <color>ffffffff</color> 
    <visibility>1</visibility> 
    <Icon> 
    <href>${fromLogo}</href> 
    </Icon> 
    <overlayXY x="1" y="1" xunits="fraction" yunits="fraction"/> 
    <screenXY x="1" y="1" xunits="fraction" yunits="fraction"/> 
    <size x="-1" y="-1" xunits="fraction" yunits="fraction"/> 
    </ScreenOverlay> 
    </Folder>
    % endif
    """
    
    
# These icons (these are base64 encoded 3-pixel sized dots in a 32x32 transparent PNG)    
DotData="iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAAJOgAACToB8GSSSgAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAAAEASURBVFiF7ZY7DsIwEEQfET09Ej11lFtwK06Re3ANlCoFPQpnoGJoHClCXpOPg10wUhonnnlyvF5vJJFSRdL0P0AOANsZcwqgAkrg6MZuQANcgdckN0ljn52kWlInW537ZjfWd2z4SVIbCP5U6+ZEAThLek4I7/V0cxcBnGaGDyGCK/Htn09ZdkutAnsiBFBHCO9VWzkb+XtBAdyB/Ywy9ekBHPCUqHUQVRHDcV6V74UFUEYMD3paAEdjfIm8nsl7gQVwWyHL62kBNCsAeD2zLcMXcIkUjvPyt+nASZj8KE7ejLJox1lcSIZ7IvqVzCrDkKJeSucARFW2veAP8DO9AXV74Qmb/4vgAAAAAElFTkSuQmCC"
XData="iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAAN1wAADdcBQiibeAAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAAAHKSURBVFiFrdXLq01hGMfx12HMzMCUU4zFQEYiROYkEkkpHbeTXI5LSDqHtomBEJGY+RMMGBlKKWVmaiDXzvExsN7admuv9azLU89k7ef5fb/ruhMSVuIy3uEOVhXH++w1mMEbnMFSpITl+Ob/+oOpHuFHiszh+oIVCbPGVx8Sh0vguaYT3lcIdJU4VAGHtwm3agTaShysgcMgYUNAoKnEgQAcVueFqR4l9mMhkHVJ8RbkPt6DxL4g/EreGQ3oIrE3CL86vFd2FidaSOzBfGDn+ihv3KU82UBidxB+o4xV9TBFJSKX/eY4Tt0TfSooUVWzVYzIO326A3yuLj/6YWkjcTuSHRVImG4AH0RzJ1K8PqSUFoKzn8KpQdNd+N3wFoT+OyLwnfjVEB6WqIPv6AAPSVTBt+NnR3itxDj4tiD8Hs52kSiDb8WPQOB9LCp2WkuMwrcE4Q8xMbJ7ro3EcMBmfA8EPCqBt5bIi5uC8McV8Nznm0gkLMPXwMKTADz3haDExoRjgcGnWByEN5EYJLyuGXrWAp57pib7Y8K1ioHnHeC5L1bkP0iYHPPjCyzpCK+SmMdkHliLl8XBVzjaIzz3Ov++H59xF+uR/gJmOo2+fdNArAAAAABJRU5ErkJggg=="

################################
# main Template
################################

LE_Template = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>${name}</name>
    <open>1</open>
      % if description :
      <description><![CDATA[${description}]]></description>
      % endif
      
""" + LeStyles + """
    
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
   LogoOverlay
)  +  """  
    
  </Document>
</kml>
"""



