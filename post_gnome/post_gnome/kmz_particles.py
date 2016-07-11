#!/usr/bin/env python

"""
Code for workign with particle fiels in mkz

Only handles reading for now
"""
# for py2/3 compatibility
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import zipfile
import base64
from datetime import datetime

import numpy as np

from post_gnome import nc_particles

file_attributes = nc_particles.file_attributes
var_attributes = nc_particles.var_attributes


class Writer(object):
    """
    class to write kmz files sutable for Google Earth
    """
    time_formatter = '%m/%d/%Y %H:%M'

    def __init__(self,
                 filename,
                 num_timesteps=None,
                 ref_time=None,
                 file_attributes=file_attributes,
                 var_attributes=var_attributes
                 ):

        """
        create a kmz_particle file Writer

        Creates the kml / kmz file, opens it for writing,
        writes the global attributes, creates required variables etc.

        :param filename: name of kmz file to open - if it exists,
                         it will be written over!

        :param num_timesteps=None: number of timesteps that will be output. Must be defined for netcdf3.
                                   Can be None for netcdf4
        :type num_timesteps: integer

        :param ref_time=None: reference time for time units (i.e. seconds since..).
                              If None, the first time used will be used.
        :type ref_time: datetime object

        :param file_attributes: keys and values for teh file-level attributes.
                                Defaults to the set defined in this module.
        :type file_attributes: dict

        :param var_attributes: dist of variable names, and the keys and values for variable
                               attributes.
                               Defaults to the set defined in this module.
        :type var_attributes: dict

        :param nc_version=3: version of netcdf to use -- must be 3 or 4. If 4, some extra
                             features are enabled.
        :type nc_version: integer
        """

        # strip off the .kml or .kmz
        filename = filename.rstrip(".kml").rstrip(".kmz")

        self.filename = filename + ".kmz"
        self.kml_name = os.path.split(filename)[-1] + ".kml"

        self.num_timesteps = num_timesteps
        self.ref_time = ref_time

        self.file_attributes = file_attributes
        self.var_attributes = var_attributes

        # create a list to hold what will be the contents of the kml
        self.kml = [header_template.format(caveat=caveat,
                                           kml_name=self.kml_name,
                                           # fixme: this is only doing now!
                                           valid_timestring=datetime.now().strftime(self.time_formatter),
                                           issued_timestring=datetime.now().strftime(self.time_formatter),
                                           )]
        # # fixme: put real data in here from nc_file?
        # self.kml = [header_template.format(caveat=caveat,
        #                                    kml_name=self.kml_name,
        #                                    valid_timestring=model_start_time.strftime(self.time_formatter),
        #                                    issued_timestring=datetime.now().strftime(self.time_formatter),
        #                                    )]

        # # Global attributes
        # # put some of this in the kml file?
        # for (name, value) in self.file_attributes.items():
        #     setattr(nc, name, value)
        self.closed = False





    def write_timestep(self, timestamp, timestep, data, uncertain=False):
        """
        write the data for a timestep

        :param timestamp: the time stamp of the timestep
        :type timestamp: datetime object

        :param timestep: the timestep between this and the next data point
        :type timestep: datetime object

        :param data: dict of data arrays -- all parameters for a single time step
        :type data: dict

        :param uncertain=False: Is this an uncertaintly run?
        :type uncertain: bool
        """

        start_time = timestamp.isoformat()
        end_time = (timestamp + timestep).isoformat()


        positions = np.c_[data['longitude'], data['latitude']]

        try:
            in_water = 2
            on_land = 3

            water_positions = positions[data['status_codes'] == in_water]
            beached_positions = positions[data['status_codes'] == on_land]
        except KeyError:
            water_positions = positions
            beached_positions = np.zeros((0, 2), dtype=np.float64)

        self.kml.append(build_one_timestep(water_positions,
                                           beached_positions,
                                           start_time,
                                           end_time,
                                           uncertain
                                           ))

    def close(self):
        """
        close the kmz file

        This forces the write of the file
        """
        if not self.closed:
            self.kml.append(footer)
            with zipfile.ZipFile(self.filename, 'w', compression=zipfile.ZIP_DEFLATED) as kmzfile:
                kmzfile.writestr('dot.png', base64.b64decode(DOT))
                kmzfile.writestr('x.png', base64.b64decode(X))
                # write the kml file
                kmzfile.writestr(self.kml_name, "".join(self.kml).encode('utf8'))
            self.closed = True
            return True
        else:
            return False

    def __del__(self):
        """ make sure to close the netcdf file """
        # anything to close?
        self.close()



class Reader(object):
    """
    Class to handle reading a nc_particle file

    (such as those written by GNOME or the Writer class above)
    """
    def __init__(self, kml_file):
        """
        initialize a file reader.

        :param kml_file: the kml/kmz file to read.
        :type kml_file: string

        """
        raise NotImplementedError

    @property
    def variables(self):
        """
        return the names of all the variables associated with the particles
        """
        raise NotImplementedError

    def __str__(self):
        return ("kml_particles Reader object:\n"
                "variables: {}\n"
                "number of timesteps: {}\n"
                ).format(self.variables, len(self.times))

    def get_all_timesteps(self, variables=['latitude', 'longitude']):
        """
        returns the requested variables data from all timesteps as a
        dictionary keyed by the variable names

        :param variables: the variables desired as a list string names.
                          Defaults to ['latitude','longitude']
        :type variables: list of strings

        :returns data: returns a dict of arrays -- the keys are the
                       variable names, and the values are numpy arrays
                        of the data. The arrays are the flattened ragged
                        array of data.
        """
        raise NotImplementedError

    def get_units(self, variable):
        """
        return the units of the given variable

        :param variable: name of the variable for which the units are required
        :type variable: string
        """
        raise NotImplementedError

    def get_attributes(self, variable):
        """
        return all the attributes of the given variable

        :param variable: name of the variable for which the attributes are required
        :type variable: string
        """
        raise NotImplementedError

    def get_timestep(self, timestep, variables=['latitude', 'longitude']):
        """
        returns the requested variables data from a given timestep as a
        dictionary keyed by the variable names

        :param variables: The variables desired as a list string names.
                          Defaults to ['latitude','longitude']
        :type variables: list of strings

        :returns data: returns a dict of arrays -- the keys are the
                       variable names, and the values are numpy arrays
                       of the data.
        """
        raise NotImplementedError

    def get_individual_trajectory(self, particle_id, variables=['latitude', 'longitude']):
        """
        returns the requested variables from trajectory of an individual particle

        note: this is inefficient -- it has to read the entire file to get it.
        """
        raise NotImplementedError

    def close(self):
        """
        close the kml file

        -- anything to be done?
        """
        pass


    def __del__(self):
        """ make sure to close the file """
        self.close()


def nc2kmz(nc_file, kmz_file=None):
    """
    convert a nc_particles file to kmz

    :param nc_file: name of nertcdf file to read

    :param kmz_file=None: name of kmz file to write. If None, the nc_file's name wil be used, with .kmz as teh extansion.

    """

    if kmz_file is None:
        root = nc_file
        root = root[:-3] if root.endswith(".nc") else root
        kmz_file = root + ".kmz"

    reader = nc_particles.Reader(nc_file)

    # create a kmz writer:
    writer = Writer(kmz_file)
    variables = reader.variables
    # loop to read / write the data
    for step, time in enumerate(reader.times):
        try:
            timestep = reader.times[step + 1] - time
        except IndexError:
            timestep = reader.times[-1] - reader.times[-2]
        # get the data
        data = reader.get_timestep(step, variables)
        writer.write_timestep(time, timestep, data, uncertain=False)
    writer.close()

    return kmz_file


# Templates for the kmz files

caveat = ("This trajectory was produced by GNOME (General NOAA Operational Modeling",
          " Environment), and should be used for educational and planning purposes only",
          "--not for a real response. In the event of an oil or chemical spill in U.S.",
          "waters, contact the U.S. Coast Guard National Response Center at 1-800-424-8802."
          )

# The kml templates:
header_template = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>{kml_name}</name>
    <open>1</open>
    <description><![CDATA[<b>Valid for:</b> {valid_timestring}<br>
                          <b>Issued:</b>{issued_timestring} <br>
                          {caveat}]]>
    </description>

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
         <color>ff00ffff</color>
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
         <color>ff00ffff</color>
         <hotSpot x="0.5"  y="0.5" xunits="fraction" yunits="fraction"/>
      </IconStyle>
      <LabelStyle>
         <color>00000000</color>
      </LabelStyle>
    </Style>
"""


point_template = """             <Point>
                     <altitudeMode>relativeToGround</altitudeMode>
                     <coordinates>{:.6f},{:.6f},1.000000</coordinates>
             </Point>
"""


timestep_header_template = """<Folder>
  <name>{date_string}:{certain}</name>
"""

one_run_header = """    <Placemark>
      <name>{certain} {status} Splots </name>
      <styleUrl>{style}</styleUrl>
      <TimeSpan id="ID">
        <begin>{start_time}</begin>     <!-- kml:dateTime -->
        <end>{end_time}</end>         <!-- kml:dateTime -->
      </TimeSpan>
      <MultiGeometry>
"""
one_run_footer = """      </MultiGeometry>
    </Placemark>
"""
timestep_footer = """
</Folder>
"""


def build_one_timestep(floating_positions,
                       beached_positions,
                       start_time,
                       end_time,
                       uncertain,
                       ):

    data = {'certain': "Uncertainty" if uncertain else "Best Guess",
            'start_time': start_time,
            'end_time': end_time,
            'date_string': start_time,
            }
    kml = []
    kml.append(timestep_header_template.format(**data))

    for status, positions in [('Floating', floating_positions),
                              ('Beached', beached_positions)]:
        color = "Red" if uncertain else "Yellow"
        data['style'] = "#" + color + "DotIcon" if status == "Floating" else "#" + color + "XIcon"

        data['status'] = status
        kml.append(one_run_header.format(**data))

        for point in positions:
            kml.append(point_template.format(*point[:2]))
        kml.append(one_run_footer)
    kml.append(timestep_footer)

    return "".join(kml)

footer = """
  </Document>
</kml>
"""
# These icons (these are base64 encoded 3-pixel sized dots in a 32x32 transparent PNG)
#   these were encoded by the "build_icons" script
DOT = "iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAAJOgAACToB8GSSSgAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAAAEASURBVFiF7ZY7DsIwEEQfET09Ej11lFtwK06Re3ANlCoFPQpnoGJoHClCXpOPg10wUhonnnlyvF5vJJFSRdL0P0AOANsZcwqgAkrg6MZuQANcgdckN0ljn52kWlInW537ZjfWd2z4SVIbCP5U6+ZEAThLek4I7/V0cxcBnGaGDyGCK/Htn09ZdkutAnsiBFBHCO9VWzkb+XtBAdyB/Ywy9ekBHPCUqHUQVRHDcV6V74UFUEYMD3paAEdjfIm8nsl7gQVwWyHL62kBNCsAeD2zLcMXcIkUjvPyt+nASZj8KE7ejLJox1lcSIZ7IvqVzCrDkKJeSucARFW2veAP8DO9AXV74Qmb/4vgAAAAAElFTkSuQmCC"
X = "iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAAN1wAADdcBQiibeAAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAAAHKSURBVFiFrdXLq01hGMfx12HMzMCUU4zFQEYiROYkEkkpHbeTXI5LSDqHtomBEJGY+RMMGBlKKWVmaiDXzvExsN7admuv9azLU89k7ef5fb/ruhMSVuIy3uEOVhXH++w1mMEbnMFSpITl+Ob/+oOpHuFHiszh+oIVCbPGVx8Sh0vguaYT3lcIdJU4VAGHtwm3agTaShysgcMgYUNAoKnEgQAcVueFqR4l9mMhkHVJ8RbkPt6DxL4g/EreGQ3oIrE3CL86vFd2FidaSOzBfGDn+ihv3KU82UBidxB+o4xV9TBFJSKX/eY4Tt0TfSooUVWzVYzIO326A3yuLj/6YWkjcTuSHRVImG4AH0RzJ1K8PqSUFoKzn8KpQdNd+N3wFoT+OyLwnfjVEB6WqIPv6AAPSVTBt+NnR3itxDj4tiD8Hs52kSiDb8WPQOB9LCp2WkuMwrcE4Q8xMbJ7ro3EcMBmfA8EPCqBt5bIi5uC8McV8Nznm0gkLMPXwMKTADz3haDExoRjgcGnWByEN5EYJLyuGXrWAp57pib7Y8K1ioHnHeC5L1bkP0iYHPPjCyzpCK+SmMdkHliLl8XBVzjaIzz3Ov++H59xF+uR/gJmOo2+fdNArAAAAABJRU5ErkJggg=="

