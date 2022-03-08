from __future__ import print_function
from six import iteritems
import json
import os
import matplotlib.colors as colors

cnames = colors.cnames
# Add the single letter colors.
for (name, rgb) in iteritems(colors.ColorConverter.colors):
    hex_ = colors.rgb2hex(rgb)
    cnames[name] = hex_

sites = {}
sites['dwh'] = {'site_id': 1, 'zoom': 6, 'longitude': -90.42, 'latitude': 28.03}
sites['arctic'] = {'site_id': 1, 'zoom': 5, 'longitude': -150.00, 'latitude': 70.00}
sites['atlantic'] = {'site_id': 2, 'zoom': 7, 'longitude': -74.10, 'latitude': 	38.94}
sites['caribbean'] = {'site_id': 3, 'zoom': 6, 'longitude': -65.44, 'latitude': 18.20}
sites['greatlakes'] = {'site_id': 9, 'zoom': 1, 'longitude': -83.44, 'latitude': 45.60}
sites['gulfofmexico'] = {'site_id': 4, 'zoom': 4, 'longitude': -90.42, 'latitude': 28.03}
sites['northwest'] = {'site_id': 6, 'zoom': 6, 'longitude': -125.00, 'latitude': 47.00}
sites['pacific'] = {'site_id': 7, 'zoom': 4, 'longitude': -178.00, 'latitude': 11.00}
sites['southwest'] = {'site_id': 8, 'zoom': 6, 'longitude': -122.00, 'latitude': 37.00}
    
            
def particles(package_dir, fn, params):
    '''
    params is a dict with the following keys:
        site_name
        event (can be None)
        title
        shape_zipfilename
        meta_data
        attachment_file - just set up to do one at the moment
        color - matplotlib named colors (e.g. "r" or "red")
        folder_path - a list with the hierarchy of subfolders in TOC
    '''
    
    filename = os.path.join(package_dir, 'layers', fn + '.json')
    f = open(filename, 'w')
    
    if params['classitem'] == 'status':
    
        layer_classes = [
          {
            "styles": [
              {
                "angle_field": None,
                "angle": None,
                "font_character": None,
                "map_layer_class": 101964,
                "ordering": 1,
                "symbol": {
                  "style_text": None,
                  "symbol_def": "NAME 'cross'\nTYPE vector\nPOINTS\n 2 0\n 2 4\n -99 -99\n 0 2\n 4 2\nEND\n",
                  "symbol_desc": "Cross",
                  "symbol_name": "cross"
                },
                "size_field": None,
                "color": cnames[params['color_beached']],
                "style_width": None,
                "offset_x": None,
                "offset_y": None,
                "outlinesymbol": None,
                "outlinecolor": None,
                "width_field": None,
                "style_size": 5
              }
            ],
            "name": "Beached",
            "ordering": 2,
            "expression_type": {
              "expression_type": "V",
              "description": "Exact Value"
            },
            "labels": [
              
            ],
            "expression": "3"
          },
          {
            "styles": [
              {
                "angle_field": None,
                "angle": None,
                "font_character": None,
                "map_layer_class": 101965,
                "ordering": 1,
                "symbol": {
                  "style_text": None,
                  "symbol_def": "NAME 'filledcircle'\nTYPE ellipse\nFILLED true\nPOINTS\n 1 1\nEND\n",
                  "symbol_desc": "Circle (Filled)",
                  "symbol_name": "filledcircle"
                },
                "size_field": None,
                "color": cnames[params['color']],
                "style_width": None,
                "offset_x": None,
                "offset_y": None,
                "outlinesymbol": None,
                "outlinecolor": None,
                "width_field": None,
                "style_size": 5
              }
            ],
            "name": "Floating",
            "ordering": 3,
            "expression_type": {
              "expression_type": "V",
              "description": "Exact Value"
            },
            "labels": [

            ],
            "expression": "2"
          }
        ]
        
    elif params['classitem'] == 'surf_conc':
    
        layer_classes = [
          {
            "styles": [
              {
                "angle_field": None,
                "angle": None,
                "font_character": None,
                "map_layer_class": 169684,
                "ordering": 1,
                "symbol": {
                  "symbol_def": "NAME 'filledcircle'\nTYPE ellipse\nFILLED true\nPOINTS\n 1 1\nEND\n", 
                  "symbol_name": "filledcircle", 
                  "style_text": None, 
                  "symbol_desc": "Circle (Filled)"
                }, 
                "size_field": None,
                "color": "#00ffff",
                "style_width": None,
                "offset_x": None,
                "offset_y": None,
                "outlinesymbol": None,
                "outlinecolor": "#00ffff",
                "width_field": None,
                "style_size": 5
              }
            ],
            "name": "Light",
            "ordering": 1,
            "expression_type": {
              "expression_type": "M",
              "description": "Mapserver Expressions"
            },
            "labels": [],
            "expression": "[surf_conc] >= 0 AND [surf_conc] < 10"
          },
          {
            "styles": [
              {
                "angle_field": None,
                "angle": None,
                "font_character": None,
                "map_layer_class": 169683,
                "ordering": 1,
                "symbol": {
                  "symbol_def": "NAME 'filledcircle'\nTYPE ellipse\nFILLED true\nPOINTS\n 1 1\nEND\n", 
                  "symbol_name": "filledcircle", 
                  "style_text": None, 
                  "symbol_desc": "Circle (Filled)"
                }, 
                "size_field": None,
                "color": "#007fff",
                "style_width": None,
                "offset_x": None,
                "offset_y": None,
                "outlinesymbol": None,
                "outlinecolor": "#007fff",
                "width_field": None,
                "style_size": 5
              }
            ],
            "name": "Medium",
            "ordering": 2,
            "expression_type": {
              "expression_type": "M",
              "description": "Mapserver Expressions"
            },
            "labels": [],
            "expression": "[surf_conc] >= 10 AND [surf_conc] < 50"
          },
          {
            "styles": [
              {
                "angle_field": None,
                "angle": None,
                "font_character": None,
                "map_layer_class": 169682,
                "ordering": 1,
                "symbol": {
                  "symbol_def": "NAME 'filledcircle'\nTYPE ellipse\nFILLED true\nPOINTS\n 1 1\nEND\n", 
                  "symbol_name": "filledcircle", 
                  "style_text": None, 
                  "symbol_desc": "Circle (Filled)"
                }, 
                "size_field": None,
                "color": "#0000ff",
                "style_width": None,
                "offset_x": None,
                "offset_y": None,
                "outlinesymbol": None,
                "outlinecolor": "#0000ff",
                "width_field": None,
                "style_size": 5
              }
            ],
            "name": "Heavy",
            "ordering": 3,
            "expression_type": {
              "expression_type": "M",
              "description": "Mapserver Expressions"
            },
            "labels": [],
            "expression": "[surf_conc] >= 50"
          }
        ]
    
        
    else:
    
        print('no classitem specified for styling')

    layer_info = {
      "layer_type": "wms internal",
      "folder_path": ' > '.join(params['folder_path']),
      "additional_metainfo": params['metadata'],
      "ephemeral": False,
      "sensitivity": "Responders",
      "background_layer": False,
      "legend_annotation": None,
      "applicability": [
        {
          "event": params['event'],
          "site": {
            "site_path": "/" + params['site_name'],
            "site_name": params['site_name'],
            "site_id": sites[params['site_name']]['site_id'],
            "zoom": sites[params['site_name']]['zoom'],
            "longitude": sites[params['site_name']]['longitude'],
            "extent": None,
            "latitude": sites[params['site_name']]['latitude'],
            "baselayer": 16
          }
        }
      ],
      "tilecache": False,
      "modified_by": None,
      "layer_names": None,
      "refresh_rate": 0,
      "title": params['title'],  # name of layer as it appears in TOC
      "single_tile": True,
      "created_by": None,
      "legendgraphic": None,
      "marker_icon": None,
      "external_proj": None,
      "strokecolor": None,
      "opacity": 1,
      "mapfile_layer": {
        "layer_type": "point",
        "opacity": 100,
        "sort_field": "gid",
        "classitem": params['classitem'],
        "labelitem": None,
        "date_modified": None,
        "shapefile": {
          "file": "file://source_files/" + params['shape_zipfilename'],
          "has_hotlinks": False,
          "name": params['shape_zipfilename'].split('.')[0], #a "friendly" display name for the shapefile
          "description": None,
          "created_by": None,
          "keywords": None,
          "allow_download": True,
          "timezone_fields": "{\"time\": \"US/Pacific\"}",
          "content_type": "application/zip",
          "srid": 4326
        },
        "time_column": "time",
        "maxscaledenom": None,
        "created_by": None,
        "labelminscaledenom": None,
        "layer_classes": layer_classes,
        "sort_field": "surf_conc",
        "sort_order": None,
        "labelmaxscaledenom": None,
        "template": False,
        "modified_by": None,
        "layer_desc": 'Styling for trajectory points file',
        "minscaledenom": None,
        "layer_name": "trajectory_points_layer"
      },
      "metadata_url": None,
      "visibility": None,
      "strokewidth": None,
      "proxy": False,
      "animate_script": None,
      "legend": None,
      "uuid": None,
      "gutter": None,
      "graphicname": None,
      "date_modified": None,
      "url": None,
      "metadata_path": None,
      "pointradius": None,
      "fillcolor": None
    }

    if params['attachment_file'] is not None:
        layer_info["attachments"] = [
            {"ordering": 1,
             "attachment": "file://attachments/" + params['attachment_file'],
             "content_type": "application/gif"}
        ]

    json.dump(layer_info, f, indent=2)
    f.close()


def contours(package_dir, fn, params):
    '''
    params is a dict with the following keys:
        site_name
        event (can be None)
        title
        shape_zipfilename
        meta_data
        attachment_file - just set up to do one at the moment
        color - matplotlib named colors (e.g. "r" or "red")
        folder_path - a list with the hierarchy of subfolders in TOC
        uncertain - just has to have this param for uncertainty layer
    '''

    filename = os.path.join(package_dir,'layers',fn + '.json')
    f = open(filename,'w')

    if 'uncertain' in params:

        layer_classes = [
              {
                "styles": [
                  {
                    "angle_field": None,
                    "angle": None,
                    "font_character": None,
                    "map_layer_class": 102054,
                    "ordering": 1,
                    "symbol": None,
                    "size_field": None,
                    "color": None,
                    "style_width": None,
                    "offset_x": None,
                    "offset_y": None,
                    "outlinesymbol": None,
                    "outlinecolor": "#bf0000",
                    "width_field": None,
                    "style_size": None
                  }
                ],
                "name": "Uncertainty Boundary",
                "ordering": 1,
                "expression_type": {
                  "expression_type": "V",
                  "description": "Exact Value"
                },
                "labels": [

                ],
                "expression": "Uncertainty"
              }
            ]

    elif 'SinglePoly' in params:
    

        layer_classes = [
              {
                "styles": [
                  {
                    "angle_field": None,
                    "angle": None,
                    "font_character": None,
                    "map_layer_class": 102031,
                    "ordering": 1,
                    "symbol": None,
                    "size_field": None,
                    "color": "#00ffff",
                    "style_width": None,
                    "offset_x": None,
                    "offset_y": None,
                    "outlinesymbol": None,
                    "outlinecolor": "#000000",
                    "width_field": None,
                    "style_size": None
                  }
                ],
                "name": "Most Probable Location",
                "ordering": 1,
                "expression_type": {
                  "expression_type": "V",
                  "description": "Exact Value"
                },
                "labels": [
                ],
                "expression": "Light"
              }
            ]
    else:

        layer_classes = [
           {
                "styles": [
                  {
                    "angle_field": None,
                    "angle": None,
                    "font_character": None,
                    "map_layer_class": 102030,
                    "ordering": 1,
                    "symbol": None,
                    "size_field": None,
                    "color": "#0000ff",
                    "style_width": None,
                    "offset_x": None,
                    "offset_y": None,
                    "outlinesymbol": None,
                    "outlinecolor": "#000000",
                    "width_field": None,
                    "style_size": None
                  }
                ],
                "name": "Heavy",
                "ordering": 4,
                "expression_type": {
                      "expression_type": "V",
                      "description": "Exact Value"
                },
                "labels": [
                ],
                "expression": "Heavy"
              },
              {
                "styles": [
                  {
                    "angle_field": None,
                    "angle": None,
                    "font_character": None,
                    "map_layer_class": 102029,
                    "ordering": 1,
                    "symbol": None,
                    "size_field": None,
                    "color": "#007fff",
                    "style_width": None,
                    "offset_x": None,
                    "offset_y": None,
                    "outlinesymbol": None,
                    "outlinecolor": "#000000",
                    "width_field": None,
                    "style_size": None
                  }
                ],
                "name": "Medium",
                "ordering": 3,
                "expression_type": {
                    "expression_type": "V",
                    "description": "Exact Value"
                },
                "labels": [
                ],
                "expression": "Medium"
              },
              {
                "styles": [
                  {
                    "angle_field": None,
                    "angle": None,
                    "font_character": None,
                    "map_layer_class": 102031,
                    "ordering": 1,
                    "symbol": None,
                    "size_field": None,
                    "color": "#00ffff",
                    "style_width": None,
                    "offset_x": None,
                    "offset_y": None,
                    "outlinesymbol": None,
                    "outlinecolor": "#000000",
                    "width_field": None,
                    "style_size": None
                  }
                ],
                "name": "Light",
                "ordering": 2,
                "expression_type": {
                  "expression_type": "V",
                  "description": "Exact Value"
                },
                "labels": [
                ],
                "expression": "Light"
              }
              ]
              
    layer_info = {
        "layer_type": "wms internal",
        "folder_path": ' > '.join(params['folder_path']),
        "additional_metainfo": params['metadata'],
        "ephemeral": False,
        "sensitivity": "Responders",
        "background_layer": False,
        "legend_annotation": None,
        "applicability": [
          {
            "event": params['event'],
            "site": {
                "site_path": "/" + params['site_name'],
                "site_name": params['site_name'],
                "site_id": sites[params['site_name']]['site_id'],
                "zoom": sites[params['site_name']]['zoom'],
                "longitude": sites[params['site_name']]['longitude'],
                "extent": None,
                "latitude": sites[params['site_name']]['latitude'],
                "baselayer": 16
            }
         }
      ],
      "tilecache": True,
      "modified_by": None,
      "layer_names": None,
      "refresh_rate": 0,
      "title": params['title'], #name of layer as it appears in TOC
      "single_tile": False,
      "created_by": None,
      "legendgraphic": None,
      "marker_icon": None,
      "external_proj": None,
      "strokecolor": None,
      "opacity": 1,
      "mapfile_layer": {
        "layer_type": "polygon",
        "opacity": 100,
        "sort_field": "type",
        "classitem": "type",
        "labelitem": None,
        "date_modified": None,
        "shapefile": {
          "file": "file://source_files/" + params['shape_zipfilename'],
          "has_hotlinks": False,
          "name": params['shape_zipfilename'].split('.')[0], #a "friendly" display name for the shapefile
          "description": None,
          "created_by": None,
          "keywords": None,
          "allow_download": True,
          "timezone_fields": None,
          "content_type": "application/zip",
          "srid": 4326
        },
        "time_column": None,
        "maxscaledenom": None,
        "created_by": None,
        "labelminscaledenom": None,
        "layer_classes": layer_classes,
        "sort_order": None,
        "labelmaxscaledenom": None,
        "template": False,
        "modified_by": None,
        "layer_desc": 'Styling for trajectory contours',
        "minscaledenom": None,
        "layer_name": "trajectory_contour_layer"
      },
      "metadata_url": None,
      "visibility": None,
      "strokewidth": None,
      "proxy": False,
      "animate_script": None,
      "legend": None,
      "uuid": None,
      "gutter": None,
      "graphicname": None,
      "date_modified": None,
      "url": None,
      "metadata_path": None,
      "pointradius": None,
      "fillcolor": None
    }

    if params['attachment_file'] is not None:
        layer_info["attachments"] = [
            {"ordering": 1,
             "attachment": "file://attachments/" + params['attachment_file'],
             "content_type": "application/gif"}
        ]

    json.dump(layer_info, f, indent=2)
    f.close()