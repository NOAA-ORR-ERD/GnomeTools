import numpy as np
import shapely.geometry as sgeo
import fiona
   
def format_for_bna(polys,levs,dl=0):
    data = []
    ii = 1
    for lev,poly in zip(levs,polys):
        x = np.array(poly[0])
        y = np.array(poly[1])
        if dl:
            x = (x < 0).choose(x,x+360)
        if lev == 0:
            data.append('"Map Bounds", "2", ' + str(len(x)))
        else:
            data.append('"' + str(ii) + '","' + str(int(lev)) + '",' + str(len(x)))
        for lon,lat in zip(x,y):
            data.append("{0:.6f}".format(lon) + ', ' + "{0:.6f}".format(lat))
        ii = ii+1
    data = '\n'.join(data)
    return data

def load_polys_and_clip(shpfile,wl,el,nl,sl,dl=0,lev=1):
    
    '''
    Loop through polygons in shapefile, if they overlap with user
    specified rectangle, load the polygon and clip 
    Note that GSHHS polygons are -180-->180, shorelines that cross the 180 line are
    clipped into west and east polygons, so there is some special casing to handle
    clip rectangles that cross the dateline
    Change lev to 2 if loading lakes
    These *may not* be able to be used to clip trajectories from GAn in Arc?
    '''
    
    bound_box=(wl,sl,el,nl) 
    clip_box = sgeo.box(*bound_box)
    
    if dl:
        bound_box=(-180,sl,180,nl) 

    bna_polys = []
    
    with fiona.open(shpfile) as source:
        for s in source.filter(bbox=bound_box):
            geo = s['geometry']['coordinates']
            for c in geo:
                poly = sgeo.Polygon(c)
                if dl:
                    poly_bnds = poly.bounds
                    if poly_bnds[0] > -10 and poly_bnds[2] <= 180:
                        clip_box = sgeo.box(wl,sl,180,nl)
                    else:
                        clip_box = sgeo.box(-180,sl,el,nl)

                clipped_polys = poly.intersection(clip_box)
                if type(clipped_polys) is sgeo.multipolygon.MultiPolygon:
                    for clipped_poly in clipped_polys:
                        try:
                            bna_polys.append(clipped_poly.exterior.xy)
                        except AttributeError:
                            pass
                elif type(clipped_polys) is sgeo.polygon.Polygon:
                    try:
                        bna_polys.append(clipped_polys.exterior.xy)
                    except AttributeError:
                        pass    
    
    if len(bna_polys) < 1 and lev == 1: #all water
        bna_polys.append(clip_box.exterior.xy)
        levels = [2,]
    else: #define levels and add lakes if specified
        levels = [lev for x in range(len(bna_polys))]

    if lev == 1:   
        bna_polys.insert(0,([wl,wl,el,el],[sl,nl,nl,sl]))
        levels.insert(0,0)
        
    return bna_polys, levels