from libgoods import get_map
import os

coast_data = 'GSHHS' #GSHHS or NOS
nl = 50
sl = 40
wl = -130
el = -120

if coast_data == 'GSHHS':
    res = 'h' # can be f/h/i/l/c (from highest res to lowest)
    includelakes = True 
    shpfile_dir = 'C:\\Users\\amy.macfadyen\\PyProjects\\GOODS\\trunk\\static\\coastline_data\\gshhs_shp'
    shpfile = os.path.join(shpfile_dir,res,'GSHHS_' + res + '_L1.shp')
    bna_polys,levs = get_map.load_polys_and_clip(shpfile,wl,el,nl,sl,dl=0)             
    if includelakes:
        shpfile = os.path.join(shpfile_dir,res,'GSHHS_' + res + '_L2.shp')
        bna_polys2,levs2 = get_map.load_polys_and_clip(shpfile,wl,el,nl,sl,lev=2)
        bna_polys.extend(bna_polys2)
        levs.extend(levs2)
elif coast_data == 'NOS':
    shpfile_dir = 'C:\\Users\\amy.macfadyen\\PyProjects\\GOODS\\trunk\\static\\coastline_data\\nos_shoreline'
    shpfile = os.path.join(shpfile_dir,'nosshore.shp')
    bna_polys,levs = get_map.load_polys_and_clip(shpfile,wl,el,nl,sl)
    

data = get_map.format_for_bna(bna_polys,levs)

f = open(coast_data + '_coast.bna','w')                   
f.write(data)
f.close()