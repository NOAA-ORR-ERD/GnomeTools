from __future__ import print_function
from libgoods import curv_grid, data_files_dir
import numpy as np
import os

acnfs_dir = os.path.join(data_files_dir,'acnfs')
download_file = 0

if not os.path.exists(acnfs_dir):
    os.mkdir(acnfs_dir)
if download_file:
    #Get latest file from ACNFS FTP site
    import ftplib
    import tarfile
    ftp_path = 'tgftp.nws.noaa.gov'
    ftp = ftplib.FTP(ftp_path)   # connect to host, default port
    ftp.login()
    ftp.cwd('/SL.us008001/CS.navo/MT.arctic/')               
    try:
        files = ftp.nlst()
    except ftplib.error_perm as resp:
        if str(resp) == "550 No files found":
            print("No files in this directory")
        else:
            print("Failed to download ACNFS file")
            raise
    ofn = os.path.join(acnfs_dir,files[0])
    ftp.retrbinary("RETR " + files[0] , open(ofn,'wb').write)
    ftp.quit()
    tar = tarfile.open(ofn, "r:gz")
    tar.extractall(path=acnfs_dir)
    tar.close()
else:
    pass

fns = os.listdir(acnfs_dir)
flist = []
for fn in fns:
    if fn.split('_')[-1] in ['t000.nc','t024.nc','t048.nc','t072.nc']:
        flist.append(os.path.join(acnfs_dir,fn))
print(flist)


bbox = [69,-172,72.5,-155] #Geographic domain [South Lat, West Lon, North Lat, East Lon]
out_dir = data_files_dir #Where to write files (in this case libgoods/data_files )

var_map = { 'time':'time',
            'longitude': 'lon',
            'latitude': 'lat',
            'u_velocity': 'uocn',
            'v_velocity': 'vocn',
            } 
            
acnfs = curv_grid.cgrid(flist)
acnfs.get_dimensions(var_map)

#meshgrid lon/lat for curvilinear grid
acnfs.data['lon'],acnfs.data['lat'] = np.meshgrid(acnfs.data['lon'],acnfs.data['lat'])

#Determine geographic subset indices and get data
acnfs.subset(bbox) #south lat, west lon, north lat, east lon
acnfs.get_data(var_map,yindex=acnfs.y,xindex=acnfs.x,zindex=0,is3d=False,extra_2dvars=['hi','aice','uvel','vvel'])     

#make mask
mask = (acnfs.data['u'][0,:,:] == acnfs.atts['u']['missing_value']).choose(1,0)
acnfs.grid['mask'] = mask

#rename ice vars
rename_dict = {'hi':'ice_thickness','aice':'ice_fraction','uvel':'ice_u','vvel':'ice_v'}
for key, val in rename_dict.items():
    acnfs.data[val] = acnfs.data[key]
    acnfs.atts[val] = acnfs.atts[key]

acnfs.write_nc(os.path.join(out_dir,'acnfs_example.nc'),is3d=False,extra_2dvars=rename_dict.values())
  
