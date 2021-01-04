from __future__ import print_function
from libgoods import rect_model, data_files_dir
import numpy as np
import os
reload(rect_model)

acnfs_dir = os.path.join(data_files_dir,'acnfs_nov15')
download_file = 1

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
    if fn.split('_')[-1] in ['t000.nc','t012.nc','t024.nc','t036.nc','t048.nc','t060.nc','t072.nc']:
        flist.append(os.path.join(acnfs_dir,fn))

bbox = [69.8,-150.7,71.4,-143.4] #Geographic domain [South Lat, West Lon, North Lat, East Lon]
out_dir = data_files_dir #Where to write files (in this case libgoods/data_files )

var_map = { 'time':'time',
            'lon': 'lon',
            'lat': 'lat',
            'u_velocity': 'uocn',
            'v_velocity': 'vocn',
            'ice_u': 'uvel',
            'ice_v': 'vvel',
            'ice_thickness': 'hi',
            'ice_fraction': 'aice',
            } 
            
acnfs = rect_model.rect(flist)
acnfs.get_dimensions(var_map)


#Determine geographic subset indices and get data
acnfs.subset(bbox) #south lat, west lon, north lat, east lon

acnfs.write_nc(var_map,os.path.join(acnfs_dir,'acnfs.nc'))


  
