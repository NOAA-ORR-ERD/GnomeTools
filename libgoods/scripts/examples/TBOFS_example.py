#!/usr/bin/env python
from __future__ import print_function
from libgoods import roms_model, nctools, noaa_coops, data_files_dir
import datetime as dt
import os 
reload(roms_model)

out_dir = os.path.join(data_files_dir,'tbofs')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
# specify local file or opendap url

sdate = dt.date(2020,1,18)
edate = dt.date(2020,1,22)
flist = noaa_coops.make_server_filelist('tbofs',sdate)
print("Number of files to aggregate: ", len(flist))

model = roms_model.roms(flist[0])
print('getting dimensions')
model.get_dimensions()


for f in flist:
    print(f)
    if f != flist[0]:
        model.update(f)
        model.get_dimensions(get_xy=False)
    model.when()
    print('writing file')
    coops_name = f.split('/')[-1].split('.')
    ofn = coops_name[1] + '_' + coops_name[4] + '_' + coops_name[5] + '_' + coops_name[3] + '.'+ coops_name[-1]
    print(ofn)
    model.write_nc(os.path.join(data_files_dir,'tbofs',ofn))
    shift_time(os.path.join(data_files_dir,'tbofs',ofn), -5, tvar='ocean_time')
