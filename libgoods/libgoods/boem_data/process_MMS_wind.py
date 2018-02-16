#!/usr/bin/env python

"""
script to process the MMS GOM data

"""

import os, sys, time

import read_gom
    

#for file in files: print file

#raise Exception("stopping here")

def convert_files(mms_path, nc_path, size, filelist, resolution):
    AllDone = True

    files = [line.strip().split(',') for line in file(filelist,'r').readlines()]

    for mms, nc in files:
        mms_name = os.path.join(mms_path, mms)
        nc_name  = os.path.join(nc_path, nc)

        try:
            file_size = os.path.getsize(mms_name)
        except WindowsError:
            AllDone = False
            print "file: %s doesn't exist"%(mms_name)
            continue        
        #print "file: %s is %i bytes"%(mms_name, file_size)
        if file_size == size:
            print "valid file"
            if os.path.exists(nc_name):
                print "%s already exists"%nc_name
            else:
                print "need to process: "
                try:
                    convert_file(mms_name, nc_name, resolution)
                except:
                    AllDone = False
                    print "an exception occuredying to process file"
                    #raise
                    print "trying the next one"
                    pass # go on, don't want to stop the while thing...
        else:
            AllDone = False
            print "partially downloaded"
    return AllDone

def convert_file(infilename, outfilename, resolution):
    print "processing:", infilename
    time.sleep(2) # seconds
    dataset = read_gom.data_reader(infilename, resolution)
    dataset.read_all()
    dataset.to_netcdf_wind(outfilename)


if __name__ == "__main__":
    try:
        if sys.argv[1] == 'course':
            mms_path = r"t:\DeepWaterHorizon\\MMS_DATA\Oey_gom_2003"
            nc_path = r"T:\DeepWaterHorizon\MMS_DATA\nc_wind_2003"
            size = 384807696
            filelist = 'course_filename.txt'
        elif sys.argv[1] == 'fine':
            mms_path   = r"T:\DeepWaterHorizon\MMS_DATA\Oey_gom_2009"
            nc_path    = r"T:\DeepWaterHorizon\MMS_DATA\nc_wind_2009"
            size = 660151056
            filelist = 'fine_filename.txt'
        else:
            raise IndexError
        while True:
            AllDone = convert_files(mms_path, nc_path, size, filelist, sys.argv[1])
            if AllDone:
                break
            else:
                print "not done yet, sleeping for a bit, then trying again"
                time.sleep(300)
                print "Trying again"
    except IndexError:
        print 'you must pass either "course" or "fine" in on the command line'

