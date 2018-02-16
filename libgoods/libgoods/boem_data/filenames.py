#!/usr/bin/env python

"""
computes the file names we need
"""

import sys, datetime

days_60 = datetime.timedelta(60)
days_30 = datetime.timedelta(30)

#months = (3,4,5,6,7,8,9)
months = range(1,13)

#res = 'course'
#res = 'fine'
res = sys.argv[1]

if res == 'course':
    start_date = datetime.datetime(1992, 8, 29)
    outfile = file("course_filename.txt", 'w')
    for num in range(5, 50):
        filename = "SURF_GOM27_%02i"%(num)
        date = start_date + (num-5)*days_60    
        print "file: %s -- date: %s"%(filename, date)
        if date.month in months:
            nc_filename = "%i_%s.nc"%(date.year,filename)
            print "we need: %s, %s"%(filename, nc_filename)
            outfile.write("%s,%s\n"%(filename, nc_filename))

elif res == 'fine':
    start_date = datetime.datetime(1999, 10, 22)

    outfile = file("fine_filename.txt", 'w')
    outfile2 = file("fine_filetimes.txt",'w')
    for num in range(47, 99):
        for num2 in (2, 3):
            filename = "SURF_GOM27_%i_0%s"%(num, num2)
            date = start_date + (num-47)*days_60 + (num2-2)*days_30    
            print "file: %s -- date: %s"%(filename, date)
            outfile2.write("%s -- %s\n"%(filename, date))
            if date.month in months:
                nc_filename = "%i_%s.nc"%(date.year,filename)
                print "we need: %s, %s"%(filename, nc_filename)
                outfile.write("%s,%s\n"%(filename, nc_filename))
    
    
