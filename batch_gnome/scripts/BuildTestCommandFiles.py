#!/usr/bin/env python

import os
from datetime import datetime, timedelta

import BatchGnome

cf = open(r"c:\DeepWaterHorizon\Test_2009\Machine1\command.txt", 'w')
starttimes = open(r"c:\DeepWaterHorizon\Test_2009\start_times.txt", 'w')


cf.write(r"""[GNOME COMMAND FILE]
-- This Command File was written by the BatchGnome.py module --  
-- It is set up to use a GNOME save file that had been pre set up
--
-- this would clear all maps etc in case we are not launching 
MESSAGE close; TO model; 
--
--load the command file
Message open; TO model; PATH c:\DeepWaterHorizon\Test_2009\nowind2009.SAV
--
-- Now the actual runs
MESSAGE createSpill;TO model; NAME TEST_spill; runDurationInHrs 2880; numLEs 10000; windageA 0.0001; windageB 0.0300; startRelTime   6,   4,  1996,  12,   0; endRelTime   5,   7,  1996,  12,   0; startRelPos 86.5 W 27.0 N; endRelPos 86.5 W 28.75 N;
""")

start_date = datetime(1999, 12, 1)
end_date = datetime(2008, 5, 5)
run_length = timedelta(days = 60)
delta = timedelta(days=14)

path = r"c:\DeepWaterHorizon\Test_2009\trajectories\Spring\001"

start = start_date - delta
run_num = 0
while start < end_date-run_length:
    start += delta
    CS = BatchGnome.ChangeSpill('TEST_spill',
                                StartTime = start,
                                EndTime = start,
                                )
    cf.write(CS.CommandString())
    r = BatchGnome.JustRun(run_length,
                           start,
                           os.path.join(path, "time%03i.bin"%run_num),
                           )
    cf.write(r.CommandString())
    run_num += 1
    starttimes.write(str(start) + '\n')

cf.write("--MESSAGE quit; TO model;\n")


