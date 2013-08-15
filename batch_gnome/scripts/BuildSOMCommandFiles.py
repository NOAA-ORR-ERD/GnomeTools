#!/usr/bin/env python

"""
Special purpose script for DeepWater Horizon event operations.

Probably no other use for this.
"""

import os
from datetime import datetime, timedelta
import BatchGnome


run_length = timedelta(days = 60)
traj_path = r"c:\DeepWaterHorizon\EndRun\trajectories\Spring\001"

init_le_file = r"c:\DeepWaterHorizon\EndRun\Splots_07_13_FORCST"
spill_name = "Init_07_13"

def build_cmd_file(cfpath, startsfile, savefile, start_run_num):
    cf = open(cfpath, 'w')

    starttimes = open(startsfile)

    cf.write(r"""[GNOME COMMAND FILE]
-- this clears all maps etc in case we are not launching 
MESSAGE close; TO model; 
--
-- load the command file
""")
    cf.write("Message open; TO model; PATH %s\n--\n"%savefile)
    cf.write("--\nMESSAGE clearSpills; TO model;\n--\n")
    cf.write("-- Now the actual runs\n")

    spill = BatchGnome.SpillFromFile(spill_name,
                                     init_le_file,
                                     Windage=(0.00001, 3.0),
                                     )
    
    cf.write(spill.CommandString())

    run_num = start_run_num
    
    for line in starttimes:
        start = datetime.strptime(line.split()[0], "%Y-%m-%d")
        CS = BatchGnome.ChangeSpill(spill_name,
                                    StartTime = start,
                                    )
        cf.write(CS.CommandString())
        r = BatchGnome.JustRun(run_length,
                               start,
                               os.path.join(traj_path, "time%03i.bin"%run_num),
                               )
        cf.write(r.CommandString())
        run_num += 1
    cf.write("MESSAGE quit; TO model;\n")
    
    return run_num

run_num = build_cmd_file(r"c:\DeepWaterHorizon\EndRun\machine1\command.txt",
                         r"c:\DeepWaterHorizon\EndRun\SOM_starts2003.txt",
                         r"c:\DeepWaterHorizon\EndRun\EastFlorida_2003_C.SAV",
                         0)
build_cmd_file(r"c:\DeepWaterHorizon\EndRun\machine2\command.txt",
               r"c:\DeepWaterHorizon\EndRun\SOM_starts2009.txt",
               r"c:\DeepWaterHorizon\EndRun\EastFlorida_2009_C.SAV",
               run_num)
