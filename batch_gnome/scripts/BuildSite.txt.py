#!/usr/bin/env python

from TAP_Setup import setup
import os

print "starting"

file = open(os.path.join(setup.RootDir,"site.txt"),'w')

file.write('"%s"  //Map Name for User\n'%setup.MapName)
file.write('"%s"  "%s"  // "map file name"  "map file type" -- either "BNA" or "PICT" \n'%(setup.MapFileName, setup.MapFileType) )
file.write('%i SPILLS // number of spills per site and time in every cube  \n'%setup.NumStarts )
file.write('%i SEASONS\n'%len(setup.StartTimeFiles) )
for i in range(len(setup.Seasons)):
    Season = setup.Seasons[i][0]
    file.write('"%s" "%s" "%s" // Name for menu, prefix of cube names, name of directory\n'%
               (Season, setup.CubesRootNames[i], Season) )
file.write('%i TIMES \n'%(len(setup.OutputTimes)) )
#reverse output times:
setup.OutputTimes.reverse()
setup.OutputUserStrings.reverse()
for OutputTime, OutputUserString in zip(setup.OutputTimes, setup.OutputUserStrings):
    file.write('%i "%s" // numHours  "user string" blank user string means not in menu \n'%(
               OutputTime, OutputUserString) )
file.write('%i AMOUNTS // amounts of oil preset for the user\n'%len(setup.PresetSpillAmounts) )
for s in setup.PresetSpillAmounts:
    file.write('%s  \n'%s )
               
file.write('%i LOCS // levels of concern (LOC) preset for the user \n'%len(setup.PresetLOCS) )
for s in setup.PresetLOCS:
    file.write('%s  \n'%s )

# now write the BNA of the receptor sites:
if setup.ReceptorType == "Grid":
    #file.write("This is where the receptors go\n")
    from batch_gnome import tap_mod
    Grid = setup.Grid
    Receptors = tap_mod.Grid(Grid.min_long, Grid.max_long, Grid.min_lat, Grid.max_lat,Grid.num_lat,Grid.num_long)
    file.write("%i SITES // number of sites\n"%Receptors.num_cells)
    Receptors.WriteBNA(file)

# write cube locations:
file.write("%i CUBES\n"%len(setup.CubeStartSites))

print "Writing the start sites"
for site in setup.CubeStartSites:
    file.write(site + "\n")
file.close()    







