= lib_tap =

A library of scripts, etc. for buiding a  "TAP" location with py_gnome.

This is still in the protype stage.

It has been (lightly) tested through the RunPyGnome stage. The processing of the trajectory files to make cubes needs to be updated.

Moving on, other stuff will need to be pulled from batch_gnome

Amy's notes on the next steps:

 - Needs to be updated to use the recent nc_particles (in post_gnome)
 - Change line 153 (in ???) time -->> time[:]
 - In TAP_mod.py:
   - use tap_comp_volume, rather than cy_tap (line 666) -- or get cy_tap to work!
   - line 754 -- flag --> status-codes

