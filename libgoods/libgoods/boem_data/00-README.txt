README for MMS data:

Info about the MMS data files:



I uploaded to your ftp site the first file SURF_GOM27_05. About
384807696 bytes.

 I also put three simple programs, and two txt and dat files.

The output of hdr_only is a list of the long, lat, depth at each grid
pt. I did not put this on this email, but it is about 1.5Mb.

The output of gom-time is a list of the 480 times as yr month day hour.

In the source code there a couple of switches specified /assume:byterecl
means the FORTRAN OPEN statement will use bytes as the recl size.  The
/convert:big_endian changes the endian issue. If the output looks bad,
try it the other way.
I have similar switches in Linux Portland Group compiler, if you need
them.

