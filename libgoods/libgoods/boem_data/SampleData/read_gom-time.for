       program read_gom_time
c!  Little endian
C   compile ifort /assume:byterecl /O2 /convert:big_endian read_gom-time.for
c!       implicit none
       character*35 filename, fileout
       character*3  arg1
       CHARACTER*2 fnum
       parameter  nx=200,ny=200
       PARAMETER irecord=4*4+4*5*nx*ny
       real uo(nx,ny),vo(nx,ny),uw(nx,ny),vw(nx,ny),el(nx,ny)
       real xlat(nx,ny), xlon(nx,ny)
       integer   itime1(4)
c
      common /inter/ xlat,xlon,alat,alon,iweights,pweights
c
c     Input grid is 200 x 200
      Print *, ' Input Grid is nx=',nx, ' ny= ',ny
c
c     read lat-lon grid
c
      Print *, 'Enter file to process (05 thru 49)'
       read(5,2,end=999) fnum
2      format(A2)
c
       filename='SURF_GOM27_'//fnum
        print *, 'Input filename  ', filename
        OPEN (15,file=filename,ACCESS='DIRECT',FORM='UNFORMATTED',
     &        RECL=irecord,status='UNKNOWN')

       INQUIRE(file=filename ,exist=lexist)
        IF (.not.lexist) THEN
          PRINT *,filename,' does not exist'
          stop 7
        ELSE

          fileout='GOM27_times'//fnum//'.dat'
          open(16,file=fileout,form='FORMATTED',status='NEW')

        print *, 'Output filename ', fileout
       ENDIF
c
      k=0
C
      do iday=1,480
c
      k = k+1
      read(15,rec=k) itime1,uo,vo,uw,vw,el
      print *, ' Read  step ',k, ' Time ',itime1
c
	write(16,301) (itime1(ii),ii=1,4)
301   format('  Time =   ', 4I12 )
c
c	do jp = 1, ny
c		do ip = 1, nx
c			write(16,310)uw(ip,jp),vw(ip,jp),uo(ip,jp),vo(ip,jp)
c310   format(2f7.3,2f8.3)
c		end do
c	end do
c
      print *, ' Wrote step ',k
c       iday loop
	end do
c
      close(15)
c
999   print *, 'Done '
c
      stop
      end
c      end main program
