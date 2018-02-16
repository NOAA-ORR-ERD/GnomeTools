       PROGRAM rdgom_27_hdr_only
       IMPLICIT NONE
       CHARACTER*80 filename
       LOGICAL LEXIST
       PARAMETER is=1,ie=200,im=ie-is+1,js=1,je=200,jm=je-js+1
       PARAMETER kb=1,irecord=4*4+4*5*im*jm
       INTEGER  itime1(4)
       REAL TIME
       REAL z(kb),zz(kb),dz(kb),dzz(kb)
       REAL h(im,jm),dx(im,jm),dy(im,jm),alon(im,jm),
     &      alat(im,jm),dum(im,jm),dvm(im,jm),fsm(im,jm),cor(im,jm)
       INTEGER i,igs,jgs,ige,jge,ip,jp
c
c      compile with
c      [...] ifort /assume:byterecl /convert:big_endian rdgom_hdr_only.f
c
       filename='SURF_GOM27_05'

        OPEN (10,file='gom_hdr_list.txt',
     &     FORM='FORMATTED',
     &     status='NEW')

        OPEN (11,file=filename,
     &     ACCESS='DIRECT',FORM='UNFORMATTED',
     &     RECL=irecord,status='UNKNOWN')

       INQUIRE(file=filename ,exist=lexist)
        IF (.not.lexist) THEN
          PRINT *,filename,'is not exist'
        ELSE

          READ(11,rec=481)itime1,alon,alat,h,cor,z,zz,dz,dzz

		do jp=1,jm
		do ip=1,im
		write(10,107) alon(ip,jp),alat(ip,jp),h(ip,jp)
		end do
		end do
107    format(f12.4,',',f12.4,',',f12.4)
        ENDIF

      END
