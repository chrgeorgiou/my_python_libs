
      subroutine openfits(imagename,nunit)
      character*(*) imagename
      irw=0
      iblocksize=0
      istat=0
      call ftnopn(nunit,imagename,irw,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Open failed: status=',istat
         write(0,*) imagename
         stop
      endif
      write(0,*) 'Opened fits image ',imagename
      return
      end

      subroutine closefits(nunit)
      istat=0
      call ftclos(nunit,istat)
      call ftrprt('stderr',istat)
      return
      end
      

      function igetsec(nunit,ix1,iy1,nx,ny,zz,maxx)

c     from previously opened fits file with unit nunit, 
c     extract nx by ny pixel rectangle starting at ix1,iy1. 
c     Result gets put in zz array (which can be
c     declared as a 2d nx by ny array elsewhere).

c     igetsec takes a diagnostic value on return:
c     0 if all OK
c     -999 if nx*ny too large for compiled size of zz array
c     -1 if requested rectangle spills out off the fits image
c     istat code from cfitsio otherwise.
c     Programme stops if image is not 2D.

      dimension zz(maxx),naxlen(7)
      integer fpix(7),lpix(7),incpix(7)
      data fpix,lpix,incpix /21*1/
      data naxlen /7*0/

c      irw=0
c      iblocksize=0
c      istat=0
c      call ftnopn(nunit,'inimage.fits',irw,istat)

c      call openfits('inimage.fits',nunit)

      igetsec=-999
      naxis=0
      istat=0
c      write(*,*) nunit
      call ftgidm(nunit,naxis,istat)
c      write(*,*) 222,nunit,naxis,istat
      call ftrprt('stderr',istat)
      call ftgisz(nunit,naxis,naxlen,istat)
c      write(*,*) 333,nunit,naxis,naxlen,istat
      call ftrprt('stderr',istat)
c      write(*,*) 333
      if (naxis.ne.2) then
         write(0,*) 'iGET2SEC: Not a 2-D fits file: ',imagename
         stop
      endif
      if (nx*ny.gt.maxx) then
         write(0,*) 'iGET2SEC: Image dimensions',nx,ny,' too large.'
         return
      endif
      fpix(1)=ix1
      fpix(2)=iy1
      lpix(1)=ix1+nx-1
      lpix(2)=iy1+ny-1
      if (fpix(1).lt.1 .or. lpix(1).gt.naxlen(1) .or.
     :     fpix(2).lt.1 .or. lpix(2).gt.naxlen(2)) then
         write(0,*) 'Image section exceeds boundary'
         igetsec=-1
         return
      endif
c      write(*,*) fpix,lpix,incpix
      call ftgsve(nunit,0,naxis,naxlen,fpix,lpix,incpix,0,zz,anyf,istat)
      call ftrprt('stderr',istat)
      write(0,*) 'Read',nx,' x ',ny,' image section.'
      igetsec=istat
      call ftclos(nunit,istat)
      call ftrprt('stderr',istat)
      return
      end
