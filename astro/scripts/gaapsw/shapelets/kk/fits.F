

      subroutine openfits(imagename,nunit,nx,ny)
      character*(*) imagename
      dimension naxlen(7)
      irw=0
      iblocksize=0
      istat=0
      call ftnopn(nunit,imagename,irw,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Open failed: status=',istat,', ',imagename
         stop
      endif
      call ftgidm(nunit,naxis,istat)
      call ftrprt('stderr',istat)
      call ftgisz(nunit,naxis,naxlen,istat)
      call ftrprt('stderr',istat)
      if (naxis.ne.2) then
         write(0,*) 'Not a 2D image: ',imagename
         stop
      endif
      nx=naxlen(1)
      ny=naxlen(2)
      write(0,*) 'Opened ',nx,' by ',ny,' image ',imagename
      return
      end

      subroutine closefits(nunit)
      istat=0
      call ftclos(nunit,istat)
      call ftrprt('stderr',istat)
      return
      end
      
      subroutine copyfits(nunit1,nunit2,image2,morekeys)
      character image2*(*)
c copy opened fits file (unit unit1) to a new file, image2 with unit 2.
c open this image. Headers copied as well.
c morekeys specifies for how many new keywords space is reserved
      open(unit=nunit2,file=image2,status='old',err=1)
      close(nunit2,status='delete')
      write(0,*) 'Deleted existing file ',image2
    1 continue
      istat=0
      iblocksize=1
      call ftinit(nunit2,image2,iblocksize,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Open new image failed: status=',istat
         write(0,*) image2
         stop
      endif
      call ftcopy(nunit1,nunit2,morekeys,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Copy to new image failed: status=',istat
         write(0,*) image2
         stop
      endif
      return
      end

      subroutine fitscomment(nunit,comment)
      character comment*(*)
      istat=0
      call ftpcom(nunit,comment,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Write comment failed: status=',istat
         write(0,*) 'image unit ',nunit
         stop
      endif
      return
      end

      subroutine fitshistory(nunit,comment)
      character comment*(*)
      istat=0
      call ftphis(nunit,comment,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Write history failed: status=',istat
         write(0,*) 'image unit ',nunit
         stop
      endif
      return
      end

      function fitskeyr(nunit,keyword)
c read value of a REAL type keyword
      character*(*) keyword,comment*80
      istat=0
      call ftgkye(nunit,keyword,value,comment,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Read REAL keyword failed: status=',istat
         write(0,*) 'image unit ',nunit
         write(0,*) 'keyword ',keyword
         write(0,*) 'Value ',value
         write(0,*) 'Comment ',comment
         stop
      endif
      fitskeyr=value
      return
      end


      subroutine fitswritekey(nunit,keyword,value,comment)
      character*(*) keyword,comment
      istat=0
      call ftpkye(nunit,keyword,value,-6,comment,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Write keyword failed: status=',istat
         write(0,*) 'image unit ',nunit
         write(0,*) 'keyword ',keyword
         write(0,*) 'Value ',value
         write(0,*) 'Comment ',comment
         stop
      endif
      return
      end

      function igetsec(nunit,nxtot,nytot,ix1,iy1,nx,ny,zz,maxx)

c     from previously opened fits file with unit nunit, size nxtot x nytot,
c     extract nx by ny pixel rectangle starting at ix1,iy1. 
c     Result gets put in 1D zz array (which can be
c     declared as a 2d nx by ny array elsewhere).

c     igetsec takes a diagnostic value on return:
c     0 if all OK
c     -999 if nx*ny too large for compiled size of zz array
c     -1 if requested rectangle spills out off the fits image
c     istat code from cfitsio otherwise.

      dimension zz(maxx),naxlen(7)
      integer fpix(7),lpix(7),incpix(7)
      data fpix,lpix,incpix /21*1/
      data naxlen /7*0/
      logical anyf

      istat=0
      igetsec=0
      naxlen(1)=nxtot
      naxlen(2)=nytot
      if (nx*ny.gt.maxx) then
         write(0,*) 'iGETSEC: patch dimensions',nx,ny,' too large.'
         igetsec=-999
         return
      endif
      fpix(1)=ix1
      fpix(2)=iy1
      lpix(1)=ix1+nx-1
      lpix(2)=iy1+ny-1
      if (fpix(1).lt.1 .or. lpix(1).gt.naxlen(1) .or.
     :     fpix(2).lt.1 .or. lpix(2).gt.naxlen(2)) then
         write(0,*) 'iGETSEC: Image section exceeds boundary',ix1,nx,nxtot,iy1,ny,nytot
         igetsec=-1
         return
      endif
c      write(0,*) fpix,lpix,incpix
      call ftgsve(nunit,0,2,naxlen,fpix,lpix,incpix,0,zz,anyf,istat)
      call ftrprt('stderr',istat)
c      write(0,*) 'Read',nx,' x ',ny,' image section.'
      igetsec=istat
      return
      end


      function guess_errors_sec(nunit,ntot)
      parameter (nprime=10099)
      real sky(nprime)
      logical anyf
c take prime number to avoid many samples form same col or row.

c sample image at 100x100 points, look at the median and lower quartile
c assign a uniform sigma based on that

c subsample from the ntot pixel 1-d data stream
      k=0
      istat=0
      if (ntot.le.nprime) then
         i1=1
         istep=1
         nstep=ntot
      else
         istep=ntot/nprime 
         i1=(ntot-istep*(nprime-1))/2 + 1
         nstep=nprime
      endif
      do i=i1,i1+(nstep-1)*istep,istep
         k=k+1
         call ftgpve(nunit,0,i,1,0.,sky(k),anyf,istat)
c do not include values exactly equal to zero
c (should do something better here with masks/weights...)
         if (sky(k).eq.0.) k=k-1
c         call ftrprt('stderr',istat)
c         write(0,*) i,k,sky(k)
      enddo
      call sort(k,sky)
      bg=0.5*(sky(k/2+1)+sky((k+1)/2))
      q1=sky(nint(0.75+0.25*k))
      qdev=bg-q1
      sig=qdev/0.67
      write(0,*) 'Estimated bg, noise level from sample of ',k,
     :    ' pixels: ',bg,sig
      guess_errors_sec=sig
      call pgbeg(0,'noise-est.ps/ps',1,1)
      plotmin=sky(nint(0.05*k)+1)
      plotmax=sky(nint(0.95*k))
      plotwidth=plotmax-plotmin
      plotmin=plotmin-0.1*plotwidth
      plotmax=plotmax+0.1*plotwidth
      call pghist(k, sky, plotmin,plotmax, 100, 2)
      call pgsls(3)
      call pgmove(bg,0.)
      call pgdraw(bg,real(k))
      call pgmove(bg-sig,0.)
      call pgdraw(bg-sig,real(k))
      call pgsls(1)
      call pglab('PIXEL VALUES','#','BG AND NOISE ESIMATE')
      call pgend
      return
      end


      function get_bg_sec(nunit,ntot)
      parameter (nprime=10099)
      real sky(nprime)
      logical anyf
c take prime number to avoid many samples form same col or row.

c sample image at 100x100 points, look at the median and lower quartile
c assign a uniform sigma based on that

c subsample from the ntot pixel 1-d data stream
      k=0
      istat=0
      if (ntot.le.nprime) then
         i1=1
         istep=1
         nstep=ntot
      else
         istep=ntot/nprime 
         i1=(ntot-istep*(nprime-1))/2 + 1
         nstep=nprime
      endif
      do i=i1,i1+(nstep-1)*istep,istep
         k=k+1
         call ftgpve(nunit,0,i,1,0.,sky(k),anyf,istat)
c do not include values exactly equal to zero
c (should do something better here with masks/weights...)
         if (sky(k).eq.0.) k=k-1
c         call ftrprt('stderr',istat)
c         write(0,*) i,k,sky(k)
      enddo
      call sort(k,sky)
      bg=0.5*(sky(k/2+1)+sky((k+1)/2))
      write(0,*) 'Estimated bg level from sample of ',k,
     :    ' pixels: ',bg
      get_bg_sec=bg
      return
      end




c get image and covariance matrix; estimate the latter
      subroutine get2cov (im,nx,ny,zz,maxx,maxy,cov)
      character im*(*)
      real zz(maxx,maxy),cov(3,3)

c first open the image with pixel data, read into zz array
      call get2(im,nx,ny,zz,maxx,maxy)
      write(0,*) 'Opened image'
c then guess variance and covariances
      call guess_cov(zz,cov,nx,ny,maxx,maxy)
      write(0,*) 'got covariances'

      return
      end

      subroutine guess_cov(zz,cov,nx,ny,maxx,maxy)
      dimension zz(maxx,maxy),sky(10000),sky10(10000),
     :    sky01(10000),sky11(10000),cov(3,3)

c sample image at 100x100 points, look at the median and lower quartile
c assign a uniform sigma based on that

      write(0,*) 'Guessing errors and covar from pix value histogram'
      k=0
      do i=1+nx/100,nx-nx/100,1+nx/100
          do j=1+ny/100,ny-ny/100,1+ny/100
              k=k+1
              sky(k)=zz(i,j)
              sky10(k)=zz(i+1,j)-zz(i,j)
              sky01(k)=zz(i,j+1)-zz(i,j)
              sky11(k)=zz(i+1,j+1)-zz(i,j)
          enddo
      enddo
      call sort(k,sky)
      bg=0.5*(sky(k/2+1)+sky((k+1)/2))
      q1=sky(nint(0.75+0.25*k))
      qdev=bg-q1
      sig=qdev/0.67
      cov(2,2)=sig**2
c calculate covariances from var(a-b)=var(a)+var(b)-2cov(a,b); use sqir to
c estimate variances. Median should be zero for a-b vectors.
      call sort(k,sky10)
      qdev=sky10(nint(0.75*k))-sky10(nint(0.25*k))
      vardiff=(qdev/0.67/2)**2
      write(0,*) 'var 10',vardiff
      cov(3,2)=sig**2-vardiff/2
      call sort(k,sky01)
      qdev=sky01(nint(0.75*k))-sky01(nint(0.25*k))
      vardiff=(qdev/0.67/2)**2
      write(0,*) 'var 01',vardiff
      cov(2,3)=sig**2-vardiff/2
      call sort(k,sky11)
      qdev=sky11(nint(0.75*k))-sky11(nint(0.25*k))
      vardiff=(qdev/0.67/2)**2
      write(0,*) 'var 11',vardiff
      cov(3,3)=sig**2-vardiff/2
      cov(1,1)=cov(3,3)
      cov(1,3)=cov(3,3)
      cov(3,1)=cov(3,3)
      cov(1,2)=cov(3,2)
      cov(2,1)=cov(2,3)
      write(0,*) 'Assigned error value of ',sig,' to all pixels.'
      write(0,*) 'Correlations with 10,01,11 pixels:',
     :    cov(3,2)/sig**2,cov(2,3)/sig**2,cov(3,3)/sig**2
      write(0,1010) cov
 1010 format(3f10.2)
      return
      end



c get image and error image; if error image does not exist, fake one
      subroutine get2err (im,imerr,nx,ny,zz,zerr,maxx,maxy)
      character im*(*),imerr*(*)
      real zz(maxx,maxy),zerr(maxx,maxy)

c first open the image with pixel data, read into zz array
      call get2(im,nx,ny,zz,maxx,maxy)

c then try to open the error image
      open(49,file=imerr,status='old',iostat=ios)
      if (ios.ne.0) then
c error image cannot be opened - fake errors
          call guess_errors(zz,zerr,nx,ny,maxx,maxy)
      else
c error image can be opened - read it.
          close(49)
          call get2(imerr,nxerr,nyerr,zerr,maxx,maxy)
c if different size from main image, ignore and estimate errors instead
          if (nxerr.ne.nx .or. nyerr.ne.ny) then
              write(0,*) '***WARNING: ERROR IMAGE SIZE <> IMG SIZE!!***'
              call guess_errors(zz,zerr,nx,ny,maxx,maxy)
          endif
      endif
      return
      end


      subroutine guess_errors(zz,zerr,nx,ny,maxx,maxy)
      real zz(maxx,maxy),zerr(maxx,maxy),sky(10000)

c sample image at 100x100 points, look at the median and lower quartile
c assign a uniform sigma based on that

      write(0,*) 
     :     'No error image - guessing errors from pix value histogram'
      k=0
      do i=1,nx,1+nx/100
          do j=1,ny,1+ny/100
              k=k+1
              sky(k)=zz(i,j)
          enddo
      enddo
      call sort(k,sky)
      bg=0.5*(sky(k/2+1)+sky((k+1)/2))
      q1=sky(nint(0.75+0.25*k))
      qdev=bg-q1
      sig=qdev/0.67
      do i=1,nx
          do j=1,ny
              zerr(i,j)=sig
          enddo
      enddo
      write(0,*) 'Assigned error value of ',sig,' to all pixels.'
      return
      end



      subroutine get2 (imagename,nx,ny,zz,maxx,maxy)
c this is an interface to the fitsio routines.
c open image `imagename' (fits or imh),
c and put it into the array data, with declared
c dimensions (maxx,maxy). Also return image size as nx,ny.
      character*(*) imagename
      dimension zz(maxx,maxy),naxlen(7)
      logical anyf

      nunit=49
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
      call ftgidm(nunit,naxis,istat)
      call ftgisz(nunit,naxis,naxlen,istat)
      call ftrprt('stderr',istat)
      if (naxis.eq.1) then
         nx=naxlen(1)
         ny=1
         if (nx.gt.maxx) then
            write(0,*) 'GET2: Image dimensions',nx,' too large.'
            stop
         endif
         call ftgpve(nunit,0,1,nx,0.,zz,anyf,istat)
         call ftrprt('stderr',istat)
      else
         nx=naxlen(1)
         ny=naxlen(2)
         if (nx.gt.maxx) then
            write(0,*) 'GET2: Image X-dimension',nx,' too large.'
            stop
         endif
         if (ny.gt.maxy) then
            write(0,*) 'GET2: Image Y-dimension',ny,' too large.'
            stop
         endif
         write(0,*) 'Read',nx,' x ',ny,' image ',imagename
         call ftg2de(nunit,0,0.,maxx,nx,ny,zz,anyf,istat)
         call ftrprt('stderr',istat)
      endif
      call ftclos(nunit,istat)
      call ftrprt('stderr',istat)
      return
      end


      subroutine newfits(image,nunit,nx,ny)
      integer naxlen(7)
      data naxlen /7*0/
      character image*(*)
      logical simple, extend
      open(unit=nunit,file=image,status='old',err=1)
      close(nunit,status='delete')
      write(0,*) 'Deleted existing file ',image
    1 continue
      irw=0
      iblocksize=1
      istat=0
      simple=.true.
      extend=.false.
      naxis=2
      naxlen(1)=nx
      naxlen(2)=ny
c open new fits file
      call ftinit(nunit,image,iblocksize, istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Open new image failed: status=',istat
         write(0,*) image
         stop
      endif
c write required header info
c      call ftphpr(nunit,simple,-32,naxis,naxlen,0,1,extend,istat)
c      if (istat.ne.0) then
c         call ftrprt('stderr',istat)
c         write(0,*) 'FITS: Write header failed: status=',istat
c         write(0,*) image
c         stop
c      endif
c create primary array
      call ftiimg(nunit,-32,naxis,naxlen,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Create primary array failed: status=',istat
         write(0,*) image
         stop
      endif
c close the image and reopen it - this forces data section to be written out
      call ftclos(nunit,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Close new image failed: status=',istat
         write(0,*) image
         stop
      endif
      call ftopen(nunit,image,1,iblocksize,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Reopen new image failed: status=',istat
         write(0,*) image
         stop
      endif
      write(0,*) 'Succesfully opened new image ',image
      return
      end


      subroutine putsec(nunit,pix,nx,ny,ix1,iy1,naxlen)
c write nx,ny array pix into open fits image with unit nunit.
c bottom left pixel goes to ix1,iy1
c array pix must have dimension nx,ny ! No gaps or borders
      dimension pix(nx,ny)
      dimension naxlen(7),ifpix(7),ilpix(7)
      data ifpix,ilpix /14*0/

c write the data.
c      write(0,*) 'putsec: ',nunit,naxlen
c      write(0,*) 'putsec: ',pix(1,1),pix(nx,ny)
      ndim=2
      igroup=1
      ifpix(1)=ix1
      ilpix(1)=ix1+nx-1
      ifpix(2)=iy1
      ilpix(2)=iy1+ny-1
      istat=0
      call ftpsse(nunit,igroup,ndim,naxlen,ifpix,ilpix,pix,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Putsec: Write data failed: status=',istat
         stop
      endif
      end


      subroutine put2(image,nx,ny,pix,mx,my)
      real pix(mx,my)
      integer naxlen(7)
      character image*(*)
      logical simple,extend

c delete the file if it already exists

c      call deletefile(image,istat)

      nunit=49
      open(unit=nunit,file=image,status='old',err=1)
      close(nunit,status='delete')
      write(0,*) 'Deleted existing file ',image
    1 continue
      irw=0
      iblocksize=1
      istat=0
      simple=.true.
      extend=.false.
      naxis=2
      naxlen(1)=nx
      naxlen(2)=ny
c open new fits file
      call ftinit(nunit,image,iblocksize, istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Put2: Open new image failed: status=',istat
         write(0,*) image
         stop
      endif
c write required header info
      call ftphpr(nunit,simple,-32,naxis,naxlen,0,1,extend,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Put2: Write header failed: status=',istat
         write(0,*) image
         stop
      endif
c write the data.
      igroup=1
      call ftp2de(nunit,igroup,mx,nx,ny,pix,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Put2: Write data failed: status=',istat
         write(0,*) image
         stop
      endif
c close the image
      call ftclos(nunit, istat)
      end




      subroutine put3(image,nx,ny,nz,pix,mx,my,mz)
      real pix(mx,my,mz)
      integer naxlen(7)
      character image*(*)
      logical simple,extend

c delete the file if it already exists

c      call deletefile(image,istat)

      nunit=49
      open(unit=nunit,file=image,status='old',err=1)
      close(nunit,status='delete')
      write(0,*) 'Deleted existing file ',image
    1 continue
      irw=0
      iblocksize=1
      istat=0
      simple=.true.
      extend=.false.
      naxis=3
      naxlen(1)=nx
      naxlen(2)=ny
      naxlen(3)=nz
c open new fits file
      call ftinit(nunit,image,iblocksize, istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Put3: Open new image failed: status=',istat
         write(0,*) image
         stop
      endif
c write required header info
      call ftphpr(nunit,simple,-32,naxis,naxlen,0,1,extend,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Put3: Write header failed: status=',istat
         write(0,*) image
         stop
      endif
c write the data.
      igroup=1
      call ftp3de(nunit,igroup,mx,my,nx,ny,nz,pix,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Put3: Write data failed: status=',istat
         write(0,*) image
         stop
      endif
c close the image
      call ftclos(nunit, istat)
      end





      subroutine get3 (imagename,nx,ny,nz,zz,maxx,maxy,maxz)
c this is an interface to the fitsio routines.
c open image `imagename' (fits or imh),
c and put it into the 3D array data zz, with declared
c dimensions (maxx,maxy,maxz). Also return image size as nx,ny,nz.
      character*(*) imagename
      dimension zz(maxx,maxy,maxz),naxlen(7)
      logical anyf

      nunit=49
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
      call ftgidm(nunit,naxis,istat)
      call ftgisz(nunit,naxis,naxlen,istat)
      call ftrprt('stderr',istat)
      if (naxis.eq.1) then
         nx=naxlen(1)
         ny=1
         nz=1
         if (nx.gt.maxx) then
            write(0,*) 'GET3: Image dimensions',nx,' too large.'
            stop
         endif
         write(0,*) 'Read',nx,' x ',1,' image ',imagename
         call ftgpve(nunit,0,1,nx,0.,zz,anyf,istat)
         call ftrprt('stderr',istat)
      elseif (naxis.eq.2) then
         nx=naxlen(1)
         ny=naxlen(2)
         nz=1
         if (nx.gt.maxx) then
            write(0,*) 'GET3: Image X-dimension',nx,' too large.'
            stop
         endif
         if (ny.gt.maxy) then
            write(0,*) 'GET3: Image Y-dimension',ny,' too large.'
            stop
         endif
         write(0,*) 'Read',nx,' x ',ny,' image ',imagename
         call ftg2de(nunit,0,0.,maxx,nx,ny,zz,anyf,istat)
         call ftrprt('stderr',istat)
      else
         nx=naxlen(1)
         ny=naxlen(2)
         nz=naxlen(3)
         if (nx.gt.maxx) then
            write(0,*) 'GET3: Image X-dimension',nx,' too large.'
            stop
         endif
         if (ny.gt.maxy) then
            write(0,*) 'GET3: Image Y-dimension',ny,' too large.'
            stop
         endif
         if (nz.gt.maxz) then
            write(0,*) 'GET3: Image Y-dimension',ny,' too large.'
            stop
         endif
         write(0,*) 'Read',nx,' x ',ny,' x ',nz,' image ',imagename
         call ftg3de(nunit,0,0.,maxx,maxy,nx,ny,nz,zz,anyf,istat)
         call ftrprt('stderr',istat)
      endif
      call ftclos(nunit,istat)
      call ftrprt('stderr',istat)
      return
      end



      subroutine tabget(imagename,extname,nrow,coltempl,zz,mrow)
      real zz(mrow)
      character*(*) coltempl,extname
      character*(*) imagename
c
c find the table in the fits file imagename with extension named extname,
c and read column with name that matched coltempl into array zz.
c also return number of rows in the column.
c NB this routine opens & closes the file each time it is called, not very
c efficient.
c
      nunit=47
      irw=0
      iblocksize=0
      istat=0
      iextver=1
      call ftopen(nunit,imagename,irw,iblocksize,istat)
      call ftrprt('stderr',istat)
      call ftmnhd(nunit,-1,extname,iextver,istat)
      call ftrprt('stderr',istat)
      call ftgnrw(nunit,nrow,istat)
      call ftrprt('stderr',istat)
      call ftgcno(nunit,.false.,coltempl,icol,istat)
      call ftrprt('stderr',istat)
      call ftgcve(nunit,icol,1,1,nrow,-999.,zz,anyf,istat)
      call ftrprt('stderr',istat)
      call ftclos(nunit,istat)
      return
      end
