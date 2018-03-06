

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

      subroutine copyfits_nodata(image1,nunit2,image2)
      character image2*(*),image1*(*)
c copy fits file (file image1) to a new file, image2 and open as unit 2.
c Headers of image1, not its data, are copied.
c delete existing image2 if necessary
      open(unit=nunit2,file=image2,status='old',err=1)
      close(nunit2,status='delete')
      write(0,*) 'Deleted existing file ',image2
    1 continue
      istat=0
      call fttplt(nunit2,image2,image1,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Creating new image failed: status=',istat
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

      logical function fitshavekey(nunit,keyword)
c asks whether a key exists
      character*(*) keyword,comment*80
      istat=0
      call ftgkye(nunit,keyword,value,comment,istat)
      fitshavekey=(istat.eq.0)
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


      function getonepix(nunit,ix,iy,nx,ny)
c read a single pixel value at (ix,iy) from a 2d nx*ny fits image opened at unit nunit
      logical anyf
      istat=0
      if (ix.le.nx .and. iy.le.ny .and. ix.gt.0 .and. iy.gt.0) then
         call ftgpve(nunit,0,ix+(iy-1)*nx,1,0.,z,anyf,istat)
         getonepix=z
      else
         write(0,*) '*** getonepix attempted to read outside image ***'
         getonepix=0
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


      function igetsecpad(nunit,nxtot,nytot,ix1,iy1,nx,ny,zz,maxx,bg)

c     from previously opened fits file with unit nunit, size nxtot x nytot,
c     extract nx by ny pixel rectangle starting at ix1,iy1.
c     Background subtract all pixels with value bg
c     Result gets put in 1D zz array (which can be
c     declared as a 2d nx by ny array elsewhere).
c     
c     differs from igetsec: if the extracted rectangle overflows the
c     border of the fits file the patch is padded with zeroes, rather
c     than generating an error
c
c     igetsec takes a diagnostic value on return:
c     0 if all OK
c     -999 if nx*ny too large for compiled size of zz array
c     istat code from cfitsio otherwise.

      dimension zz(maxx),naxlen(7)
      integer fpix(7),lpix(7),incpix(7)
      data fpix,lpix,incpix /21*1/
      data naxlen /7*0/
      logical anyf

      istat=0
      igetsecpad=0
      naxlen(1)=nxtot
      naxlen(2)=nytot
      if (nx*ny.gt.maxx) then
         write(0,*) 'iGETSECPAD: patch dimensions',nx,ny,' too large.'
         igetsecpad=-999
         return
      endif
c extract the part of the rectangle that fits inside the fits image
c write it in straight fortran order, no gaps
      fpix(1)=max(1,ix1)
      fpix(2)=max(1,iy1)
      lpix(1)=min(ix1+nx-1,naxlen(1))
      lpix(2)=min(iy1+ny-1,naxlen(2))
ccc      write(0,*) 'Getting section',fpix(1),fpix(2),lpix(1),lpix(2)
      call ftgsve(nunit,0,2,naxlen,fpix,lpix,incpix,0,zz,anyf,istat)
      call ftrprt('stderr',istat)
      igetsecpad=istat
c check whether the rectangle needs to be padded
      if (ix1.lt.1 .or. ix1+nx-1.gt.naxlen(1) .or.
     :     iy1.lt.1 .or. iy1+ny-1.gt.naxlen(2)) then
ccc         write(0,*) 'iGETSECPAD: Image section exceeds boundary',ix1,nx,nxtot,iy1,ny,nytot
ccc         write(0,*) 'Padding with zeros to size ',nx,' x ',ny
c last line/col to zero-pad at beginning
         ixpad1=max(0,1-ix1)
         iypad1=max(0,1-iy1)
c first line/col to zero-pad at the end
         ixpad2=min(nx+1,naxlen(1)-ix1+2)
         iypad2=min(ny+1,naxlen(2)-iy1+2)
ccc         write(0,*) 'Padding borders in x,y:',ixpad1,ixpad2,iypad1,iypad2
c now copy the zz array from high to low index, inserting zeros as we go.
c also do the background subtraction here
c keep two running indices with location where pixels should be written and where they come from
c izin should always be <= izout, so can work from high to low
         izout=nx*ny
         izin=(lpix(1)-fpix(1)+1)*(lpix(2)-fpix(2)+1)
ccc         write(0,*) 'start at izout,izin=',izout,izin
c top lines that need to be zeroed out entirely
         do iy=ny,iypad2,-1
            do ix=nx,1,-1
               zz(izout)=0
               izout=izout-1
            enddo
         enddo
c middle lines that may need to be zero-padded in x
         do iy=iypad2-1,iypad1+1,-1
            do ix=nx,ixpad2,-1
               zz(izout)=0
               izout=izout-1
            enddo
            do ix=ixpad2-1,ixpad1+1,-1
               zz(izout)=zz(izin)-bg
               izout=izout-1
               izin=izin-1
            enddo
            do ix=ixpad1,1,-1
               zz(izout)=0
               izout=izout-1
            enddo
         enddo
c bottom lines that need to be zeroed out entirely
         do iy=iypad1,1,-1
            do ix=nx,1,-1
               zz(izout)=0
               izout=izout-1
            enddo
         enddo
ccc        write(0,*) 'izout, izin:',izout,izin
         if (izout.ne.0 .or. izin.ne.0) then
            write(0,*) 'iGETSECPAD: Something went wrong, ',
     :           'not all pixels correctly copied'
            igetsecpad=-1
         endif
      else
c no need to pad; just du background subtraction
ccc         write(0,*) 'Background subtracting: nx,ny,bg=',nx,ny,bg
         do i=1,nx*ny
            zz(i)=zz(i)-bg
         enddo
      endif
c      write(0,*) fpix,lpix,incpix
c      write(0,*) 'Read',nx,' x ',ny,' image section.'
      return
      end



      function guess_errors_sec(nunit,ntot)
      parameter (nprime=10099)
      real sky(nprime)
      logical anyf
c take prime number to avoid many samples form same col or row.

c sample image at ~ 100x100 points, look at the median and lower quartile
c assign a sigma based on that

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
c find median and lower quartile
      call sort(k,sky)
      call get_med_quartile(k,sky,bg,qdev)
c for a Gaussian, quartile is at 0.675 sigma
      sig=qdev/0.675
      write(0,2000) k,bg,sig
 2000 format (' B/g, noise estimate from sample of ',i5,' pixels:',2g12.4)
c use first estimate to remove values very close to the median (empty
c regions of the image)
      zaplo=bg-0.001*sig
      zaphi=bg+0.001*sig
      ilo=1
c make ilo first item above lower limit
      do while (sky(ilo).lt.zaplo)
          ilo=ilo+1
      enddo
      nkill=0
c make ihi last item below lower limit
      do while (sky(ilo+nkill).lt.zaphi)
          nkill=nkill+1
      enddo
c excise nkill items from the list, starting at ilo
      if (nkill.gt.0) then
          do i=ilo,k-nkill
              sky(i)=sky(i+nkill)
          enddo
      endif
c repeat median and quartile determination:
      k=k-nkill
      call get_med_quartile(k,sky,bg,qdev)
      sig=qdev/0.675
      write(0,2000) k,bg,sig

      guess_errors_sec=sig

      open(27,file='bgnoise.dat',status='unknown')
      write(27,*) bg,sig
      close(27)
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


      function guess_bg_sec(nunit,ntot)
      parameter (nprime=10099)
      real sky(nprime)
      logical anyf
c estimate median value of an image by subsampling 1D pixel stream
c take prime number to avoid many samples form same col or row.

c sample image at ~ 100x100 points, look at the median and lower quartile
c assign a bg based on that

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
         if (sky(k).eq.0.) k=k-1
      enddo

      guess_bg_sec=fmedian(k,sky)
      return
      end

      subroutine get_med_quartile(n,y,ym,qdev)
c from SORTED array y(n), iteratively find median and lower quartile width, in a way that is robust to excess of positive values
c first find true median M and lower quartile Q1
c then find median of all values between Q1 and M + (M-Q1)
c then recalculate Q1 as median of all values below M
c return M and qdev=M-Q1
      dimension y(n)
      ym=(y(n/2+1) + y((n+1)/2)) /2
      iq1=nint(0.25*n+0.75)
      q1=y(iq1)
      q2=2*ym-q1
      iq2=iq1+1
      do while (y(iq2).lt.q2 .and. iq2.lt.n)
         iq2=iq2+1
      enddo
      ym=0.5*(y((iq1+iq2)/2+1)+y((iq1+iq2+1)/2))
      iq1=nint(0.25*(iq1+iq2+0.5))
      q1=y(iq1)
      qdev=ym-q1
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


      subroutine guesscov_sec (nunit,nx,ny,mcov,cov,zero)
      parameter (nprime=10099,mcovmax=10)
      real cov(-mcov:mcov,-mcov:mcov),zz(0:mcov,0:mcov,nprime)
      real a(nprime),b(nprime)
      integer naxlen(7),incs(2),ifirst(2),ilast(2)
      data incs /1,1/
      data naxlen /7*0/
      logical anyf
c zero is the maximum pixel value that is considered to be equivalent to zero.
      if (mcov.gt.mcovmax) then
          write(0,*) 'Maximum covariance matrix offset:',mcovmax
          stop
      endif
      naxlen(1)=nx
      naxlen(2)=ny
      ntot=nx*ny
      if (ntot.lt.nprime) then
          istep=1
          i1=1
          nstep=ntot
      else
        istep=ntot/nprime 
         i1=(ntot-istep*(nprime-1))/2 + 1
         nstep=nprime
      endif
      k=0
      do i=i1,i1+(nstep-1)*istep,istep
          k=k+1
          ifirst(1)=max(1,min(nx-mcov,1+mod(i-1,nx)))
          ifirst(2)=max(1,min(ny-mcov,1+(i-1)/nx))
          ilast(1)=ifirst(1)+mcov
          ilast(2)=ifirst(2)+mcov
          istat=0
          call ftgsve(nunit,0,2,naxlen,ifirst,ilast,incs,0,zz(0,0,k),anyf,istat)
c          write(0,*) ((zz(ic,jc,k),ic=0,mcov),jc=0,mcov)
          if (istat.ne.0) then
              call ftrprt('stderr',istat)
              write(0,*) 'guesscov_sec: status=',istat
              stop
          endif
c check the patch does not include value=zero pixels
          do jc=0,mcov
              do ic=0,mcov
                  if (abs(zz(ic,jc,k)).le.zero) then
                      k=k-1
c     write(0,*) 'Skipping ',k,ifirst,ilast
                      goto 77
                  endif
              enddo
          enddo
c median-subtract the patch if large enough (at least 4x4)
          if (mcov.lt.3) goto 77
          bg=fmedian((mcov+1)**2,zz(0,0,k))
c          do jc=0,mcov
c              do ic=0,mcov
c                  zz(ic,jc,k)=zz(ic,jc,k)-bg
c              enddo
c          enddo
   77 enddo
      write(0,*) 'getcov_sec sampled',k,' patches.'
c first derive background and variance from the 0,0 pixels
      do i=1,k
          a(i)=zz(0,0,i)
      enddo
      call sort(k,a)
      bg=0.5*(a(k/2+1)+a((k+1)/2))
      iq1=nint(0.75+0.25*k)
      q1=a(iq1)
      qdev=bg-q1
c refine estimate by only looking at sky values between q1 and bg+(bg-q1)
      iq2=iq1+1
      do while (a(iq2).lt.2*bg-q1)
          iq2=iq2+1
      enddo
c      write(0,*) 'Est. bg, iq1, kmed, iq2',bg,iq1,0.5*(k/2+1+(k+1)/2),iq2
      bg=0.5*(a((iq1+iq2)/2+1)+a((iq1+iq2+1)/2))
      sig=qdev/0.675
      cov(0,0)=sig**2
      write(0,*) 'guess_cov using background ',bg
c now derive covariances
c only look at those pixel pairs which are on average <0, to filter out sources
c use medians of |z1+z2| and |z1-z2| (= quartiles of full dist without sources) 
c use (0,0) pixel as reference and correlate with (ic,jc) for +x,+y quadrant;
c use (0,mcov) pixel as reference and correlate with (ic,mcov+jc) for +x,-y quadrant.
c fill rest by symmetry 
      do ic=0,mcov
          do jc=-mcov,mcov
             kk=0
             do i=1,k
                if (jc.gt.0) then 
                   z00=zz(0,0,i)-bg
                   zij=zz(ic,jc,i)-bg
                else
                   z00=zz(0,mcov,i)-bg
                   zij=zz(ic,mcov+jc,i)-bg
                endif
                xplusy=z00+zij
                if (xplusy.le.0.) then
                   kk=kk+1
                   a(kk)=-xplusy
                   b(kk)=abs(zij-z00)
                endif
             enddo
             call sort(kk,a)
             call sort(kk,b)
             qdeva=0.5*(a(kk/2+1)+a((kk+1)/2))
             qdevb=0.5*(b(kk/2+1)+b((kk+1)/2))
             vara=(qdeva/0.675)**2
             varb=(qdevb/0.675)**2
             cov(ic,jc)=0.25*(vara-varb)
             cov(-ic,-jc)=cov(ic,jc)
c             write(*,*) ic,jc,cov(ic,jc),0.25*(vara+varb)
          enddo
      enddo
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
      sig=qdev/0.675
      cov(2,2)=sig**2
c calculate covariances from var(a-b)=var(a)+var(b)-2cov(a,b); use sqir to
c estimate variances. Median should be zero for a-b vectors.
      call sort(k,sky10)
      qdev=sky10(nint(0.75*k))-sky10(nint(0.25*k))
      vardiff=(qdev/0.675/2)**2
      write(0,*) 'var 10',vardiff
      cov(3,2)=sig**2-vardiff/2
      call sort(k,sky01)
      qdev=sky01(nint(0.75*k))-sky01(nint(0.25*k))
      vardiff=(qdev/0.675/2)**2
      write(0,*) 'var 01',vardiff
      cov(2,3)=sig**2-vardiff/2
      call sort(k,sky11)
      qdev=sky11(nint(0.75*k))-sky11(nint(0.25*k))
      vardiff=(qdev/0.675/2)**2
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
     :     '  Guessing image noise level from pixel value histogram'
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
      sig=qdev/0.675
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
