c read fits file and weight image.
c for all pixels with zero weight, determine whether pixel value falls outside given limits. If so, replace with the limits.


      parameter (m=5000)
      real z(m*m),w(m*m),naxlen(7)
      logical anyf

c get bandwidth limits for the flagged pixels, and replacement values
c everything below zmin becomes zlo, everything above zmax becomes zhi
      read(*,*) zmin,zmax,zlo,zhi

      call openfits('inimage.fits',51,nx,ny)
      call openfits('weights.fits',61,nxw,nyw)
      if (nx.ne.nxw .or. ny.ne.nyw) then
          write(0,*) 'inimage.fits and weights.fits have different sizes!'
          stop
      endif

      istat=0
      call ftg2de(51,0,0.,nx,nx,ny,z,anyf,istat)
      if (istat.ne.0)  call ftrprt('stderr',istat)

      istat=0
      call ftg2de(61,0,0.,nx,nx,ny,w,anyf,istat)
      if (istat.ne.0)  call ftrprt('stderr',istat)

      call ftclos(51,istat)
      call ftclos(61,istat)

      do i=1,nx*ny
          if (w(i).eq.0.) then
              if (z(i).lt.zmin) then 
                  z(i)=zlo
              elseif (z(i).gt.zmax) then
                  z(i)=zhi
              endif
          endif
      enddo

      naxlen(1)=nx
      naxlen(2)=ny

      call copyfits_nodata('inimage.fits',52,'clipped.fits')
      call ftp2de(52,0,nx,nx,ny,z,istat)
c      call putsec(52,z,nx,ny,1,1,naxlen)
      istat=0
      call ftclos(52,istat)
      end
