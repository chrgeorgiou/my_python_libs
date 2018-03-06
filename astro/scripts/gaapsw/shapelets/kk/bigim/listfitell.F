c take list of coordinates of sources in a fits image and determine crude best-fit 2D gaussian.

c open fits image, determine noise properties (ideally noise cov matrix)


      logical fitshavekey


      call openfits('inimage.fits',51,nx,ny)
      pixerr=guess_errors_sec(51,nx*ny)

c determine PSF radius from header; set to zero if not specified
      
      if (fitshavekey(51,'GPSFSIG')) then
          gpsfsig=fitskeyr(51,'GPSFSIG')
          write(0,*) 'Found Gaussian PSF radius in header: value is ',gpsfsig
      else
         gpsfsig=0
         write(0,*) '***Did not find GPSFSIG value in fits header.'
         write(0,*) '***Will not apply PSF correction to models.'
      endif
c get CD matrix from header and convert to arcsec/pix.
c CD matrix gives transformation pix --> degrees
      if (fitshavekey(51,'CD1_1')) then
         cd11=fitskeyr(51,'CD1_1')
         cd12=fitskeyr(51,'CD1_2')
         cd21=fitskeyr(51,'CD2_1')
         cd22=fitskeyr(51,'CD2_2')
         cd11=cd11*3600
         cd12=cd12*3600
         cd21=cd21*3600
         cd22=cd22*3600
         write(0,*) 'CD: ',cd11, cd12, cd21, cd22
      else
         cd11=-0.2
         cd12=0
         cd21=0
         cd22= 0.2
         write(0,*) '*** NO CD MATRIX: USE 0.2 ARCSEC/PIX, N UP, E LEFT'
      endif
c calculate inverse CD matrix for transformation from degrees to pixels
      cddet=cd11*cd22-cd12*cd21
      cdinv11= cd22/cddet
      cdinv12=-cd12/cddet
      cdinv21=-cd21/cddet
      cdinv22= cd11/cddet
      pixsize=sqrt(abs(cddet))
      write(0,*) 'Pixel scale (arcsec/pix) = ',pixsize
c read in the object list, one by one, and do the Gaussian-model fitting
c format: x,y coordinate (pixels)
c
      do i=1,mmcat
 1       read (*,*,err=1,end=2) x,y
         ix=nint(x)
         iy=nint(y)
         call fit2dgaus(51,nx,ny,x,y,gpsfsig,f,a,b,pa,bg,iflag)
         write(*,*) x,y,a,b,pa,bg,iflag
      enddo
 2    call ftclos(51)
      end

      subroutine fit2dgaus(nu,nx,ny,x,y,g,f,a,b,pa,bg,iflag)
      parameter (maxpad=300, maxpad2=maxpad*maxpad)
c gaussian-wt aperture photometry. all in pixels.
c x,y ctr of source
c xa,ya = x,y offset to 1-sigma on major axis; xb,yb ditto on minor axis
c g = gaussian psf radius
c
      real zz(maxpad2)
      iflag=0
      ix=nint(x)
      iy=nint(y)
      ir=100
      np=min(maxpad,2*ir+1)
c first check if source is in the image or not
      if (ix-ir.le.0 .or. iy-ir.le.0 .or. ix+ir.gt.nx .or. iy+ir.gt.ny) then
         iflag=100
         f=-1
         return
      endif
c get postage stamp.
c check whether patch was OK; otherwise flag and exit
      isok=igetsec(nu,nx,ny,ix-ir,iy-ir,np,np,zz,maxpad2)
      if (isok.ne.0) then
         iflag=100
         f=-1
         return
      endif
c crude background estimate from radius beyond 50 pixels
      call fitbg(zz,np,np,np,ir+1.,ir+1.,50.,bg)
      ch2min=1e32
      gmin=0
      fmin=-1e32
c first find best-fit round gaussian on a coarse grid of radii, starting at PSF size
      do ig=0,10
         s=sqrt(g*g+ig*ig)
         ch2=fit1dgaus(zz,np,np,x-ix+ir+1,y-iy+ir+1,s,bg,f)
         if (ch2.lt.ch2min) then
            ch2min=ch2
            gmin=ig
            fmin=f
         endif
         write(0,*) ig,s,f,ch2
      enddo
c now refine this model: set the pixels to be fitted and do a proper fit to a 4-parameter gaussian (keep ctr).

c set best-fit results f,a,b,pa,bg,iflag      
      return
      end



      function fit1dgaus(zz,m,n,x,y,s,bg,f)
      dimension zz(m,m)
      double precision summ2,summd,sumd2
      summ2=0
      summd=0
      sumd2=0
      s2=s*s
      r2max=9*s2+4
      do iy=1,n
         dy=y-iy
         do ix=1,n
            dx=x-ix
            r2=dx*dx+dy*dy
            if (r2.lt.r2max) then
               gmod=exp(-0.5*r2/s2)
               dat=zz(ix,iy)-bg
               summd=summd+dat*gmod
               sumd2=sumd2+dat*dat
               summ2=summ2+gmod*gmod
            endif
         enddo
      enddo
      f=summd/summ2
      ch2=f*(f*summ2-2*summd)+sumd2
      fit1dgaus=ch2/r2max
      return
      end
