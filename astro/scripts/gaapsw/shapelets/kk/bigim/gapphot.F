c perform list-driven gaussian aperture photometry on a fits image
c for each source specify x,y pixel ctr, and Gaussian radius
c
c if the header of the image contains a GPSFSIG keyword then use its
c value as Gaussian PSF width, and interpret the apertures are pre-seeing
c
c if the GPSFSIG keyword is absent then interpret its value as zero and
c simply calculate apertures given
c

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
         write(0,*) '***Will not apply PSF correction to photometry.'
      endif

c read in the object list, one by one, and do the photometry
      do i=1,mmcat
 1       read (*,*,err=1,end=2) x,y,r
         f=gapphot(51,nx,ny,x,y,r,gpsfsig,bg,iflag,pixerr,fe)
         write(*,*) x,y,r,f,fe,bg,iflag
      enddo
 2    call ftclos()
      end

      function gapphot(nu,nx,ny,x,y,r,g,bg,iflag,pixerr,ferr)
      parameter (maxpad=400, maxpad2=maxpad*maxpad)

      real zz(maxpad2)
      iflag=0
      ix=nint(x)
      iy=nint(y)
      ir=max(15,int(5*r+1))
      np=min(maxpad,2*ir+1)
c check whether PSF is smaller than requested aperture; otherwise flag and use PSF size instead
      if (g.ge.r) then
         iflag=900
         g=r
      endif
c first check if source is in the image or not
      if (ix-ir.le.0 .or. iy-ir.le.0 .or. ix+ir.gt.nx .or. iy+ir.gt.ny) then
         iflag=100
         gapphot=0
         ferr=-1
         return
      endif
c get postage stamp.
c check whether patch was OK; otherwise flag and exit
      isok=igetsec(nu,nx,ny,ix-ir,iy-ir,np,np,zz,maxpad2)
      if (isok.ne.0) then
         iflag=100
         gapphot=0
         ferr=-1
         return
      endif
      call dogapphot(zz,np,x-ix+ir+1,y-iy+ir+1,r,g,f,bg,ferr)
      gapphot=f
      ferr=ferr*pixerr
      return
      end


      subroutine dogapphot(zz,n,x,y,r,g,f,bg,ferr)
      dimension zz(n,n)
      double precision sum
c fit background first
      call fitbg(zz,n,n,n,x,y,r*5,bg)
c photometry correction factors, PSF-corrected weight fn radius
      rw2=r*r-g*g
      pcorr=r*r/rw2
c do the weighted aperture photometry
      sum=0
      sumw2=0
      do j=1,n
         do i=1,n
            w=exp(-0.5*((i-x)**2+(j-y)**2)/rw2)
            sum=sum+(zz(i,j)-bg)*w
            sumw2=sumw2+w*w
         enddo
      enddo
      f=sum*pcorr
      ferr=sqrt(sumw2)*pcorr
      return
      end

