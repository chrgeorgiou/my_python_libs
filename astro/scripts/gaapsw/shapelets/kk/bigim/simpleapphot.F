c perform list-driven gaussian aperture photometry on a fits image
c for each source specify x,y pixel ctr, and Gaussian radius
c
c if the header of the image contains a GPSFSIG keyword then use its
c value as Gaussian PSF width, and interpret the apertures are pre-seeing
c
c if the GPSFSIG keyword is absent then interpret its value as zero and
c simply calculate apertures given
c
c     The covariance matrix of the noise is incorporated in the error
c     estimate of the photometry. This covariance matrix can either be
c     supplied as a shapelet map (eg by running the convolution kernel
c     shapelet map through keracfmap), or otherwise it is estimated from
c     the noise properties of the image itself. In the latter case the
c     covariance will be assumed to be constant over the whole image.
c

c open fits image, determine noise properties (ideally noise cov matrix)

      parameter (pi=3.1415926535)
      logical fitshavekey

      write(0,1000)
 1000 format ('SIMPLE APERTURE PHOT -------------------------')

      call openfits('inimage.fits',51,nx,ny)
      pixerr=guess_errors_sec(51,nx*ny)

c get CD matrix from header and convert to arcsec/pix.
c CD matrix gives transformation pix --> degrees
c if there is no CD matrix then read the pixel scale form CDELT1,2
c failing that, use 0.2 arcsec per pixel as default.
c this is used to transform the size of the aperture and its orientation into pixels (not for absolute astrom).
c      cd11=fitskeyr(51,'CDELT1')
c      cd22=fitskeyr(51,'CDELT2')
c      write(0,*) cd11,cd22
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
      elseif (fitshavekey(51,'CDELT1')) then
         cd11=fitskeyr(51,'CDELT1')
         cd12=0
         cd21=0
         cd22=fitskeyr(51,'CDELT2')
         cd11=cd11*3600
         cd22=cd22*3600
         write(0,*) 'CDELT: ',cd11,cd22
      else
         cd11=-0.2
         cd12=0
         cd21=0
         cd22= 0.2
         write(0,*) '*** NO CD MATRIX FOUND - ASSUME 0.2 ARCSEC/PIX'
      endif
c calculate inverse CD matrix for transformation from degrees to pixels
      cddet=cd11*cd22-cd12*cd21
      cdinv11= cd22/cddet
      cdinv12=-cd12/cddet
      cdinv21=-cd21/cddet
      cdinv22= cd11/cddet
      pixsize=sqrt(abs(cddet))
      write(0,*) 'Pixel scale (arcsec/pix) = ',pixsize
c read in the object list, one by one, and do the Gaussian-aperture photometry
c format: x,y coordinate (pixels), major, minor axis, major axis position angle (deg)
c
c x,y coordinates are in pixels; the aperture size is in arcsec and will
c be converted to pixels with the CD matrix in the FITS header.
c
c write header line
      write(*,1001)

c feed gapphot4 the x,y coordinates of the ends of the major and minor axes, in pixels, wrt center.
      do i=1,mmcat
 1       read (*,*,err=1,end=2) x,y,a,b,pa,id
         parad=pa*pi/180.
c pa is in mathematical convention, from +ve X (PA=0) to +ve Y (PA=90) axis.
c major and minor axis ends in arcsec along image X and Y axes
         xasec= a*cos(parad)
         yasec= a*sin(parad)
         xbsec=-b*sin(parad)
         ybsec= b*cos(parad)
c and converted to pixels. Change sign of x coordinate since RA runs from right to left.
         xa=-xasec*cdinv11+yasec*cdinv12
         ya=-xasec*cdinv21+yasec*cdinv22
         xb=-xbsec*cdinv11+ybsec*cdinv12
         yb=-xbsec*cdinv21+ybsec*cdinv22
c         write(0,*) xa,ya,xb,yb
c
c get the gaussian aperture photometry measurement
         f=simpleapphot(51,nx,ny,x,y,xa,ya,xb,yb,bg,iflag)
c fe is estimate of the error assuming uncorrelated noise. Negative fe indicates no flux could be obtained
         write(*,'(x,2f12.3,2f8.2,f7.1,g15.5,g13.3,g15.5,i6,i7,g15.5)') x,y,a,b,pa,f,fsig,bg,iflag,id,fe
 1001    format('#    X [pix]     Y [pix]    A ["]   B ["] PA [deg]  ',
     :        'FLX       FLX_ERR          BG          FLG   ID   ',
     :        'ERR_UNCORR')
      enddo
 2    write(0,*) 'Number of sources processed by simpleapphot: ',i-1
      write(0,*)
      call closefits(51)
      end


      function simpleapphot(nu,nx,ny,x,y,xa,ya,xb,yb,bg,iflag)
      parameter (maxpad=400, maxpad2=maxpad*maxpad)
c
c top-hat elliptical aperture photometry. all in pixels.
c x,y ctr of source
c xa,ya = x,y offset to 1-sigma on major axis; xb,yb ditto on minor axis
c g = gaussian psf radius
c maxpad is the maximum size of the image cutout
c
      real zz(maxpad2)
      iflag=0
      ix=nint(x)
      iy=nint(y)
      a=sqrt(xa*xa+ya*ya)
      b=sqrt(xb*xb+yb*yb)
      r=max(a,b)
c need pixel cutout that reaches to radius r
      ir=max(15,int(r+1))
      np=min(maxpad,2*ir+1)
c first check if source is in the image or not
      if (ix-ir.le.0 .or. iy-ir.le.0 .or. ix+ir.gt.nx .or. iy+ir.gt.ny) then
         iflag=100
         simpleapphot=0
         ferr=-1
         return
      endif
c get postage stamp.
c check whether patch was OK; otherwise flag and exit
      isok=igetsec(nu,nx,ny,ix-ir,iy-ir,np,np,zz,maxpad2)
      if (isok.ne.0) then
         iflag=100
         simpleapphot=0
         ferr=-1
         return
      endif
      call dosimpleapphot(zz,np,x-ix+ir+1,y-iy+ir+1,r,a,b,xa,ya,xb,yb,g,f,bg)
      simpleapphot=f
      return
      end


      subroutine dosimpleapphot(zz,n,x,y,r,a,b,xa,ya,xb,yb,g,f,bg)
      dimension zz(n,n)
      double precision sum
c fit background first (only if postage stamp is big enough)
      bg=0.
c      if (r/n.lt.0.3) call fitbg(zz,n,n,n,x,y,r*2,bg)
c do the aperture photometry
      sum=0
      do j=1,n
         do i=1,n
            dx=i-x
            dy=j-y
            w=(dx*xa+dy*ya)**2/a**4 + (dx*xb+dy*yb)**2/b**4
            if (w.le.1.) then
               sum=sum+(zz(i,j)-bg)
            endif
         enddo
      enddo
      f=sum
      return
      end

