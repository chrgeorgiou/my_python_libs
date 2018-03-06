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

c feed gapphot2 the x,y coordinates of the ends of the major and minor axes, in pixels, wrt center.
      do i=1,mmcat
 1       read (*,*,err=1,end=2) x,y,a,b,pa,id
         parad=pa*3.1415926535/180.
c major and minor axis ends in arcsec
         xasec= a*cos(parad)
         yasec= a*sin(parad)
         xbsec=-b*sin(parad)
         ybsec= b*cos(parad)
c and converted to pixels
         xa=xasec*cdinv11+yasec*cdinv12
         ya=xasec*cdinv21+yasec*cdinv22
         xb=xbsec*cdinv11+ybsec*cdinv12
         yb=xbsec*cdinv21+ybsec*cdinv22
         f=gapphot2(51,nx,ny,x,y,xa,ya,xb,yb,gpsfsig,bg,iflag,pixerr,fe)
         write(*,'(x,2f12.3,2f8.2,f7.1,g15.5,g13.3,g15.5,i6,i7)') x,y,a,b,pa,f,fe,bg,iflag,id
 1001    format('#    X [pix]     Y [pix]    A ["]   B ["] PA [deg]     FLX       FLX_ERR          BG          FLG   ID')
      enddo
 2    call ftclos()
      end

      function gapphot2(nu,nx,ny,x,y,xa,ya,xb,yb,g,bg,iflag,pixerr,ferr)
      parameter (maxpad=400, maxpad2=maxpad*maxpad)
c
c gaussian-wt aperture photometry. all in pixels.
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
c need pixel cutout that reaches to radius 2*sqrt(2)*r
      ir=max(15,int(2.8*r+1))
      np=min(maxpad,2*ir+1)
c check whether PSF is smaller than requested aperture; otherwise flag and exit
      if (g.ge.min(a,b)) then
         iflag=900
         gapphot2=0
         ferr=-1
         return
      endif
c first check if source is in the image or not
      if (ix-ir.le.0 .or. iy-ir.le.0 .or. ix+ir.gt.nx .or. iy+ir.gt.ny) then
         iflag=100
         gapphot2=0
         ferr=-1
         return
      endif
c get postage stamp.
c check whether patch was OK; otherwise flag and exit
      isok=igetsec(nu,nx,ny,ix-ir,iy-ir,np,np,zz,maxpad2)
      if (isok.ne.0) then
         iflag=100
         gapphot2=0
         ferr=-1
         return
      endif
      call dogapphot2(zz,np,x-ix+ir+1,y-iy+ir+1,r,a,b,xa,ya,xb,yb,g,f,bg,ferr)
      gapphot2=f
      ferr=ferr*pixerr
      return
      end


      subroutine dogapphot2(zz,n,x,y,r,a,b,xa,ya,xb,yb,g,f,bg,ferr)
      dimension zz(n,n)
      double precision sum
c fit background first (only if postage stamp is big enough)
      bg=0.
c      if (r/n.lt.0.3) call fitbg(zz,n,n,n,x,y,r*2,bg)
c photometry correction factors, PSF-corrected weight fn radius
      ra2=a*a-g*g
      rb2=b*b-g*g
      pcorr=a*b/sqrt(ra2*rb2)
c do the weighted aperture photometry
      sum=0
      sumw2=0
      do j=1,n
         do i=1,n
            dx=i-x
            dy=j-y
            w=exp(-0.5*((dx*xa+dy*ya)**2/(a*a*ra2) + ((dx*xb+dy*yb)**2/(b*b*rb2))))
            if (w.gt.1.e-5) then
               sum=sum+(zz(i,j)-bg)*w
               sumw2=sumw2+w*w
            endif
         enddo
      enddo
      f=sum*pcorr
      ferr=sqrt(sumw2)*pcorr
      return
      end

