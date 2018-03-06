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

      parameter (pi=3.1415926535)
      logical fitshavekey,havecovmap,havewtim,readpsfmap
      real noiseacfcoef(mmpoly2,msh-2),shacf(msh),cov(-10:10,-10:10)
      character*160 keracfname

      call openfits('inimage.fits',51,nx,ny)

      inquire(file='weights.fits',exist=havewtim)
      if (havewtim) then
         call openfits('weights.fits',61,nxw,nyw)
         havewtim = (nxw.eq.nx .and. nyw.eq.ny) 
      endif
      pixerr=guess_errors_sec(51,nx*ny)

c get covariance matrix of the noise.
c first check if the ACF of the convolution kernel exists: if it exists use that in combination with either
c     variance of convolved image as guessed in pixerr or
c     variance of pre-convolved image as given in the weight image 
      read(*,'(a160)') keracfname
      havecovmap=readpsfmap(keracfname,noiseacfcoef,nacf,betaacf,nfitacf,xmaxacf,ymaxacf)
      if (havecovmap) then
         call interpolpsf(0.5*nx,0.5*ny,xmaxacf,ymaxacf,noiseacfcoef,nacf,betaacf,nfitacf,shacf)
         do j=-10,10
            do i=-10,10
               cov(i,j)=shapeletpix(real(i),real(j),shacf)
            enddo
         enddo
         call put2('keracf.fits',21,21,cov,21,21)
      else
c if there is no kernel ACF, then estimate the covariance from the image itself.
         call guesscov_sec (51,nx,ny,10,cov,pixerr*0.001)
         call fitbg(cov,21,21,21,11.,11.,8.,covbg)
         do j=-10,10
            do i=-10,10
               cov(i,j)=cov(i,j)-covbg
            enddo
         enddo
         write(0,*) 'Subtracted a large-scale covariance background of',covbg
         write(0,*) 'Guessed covariance map from noise in image.'
         write(0,*) 'Peak cov, pix var:',cov(0,0),pixerr**2
         thresh=max(cov(0,0),pixerr**2)*0.04
         write(0,*) 'Zeroing all elements with |cov| less than ',thresh
         do j=-10,10
            do i=-10,10
               if (abs(cov(i,j)).lt.thresh) cov(i,j)=0
            enddo
         enddo
         call put2('cov.fits',21,21,cov,21,21)
      endif

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
         parad=pa*pi/180.
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
c get the gaussian aperture photometry measurement
         f=gapphot2(51,nx,ny,x,y,xa,ya,xb,yb,gpsfsig,bg,iflag,pixerr,fe)
c now for the error bar:
c   case where there was a shapelet map of the convolution kernel ACF:
         if (havecovmap) then
c first get the covariance matrix = normalized kernel ACF * original pixel variance
            call interpolpsf(x,y,xmaxacf,ymaxacf,noiseacfcoef,nacf,betaacf,nfitacf,shacf)
            acfpk=shapeletpix(0.,0.,shacf)
c         write(0,*) 'Normalized ACF peak =',acfpk
c         write(0,*) 'ACF integral ',seriessum(shacf)
c         write(0,*) (shacf(k),k=1,(nacf+1)*(nacf+2)/2+2)
c get pre-convolution pixel variance from original weight image if it exists, otherwise use the one estimated from the image itself
            if (havewtim) then
               pixwt=getonepix(61,nint(x),nint(y),nx,ny)
               if (pixwt.gt.0.) then
                  pixvar=1./pixwt
               else
                  pixvar=-1
                  fsig=-1
               endif
            else
               pixvar=pixerr**2/acfpk
            endif
c         write(0,*) 'x,y,pixvar =',x,y,pixvar
            if (pixvar.gt.0.) then
               apix=a/pixsize
               bpix=b/pixsize
               fsig=-1
               if (apix.gt.gpsfsig .and. bpix.gt.gpsfsig) then
c convolution of gaussian(sqrt(a^2-gpsfsig^2), ditto[b])  with itself is
c gaussian with axes sqrt(2) times as long, and prefactor (2 pi a b)^2 / (2 pi 2 a b) = pi a b 
                  afilt=sqrt(2*(apix**2-gpsfsig**2))
                  bfilt=sqrt(2*(bpix**2-gpsfsig**2))
                  var=pixvar*shintellgaus(shacf,afilt,bfilt,parad)*(pi*afilt*bfilt/2)
c gaap flux has correction factor: apply square to variance
                  var=var*(apix*bpix)**2/((apix**2-gpsfsig**2)*(bpix**2-gpsfsig**2))
c               write(0,*) 'var=',var
                  if (var.gt.0.) then
                     fsig=sqrt(var)
                  endif
               endif
            endif
         else
c case where there is no shapelet noise covariance matrix: use the estimated one (note this is constant over image)
c and calculate integral of COV(x) * pi a b exp (-0.25 (X/a)^2+(Y/b)^2) with X and Y rotated coordinates
            apix=a/pixsize
            bpix=b/pixsize
            afilt2=apix**2-gpsfsig**2
            bfilt2=bpix**2-gpsfsig**2
            sum=0
            c=cos(parad)
            s=sin(parad)
            do iy=-10,10
               do ix=-10,10
                  xr= ix*c + iy*s
                  yr=-ix*s + iy*c
                  em2=xr**2/afilt2+yr**2/bfilt2
                  if (em2.le.40.) sum=sum+exp(-0.25*em2)*cov(ix,iy)
c                  write(0,*) ix,iy,xr,yr,em2,cov(ix,iy)
               enddo
            enddo
            var=sum*pi*sqrt(afilt2*bfilt2) * (apix*bpix)**2/(afilt2*bfilt2)
c            write(0,*) 'sum,var',sum,var
            if (var.gt.0.) then
               fsig=sqrt(var)
            else
               fsig=-1
            endif
         endif
         write(*,'(x,2f12.3,2f8.2,f7.1,g15.5,g13.3,g15.5,i6,i7,g15.5)') x,y,a,b,pa,f,fe,bg,iflag,id,fsig
 1001    format('#    X [pix]     Y [pix]    A ["]   B ["] PA [deg]     FLX       FLX_ERR          BG          FLG   ID')
      enddo
 2    call ftclos(51)
      call ftclos(61)
      end



      function shintellgaus(sh,a,b,parad)
c integral of shapelet sh times elliptical Gaussian with major axis a, minor axis b, position angle parad (in radian from X axis)
c sh contains value of shapelet scale radius beta
c
c integral Hm(x/beta) Hn(y/beta) exp(-0.5(x^2+y^2)/beta^2) /sqrt[pi 2^(m+n)m!n!]  *   exp(-0.5 [X^2 + Y^2])  .dx dy
c where X and Y are rotated and stretched coordinates st position angle of major axis is parad:
c  aX=x cos(pa) + y sin(pa) ; bY = -x sin(pa) + y cos(pa)
c
c multiplying the gaussians effectively changes a^2 to A^2=1/(1/a^2+1/beta^2), ditto b^2 --> B^2, leaves PA the same.
c so absorb the exp(-0.5 r^2/beta^2) into the elliptical gaussian, and use A and B from now on.
c transform to X,Y gives 
c integral A dX B dY Hm(u1 X + u2 Y) Hn (v1 x + v2 Y) exp(-0.5[X^2+y^2])
c with 
c u1= A/beta cos PA
c u2=-B/beta sin PA
c v1= A/beta sin PA
c v2= B/beta cos PA
c
c can derive a recurrence relation
c I_m,n = (m-1) (4u1^2+4u2^2-2) I_m-2,n + 4n (u1 v1 +u2 v2) I_m-1,n-1
c and similarly m <--> n, u <--> v

      parameter (pi=3.1415926535,nmn=(mm+1)**2)
      real sh(msh)
      real imn(0:mm,0:mm)
      integer*8 mfac,nfac
      data imn /nmn*0/

      beta=sh(1)
      n=sh(2)

      biga2=1./(1./a**2+1./beta**2)
      bigb2=1./(1./b**2+1./beta**2)

c coefficients of the recurrence relation
      coefmmin2n=4*biga2/beta**2*cos(parad)**2+4*bigb2/beta**2*sin(parad)**2-2
      coefmmin1nmin1=2*(biga2-bigb2)/beta**2*sin(2*parad)
      coefmnmin2=4*biga2/beta**2*sin(parad)**2+4*bigb2/beta**2*cos(parad)**2-2

c seed the recurrence. I_00 is simple 2D elliptical Gaussian integral
      imn(0,0)=2*pi*sqrt(biga2*bigb2)
c work up the m=0 axis
      do in=2,n,2
         imn(0,in)=(in-1)*coefmnmin2*imn(0,in-2)
      enddo
c and fill in the remainder of the table:
c first do each row from left to right, then work up to row above
      do in=0,n
         do im=1,n-in
            if (im.gt.1 .and. in.gt.0) then
               imn(im,in)=(im-1)*coefmmin2n*imn(im-2,in) + in*coefmmin1nmin1*imn(im-1,in-1)
            elseif (im.gt.1) then
               imn(im,in)=(im-1)*coefmmin2n*imn(im-2,in)
            elseif (in.gt.0) then
               imn(im,in)=in*coefmmin1nmin1*imn(im-1,in-1)
            endif
         enddo
      enddo
c      write(0,*) 'pre-norm', ((im,in,imn(im,in),im=0,n-in),in=0,n)
c now modify the coefficients to include normalization of the shapelets:  divide by sqrt(pi m!n! 2^(m+n)) beta
      nfac=1
      do in=0,n
         mfac=1
         do im=0,n-in
            imn(im,in)=imn(im,in)/sqrt(pi*mfac*nfac)/beta
c            write(0,*) im,in,mfac,nfac
            mfac=mfac*2*(im+1)
         enddo
         nfac=nfac*2*(in+1)
      enddo
c      write(0,*) 'post-norm',((im,in,imn(im,in),im=0,n-in),in=0,n)
c now sum shapelet term by shapelet term.
      sum=0
      ksh=2
      do nn=0,n,2
         do j=0,nn
            i=nn-j
            ksh=ksh+1
            sum=sum+imn(i,j)*sh(ksh)
         enddo
         ksh=ksh+nn+2
      enddo
      shintellgaus=sum
      return
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

