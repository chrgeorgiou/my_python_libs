c as imxshmap, but convolves with a shapelet PSF map that describes the residual deviation of the PSF from the desired Gaussian.
c input is a PSF map that is a difference between the measured one and the target Gaussian PSF, using a scale radius that is larger than the target Gaussian.
c the program convolves the input image with a kernel derived as follows:

c I --> O = I - DEL/G 
c
c where
c I is the input image inimage.fits
c O is the output image convolved2.fits
c DEL is the residual PSF map input to the programme via stdinput
c G is the target Gaussian, identified from the GPSFSIG keyword in the input image.
c (the division should be interpreted as a deconvolution)
c
c then if the PSF in the image is G + DEL, this convolution results in a PSF equal to G, to first order in DEL.
c
c Because DEL is expressed as a shapelet map with scale radius larger than the dispersion of G, this deconvolution can be done without fear of divergence.
c
c the programme uses the same deconvolution technique as imxshmap, dividing the image in blocks with appropriate edges. The main difference is the Gaussian factor, which here includes the deconvolution by the target gaussian.

      real kercoef(mmpoly2,msh-2)
      dimension sky(20000)
      parameter (m=512)
      dimension zz(m,m),out(m,m)
      dimension work1(m*m),work2(m*m)
      dimension naxlen(7)
      complex filt(m/2,m,mmpoly2)
      logical fitshavekey
      data naxlen /7*0./ 

c convolves an image by a spatially variable kernel
c expressed as a shapelet map. the image is split
c up into smaller sections and convolved piece by piece. 


      read(*,*) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(n+1)*(n+2)/2
         read(*,*) (kercoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo

      write(0,*) 'Read shapelet map of residual nongaussianity to clean'
      write(0,*) 'Order for shape, spatial variation; beta:',n,nfit,beta

c open image to be convolved
      call openfits('inimage.fits',51,nx,ny)
      if (fitshavekey(51,'GPSFSIG')) then
          gpsfsig=fitskeyr(51,'GPSFSIG')
          goto 2
      else
          open(10,file='gpsfsig.dat',status='old',err=1)
          read(10,*) gpsfsig
          write(0,*) '***Did not find GPSFSIG value in fits header.'
          write(0,*) 'Took value from file gpsfsig.dat'
          goto 2
      endif

    1 write(0,*) '***Did not find GPSFSIG value in fits header.'
      write(0,*) '***You can specify one in text file gpsfsig.dat'
      stop
    2 continue

      write(0,*) 'Input image wants Gaussian PSF of width ',gpsfsig
      if (gpsfsig.gt.beta*0.9) then
         write(0,*) '*** CANNOT DECONVOLVE BY GPSFSIG'
         write(0,*) '*** Shapelet scale radius too small'
         stop
      endif

c set width of taper edge
      ntap=nint(4*beta)

c calculate width of anti-alias border
      nmask=nint(beta*(5.+sqrt(1.*n)))+ntap
      if (nmask.gt.m/4) then
         write(0,*) 'Shapelet map too fat for Fourier convolution.'
         stop
      endif

c prepare the filters in Fourier domain
c include the deconvolution by gaussian with dispersion gpsfsig
      call shmapftdeconv(kercoef,beta,n,nfit,filt,m,gpsfsig)

      mpatch=m-2*nmask
c process the image in sections of mpatch x mpatch, padded to m x m
c outer border of nmask pix wide cannot be used.
      nxsub=int((nx-1)/mpatch)+1
      nysub=int((ny-1)/mpatch)+1
      write(0,*) 'Will process image in ',nxsub,' x ',nysub,' sections.'
      write(0,*) 'Will pad subimages from ',mpatch,mpatch,' to ',m,m
c open output image
c      call newfits('convolved.fits',61,nx,ny)
      naxlen(1)=nx
      naxlen(2)=ny

      call copyfits(51,61,'convolved2.fits',5)
      bg=get_bg_sec(51,nx*ny)

c go through the input image, read in a patch, mask a border around it,
c convolve, write out
      do iysub=1,nysub
         if (iysub.eq.1) then
            iy1=1-nmask
         else
            iy1=iy1+mpatch
         endif
         iy2=iy1+m-1
         if (iy2.gt.ny+nmask) then
            iy2=ny+nmask
            iy1=max(1,iy2-m+1)
         endif
         do ixsub=1,nxsub
            if (ixsub.eq.1) then
               ix1=1-nmask
            else
               ix1=ix1+mpatch
            endif
            ix2=ix1+m-1
            if (ix2.gt.nx+nmask) then
               ix2=nx+nmask
               ix1=max(1,ix2-m+1)
            endif
            mx=ix2-ix1+1
            my=iy2-iy1+1
            istat=igetsecpad(51,nx,ny,ix1,iy1,mx,my,zz,mx*my,bg)
            if (istat.ne.0.) then
               write(0,*) 'igetsec failed. code = ',istat
               stop
            endif
c convolve with shmap here

c taper the array to reduce aliasing
            call taperedge(zz,mx,my,ntap)            
            call shmapxpix(filt,m,nfit,zz,out,ix1,iy1,xmax,ymax,work1,work2)
c
c convolution done, result in out
c            write(0,*) 'Writing subsection ',ixsub,iysub,ix1,iy1
c pick out the part of out that is unaffected by aliasing (strip off border)
            call stripborder(out,zz,mx-2*nmask,my-2*nmask,m,nmask)
            call subbg(zz,(mx-2*nmask)*(my-2*nmask),-bg)
            call putsec(61,zz,mx-2*nmask,my-2*nmask,ix1+nmask,iy1+nmask,naxlen)
         enddo
         write(0,*) 'Finished strip ',iysub,' / ',nysub,
     :       ' = pix ',iy1+nmask,' - ',iy1+my-nmask-1
      enddo

c      write(0,*) 'Convolved image in convolved.fits'
c      end

      call closefits(51)

      call fitshistory(61,'PSF Gaussianization tweaked')
      call closefits(61)
      end



      subroutine shmapxpix(filt,m,nfit,cpix,cout,ix1,iy1,xm,ym,ftim,ctmp)

c convolve pixel array pix with shapelet map shmap
c Use the first npix x npix pixels of the mpix x mpix pixel array pix
c pix can be called as a real npix x xpix array - inside this routine it is
c treated as a complex npix/2 x npix array cpix
c     ftim and ctmp are temporary arrays that need to be dimensioned as
c     real npix*npix for the FT routine

      parameter(twopi=2*3.1415926535, mmax=2048)
      complex filt(m/2,m,mmpoly2)
      complex speq(mmax),speq0(mmax),cpix(m/2,m),cout(m/2,m)
      complex ftim(m/2,m),ctmp(m/2,m)

c Make Fourier transform of pixel array ftim (nondestructively!)
      do iy=1,m
         do ix=1,m/2
            ftim(ix,iy)=cpix(ix,iy)
         enddo
      enddo
      call rlft3(ftim,speq,m,m,1,1)
c     preserve the nyquist freq. component (it is not affected by the added filtered part)
      do iy=1,m
          speq(iy)=speq(iy)/(m*m/2.)
      enddo
c filter pixel array with each PSF map polynomial term in turn, add in result
      icpt=0
      do kcpt=0,nfit
         do iypow=0,kcpt
            ixpow=kcpt-iypow
            icpt=icpt+1
            do iy=1,m
               do ix=1,m/2
                  ctmp(ix,iy)=filt(ix,iy,icpt)*ftim(ix,iy)
               enddo
c     zero nyquist freq. cpt of the filtered part
               speq0(iy)=0
            enddo
c     transform back
            call rlft3(ctmp,speq0,m,m,1,-1)
c     add result into the output pixel array with proper weight
c     write out full array including padded parts - could trim this
c     write(0,*) 'Adding in part for x,y powers ',ixpow,iypow
            do iy=1,m
               if (iypow.eq.0) then
                  yfac=1
               else
                  yfac=(2*(iy+iy1-1.)/ym-1)**iypow
               endif
               do ix=1,m/2
                  z1=real(ctmp(ix,iy))
                  x1=2*ix-1 + ix1-1
                  z2=aimag(ctmp(ix,iy))
                  x2=2*ix + ix1-1
                  if (ixpow.eq.0) then
                     xfac1=1
                     xfac2=1
                  else
                     xfac1=(2*x1/xm-1)**ixpow
                     xfac2=(2*x2/xm-1)**ixpow
                  endif
c     build output image section
c     add first cpt to input image; rest on top of that
                  if (icpt.eq.1) then
                     cout(ix,iy)=cpix(ix,iy)+yfac*complex(xfac1*z1,xfac2*z2)
c                     cout(ix,iy)=            yfac*complex(xfac1*z1,xfac2*z2)
                  else
                     cout(ix,iy)=cout(ix,iy)+yfac*complex(xfac1*z1,xfac2*z2)
                  endif
               enddo
            enddo
c     complete the loops over polynomial coefs of the convolution kernel
         enddo
      enddo
      return
      end

      subroutine subbg(zz,m,bg)
c subtract the constant value bg from the zz array
      dimension zz(m)
      do i=1,m
          zz(i)=zz(i)-bg
      enddo
      return
      end

      subroutine taperedge(zz,mx,my,ntap)
c taper a strip of width ntap down to zero smoothly
      dimension zz(mx,my)
      do iy=1,ntap
         filt=(iy-1.)/ntap
         filt=(1-cos(3.1415926535*filt))/2.
         do ix=1,mx
            zz(ix,iy)=zz(ix,iy)*filt
         enddo
         do ix=1,mx
            zz(ix,my-iy+1)=zz(ix,my-iy+1)*filt
         enddo
      enddo
      do ix=1,ntap
         filt=(ix-1.)/ntap
         filt=(1-cos(3.1415926535*filt))/2.
         do iy=1,my
            zz(ix,iy)=zz(ix,iy)*filt
         enddo
         do iy=1,my
            zz(mx-ix+1,iy)=zz(mx-ix+1,iy)*filt
         enddo
      enddo
      return
      end


      subroutine stripborder(full,stripped,mx,my,m,nmask)
c write the mx,my sized block of the full array starting at
c nmask+1,nmask+1 to the stripped array
      dimension stripped(mx,my),full(m,m)
      do j=nmask+1,min(m,my+nmask)
         do i=nmask+1,min(m,mx+nmask)
            stripped(i-nmask,j-nmask)=full(i,j)
         enddo
      enddo
      return
      end


      subroutine shmapftdeconv(shmap,beta,n,nfit,filt,npix,gpsfsig)
      parameter(twopi=2*3.1415926535)
      complex filt(npix/2,npix,mmpoly2),tfilt
      real xtab(n+1),ytab(n+1),shmap(mmpoly2,msh-2)

c     construct the FT of a shapelet map
c     one filter per spatial polynomial term
c     deconvolved by a gaussian of width gpsfsig
c     expressed as a filter that corresponds to kernel 1-DEL/G
c     output is a complex array that can be multiplied with output of
c     rlft3 FT of an m x m array   

      icpt=0
      do kcpt=0,nfit
         do ixpow=kcpt,0,-1
            iypow=kcpt-ixpow
            icpt=icpt+1
c     construct shapelet function's FT at frequency fxnorm, fynorm
c     Hn(x/beta)exp(-0.5*x**2/beta**2) --> i^n Hn(f*beta) exp(-0.5*f**2
c     *beta**2)
            fmax2=(sqrt(n*1.)+5)**2
            do iy=1,npix
               fynorm=(mod(iy-1+npix/2,npix)-npix/2)*twopi*beta/npix
               call htab(fynorm,n,ytab)
               do ix=1,npix/2
                  fxnorm=(mod(ix-1+npix/2,npix)-npix/2)*twopi*beta/npix
                  call htab(fxnorm,n,xtab)
                  fnorm2=fxnorm*fxnorm+fynorm*fynorm
                  tfilt=0
                  if (fnorm2*(1-(gpsfsig/beta)**2).gt.fmax2) goto 7 
c Here insert the deconvolution by gpsfsig
                  efac=exp(-0.5*fnorm2*(1-(gpsfsig/beta)**2))
                  fac=efac/sqrt(3.1415926535)*beta*twopi/(npix*npix/2)
c
                  k=0
                  do nn=0,n
                     do j=0,nn
                        i=nn-j
                        k=k+1
                        if (j.ne.0) then
                           fac=fac/sqrt(j/real(i+1))
                        else 
                           if (i.ne.0) fac=fac/sqrt(real(i))
                        endif
                        tfilt=tfilt+(0.,1.)**nn*shmap(icpt,k)*fac*xtab(i+1)*ytab(j+1)
                     enddo
                     fac=fac/sqrt(2.)
                  enddo
c     filter FT of pixel array
c here write filter as -filt since we are convolving with 1-DEL/G
 7                filt(ix,iy,icpt)=-tfilt
               enddo
            enddo
c            write(0,*) 'Made FT of PSF map poly cpt ',icpt,ixpow,iypow
         enddo
      enddo
      write(0,*) 'Made FT''s of PSF map poly cpts 1 to',icpt
      return
      end


      
