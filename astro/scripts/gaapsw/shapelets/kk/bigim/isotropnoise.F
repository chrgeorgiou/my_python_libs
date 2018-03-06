      real kercoef(mmpoly2,msh-2)
      parameter (m=512)
      dimension zz(m,m),zz2(m,m),out(m,m),ftacfiso(m/2,m)
      dimension work1(m*m),work2(m*m)
      dimension znoise(30000)
      dimension naxlen(7)
      complex filt(m/2,m,mmpoly2),ftacf(m/2,m),speq(m)
      data naxlen /7*0./ 

c based on imxshmap
c
c reads an image and a kernel map, and makes two new images:
c
c 1. convolved.fits is inimage.fits convolved with the kernel map
c
c 2. isonoise.fits is appropriately convolved noise to add to
c    convolved.fits so that the resulting covariance matrix is
c    isotropic.
c


      read(*,*) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(n+1)*(n+2)/2
         read(*,*) (kercoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo

      write(0,*) 'Read kernel shapelet map'
      write(0,*) 'Order for shape, spatial variation; beta:',n,nfit,beta

c set width of taper edge
      ntap=nint(4*beta)

c calculate width of anti-alias border
      nmask=nint(beta*(5.+sqrt(1.*n)))+ntap
      if (nmask.gt.m/4) then
         write(0,*) 'Shapelet map too fat for Fourier convolution.'
         stop
      endif

c prepare the filters in Fourier domain
      call shmapft(kercoef,beta,n,nfit,filt,m)

c open image to be convolved
      call openfits('inimage.fits',51,nx,ny)
      signoise=guess_errors_sec(51,nx*ny)
      bg=get_bg_sec(51,nx*ny)

      write(0,*) 'Calculate the kernel power spectrum and isotropize.'
      call kerftmax(filt,m,nfit,xmax,ymax,nx,ny,ftacfiso)

c max size that can be fitted into the m x m array that will be FT'd,
c accounting for borders
      mpatch=m-2*nmask

c process the image in sections of mpatch x mpatch, padded to m x m
c outer border of nmask pix wide cannot be used.
      nxsub=int((nx-1)/mpatch)+1
      nysub=int((ny-1)/mpatch)+1
      write(0,*) 'Will process image in ',nxsub,' x ',nysub,' sections.'
      write(0,*) 'Will pad subimages from ',mpatch,mpatch,' to ',m,m
c open output images
      naxlen(1)=nx
      naxlen(2)=ny
      call copyfits(51,61,'convolved.fits',5)
      call copyfits(51,71,'isonoise.fits',5)
c create noise image of same size, populate with uniform gaussian noise
c of amplitude equal to original image
      call copyfits(51,52,'tmpnoise.fits',5)
c     set random number seed based on the fractional part of the kernel
c     00 coef and on signoise
      c=kercoef(1,1)
      idum=int(10000*(c-int(c)))
      c=signoise
      idum=idum+10000*int(10000*(c-int(c)))
      idum=-idum
      do j=1,ny
          do i=1,nx
              znoise(i)=gasdev(idum)*signoise
          enddo
          call putsec(52,znoise,nx,1,1,j,naxlen)
      enddo
      write(0,*) 'White noise image tmpnoise.fits created'

c construct the target ACF of the noise, in Fourier space: maximize
c |kernel**2| across the image this k-mode by k-mode and then isotropize


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
c get image section and subtract background. Pad with zeroes if needed.
            istat=igetsecpad(51,nx,ny,ix1,iy1,mx,my,zz,mx*my,bg)
            if (istat.ne.0.) then
               write(0,*) 'igetsecpad failed. code = ',istat
               stop
            endif
c convolve inimage section with with shmap here
c taper the array to reduce aliasing
            call taperedge(zz,mx,my,ntap)               
            call shmapxpix(filt,m,nfit,zz,out,ix1,iy1,xmax,ymax,work1,work2)
c
c convolution done, result in out
c            write(0,*) 'Writing subsection ',ixsub,iysub,ix1,iy1
c pick out the part of out that is unaffected by aliasing (strip off border)
c put the background back in
            call stripborder(out,zz,mx-2*nmask,my-2*nmask,m,nmask)
            call subbg(zz,(mx-2*nmask)*(my-2*nmask),-bg)
            call putsec(61,zz,mx-2*nmask,my-2*nmask,ix1+nmask,iy1+nmask,naxlen)


c now repeat for the noise image. no need to b/g subtract here
            istat=igetsecpad(52,nx,ny,ix1,iy1,mx,my,zz2,mx*my,0.)
            call taperedge(zz2,mx,my,ntap)
c construct the kernel's FT at the center of this image section
            xc=0.5*(ix1+ix2)
            yc=0.5*(iy1+iy2)
            xnorm=2*xc/xmax-1
            ynorm=2*yc/ymax-1
            do l=1,m
                do k=1,m/2
                    ftacf(k,l)=0
                enddo
            enddo
            icpt=0
            do kcpt=0,nfit
                do iypow=0,kcpt
                    ixpow=kcpt-iypow
                    icpt=icpt+1
                    fac=xnorm**ixpow * ynorm**iypow
                    do l=1,m
                        do k=1,m/2
                            ftacf(k,l)=ftacf(k,l)+fac*filt(k,l,icpt)
                        enddo
                    enddo
                enddo
            enddo
c turn kernel FT into FT of ACF (ie square modulus of FT) and difference
c from target. Take sqrt to yield convolution kernel with that ACF.
            do l=1,m
                do k=1,m/2
                    ftacf(k,l)=sqrt(max(0.,ftacfiso(k,l)-abs(ftacf(k,l))**2))
                enddo
            enddo
c convolve white noise with this kernel, and write out to isotropization
c noise image.  can use the same convolution routine but now give it
c zero nfit since assume kernel is constant across img section
            call shmapxpix(ftacf,m,0,zz2,out,ix1,iy1,xmax,ymax,work1,work2)
            call stripborder(out,zz2,mx-2*nmask,my-2*nmask,m,nmask)
c            call secplussec(zz2,zz,(mx-2*nmask)*(my-2*nmask))
            call putsec(71,zz2,mx-2*nmask,my-2*nmask,ix1+nmask,iy1+nmask,naxlen)
         enddo
         write(0,*) 'Finished strip ',iysub,' / ',nysub,
     :       ' = pix ',iy1+nmask,' - ',iy1+my-nmask-1
      enddo

c      write(0,*) 'Convolved image in convolved.fits'
c      end

      call closefits(51)
      call closefits(52)
c delete the temporary noise image again
c      open(52,file='tmpnoise.fits',status='old')
c      close(52,status='delete')

c get target sig from somewhere
      betapsf=beta
c      open(13,file='ker.in',status='old',err=11)
c      read(13,*,err=11) betapsf
c      close(13)
c   11 continue

      call fitshistory(61,'PSF Gaussianized')
      call fitswritekey(61,'GPSFSIG',betapsf,'SIGMA OF GAUSSIANIZED PSF')
      call closefits(61)

      call closefits(71)

c     construct the real-space cov. matrix of the final noise that would
c     be obtained when adding convolved.fits + isonoise.fits
c (inverse FT of the ftacf.  Need to turn real array ftacfiso into complex)
      do l=1,m
          speq(l)=0
          do k=1,m/2
              ftacf(k,l)=abs(ftacfiso(k,l))
          enddo
      enddo
      call rlft3(ftacf,speq,m,m,1,-1)
      macf=32
      do l=1,macf
          do k=1,macf,2
              kr=k
              ki=kr+1
              kacf=mod(m+ki/2-macf/4-1,m/2)+1
              lacf=mod(m+l-macf/2-1,m)+1
              zz(kr,l)=real(ftacf(kacf,lacf))*m*m/2
              zz(ki,l)=aimag(ftacf(kacf,lacf))*m*m/2
          enddo
      enddo
      call put2('isoacf.fits',macf,macf,zz,m,m)
      write(0,*) 'Written target isotropic noise ACF to isoacf.fits'
c and write out ACF (non-isotropic) from original kernel
       do l=1,m
          speq(l)=0
          do k=1,m/2
              ftacf(k,l)=abs(filt(k,l,1))**2
          enddo
      enddo
      call rlft3(ftacf,speq,m,m,1,-1)
      macf=32
      do l=1,macf
          do k=1,macf,2
              kr=k
              ki=kr+1
              kacf=mod(m+ki/2-macf/4-1,m/2)+1
              lacf=mod(m+l-macf/2-1,m)+1
              zz(kr,l)=real(ftacf(kacf,lacf))*m*m/2
              zz(ki,l)=aimag(ftacf(kacf,lacf))*m*m/2
          enddo
      enddo
      call put2('keracf.fits',macf,macf,zz,m,m)
      write(0,*) 'Written convolved image av.noise ACF to keracf.fits'
      write(0,*) '  --ACFs in units of sigma**2 of pre-convolved image.'
      end


      subroutine shmapxpix(filt,m,nfit,cpix,cout,ix1,iy1,xm,ym,ftim,ctmp)

c convolve pixel array pix with shapelet map shmap
c Use the first npix x npix pixels of the mpix x mpix pixel array pix
c pix can be called as a real npix x xpix array - inside this routine it is
c treated as a complex npix/2 x npix array cpix
c     ftim and ctmp are temporary arrays that need to be dimensioned as
c     real npix*npix for the FT routine

      parameter(twopi=2*3.1415926535)
      complex filt(m/2,m,mmpoly2)
      complex speq(m),cpix(m/2,m),cout(m/2,m)
      complex ftim(m/2,m),ctmp(m/2,m)

c Make Fourier transform of pixel array ftim (nondestructively!)
      do iy=1,m
         do ix=1,m/2
            ftim(ix,iy)=cpix(ix,iy)
         enddo
      enddo
      call rlft3(ftim,speq,m,m,1,1)

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
c     set the nyquist freq. component to zero. (if it was not zero there is
c     a problem anyway!)
               speq(iy)=0
            enddo
c     transform back
            call rlft3(ctmp,speq,m,m,1,-1)
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
                  if (icpt.eq.1) then
                     cout(ix,iy)=yfac*complex(xfac1*z1,xfac2*z2)
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

      subroutine secplussec(z1,z2,m)
      dimension z1(m),z2(m)
c add vector z2 into vector z1
      do i=1,m
          z1(i)=z1(i)+z2(i)
      enddo
      return
      end

      subroutine shmapft(shmap,beta,n,nfit,filt,npix)
      parameter(twopi=2*3.1415926535)
      complex filt(npix/2,npix,mmpoly2),tfilt
      real xtab(n+1),ytab(n+1),shmap(mmpoly2,msh-2)

c     construct the FT of a shapelet map
c     one filter per spatial polynomial term
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
            fmax2=(sqrt(n*1.)+3)**2
            do iy=1,npix
               fynorm=(mod(iy-1+npix/2,npix)-npix/2)*twopi*beta/npix
               call htab(fynorm,n,ytab)
               do ix=1,npix/2
                  fxnorm=(mod(ix-1+npix/2,npix)-npix/2)*twopi*beta/npix
                  call htab(fxnorm,n,xtab)
                  fnorm2=fxnorm*fxnorm+fynorm*fynorm
                  tfilt=0
                  if (fnorm2.gt.fmax2) goto 7 
                  fac=exp(-0.5*fnorm2)/sqrt(3.1415926535)*beta*twopi/(npix*npix/2)
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
 7                filt(ix,iy,icpt)=tfilt
               enddo
            enddo
c            write(0,*) 'Made FT of PSF map poly cpt ',icpt,ixpow,iypow
         enddo
      enddo
      write(0,*) 'Made FT''s of PSF map poly cpts 1 to',icpt
      return
      end


      
      subroutine kerftmax(filt,m,nfit,xm,ym,nx,ny,ftacfiso)
      complex filt(m/2,m,mmpoly2)
      complex tacf(m/2,m)
      real ftacfiso(m/2,m),pacf(0:m/2),acf(m/2,m)
c 
c maximize the kernel power spectrum over position on the image, and
c then over wavenumber azimuth.  kernel is supplied as its FT filt.
c
c result is ftacfiso, the FT of the isotropised, homogenised kernel's ACF
c (ftacfiso is real)

c on an x,y grid across the image (size xm x ym) , construct the power
c spec ftacfiso(k,l) of the kernel and maximize this. Sample at 3 x nfit
c order
      do l=1,m
          do k=1,m/2
              acf(k,l)=0
          enddo
      enddo
      do ix=-3*nfit,3*nfit
          x=(0.5+0.5*ix/(3.*nfit+1))*nx
          do iy=-3*nfit,3*nfit
              y=(0.5+0.5*iy/(3.*nfit+1))*ny
              do l=1,m
                  yk=mod(l-1+m/2,m)-m/2
                  do k=1,m/2
                      xk=mod(k-1+m/2,m)-m/2
                      tacf(k,l)=0
                      icpt=0
                      do kcpt=0,nfit
                          do iypow=0,kcpt
                              ixpow=kcpt-iypow
                              icpt=icpt+1
                              if (ixpow.eq.0) then
                                  wt=1
                              else
                                  wt=(2*x/xm-1)**ixpow
                              endif
                              if (iypow.ne.0) wt=wt * (2*y/ym-1)**iypow
                              tacf(k,l)=tacf(k,l)+wt*filt(k,l,icpt)
                          enddo
                      enddo
c write out FT of ACF to as function of wavenumber
                     acf(k,l)=max(acf(k,l),abs(tacf(k,l))**2)
                  enddo
              enddo
          enddo
      enddo
c This is at each wavenumber the maximum of the kernel's power spec
c over the image. Next step is to isotropize this.
c since ACF is continuous (shapelet construction) simply tabulate the
c max ACF value on a dense grid pacf(k) of wavenumber modulus.
      do ik=0,m/2
          pacf(ik)=0
      enddo
      do k=1,m/2
          xk=mod(k-1+m/2,m)-m/2
          do l=1,m
              yk=mod(l-1+m/2,m)-m/2
              wavek=sqrt(xk*xk+yk*yk)
              ik=nint(wavek)
              if (ik.lt.m/2) then
                  pacf(ik)=max(pacf(ik),acf(k,l))
              endif
          enddo
      enddo
      open(97,file='isoacfpow1d.dat',status='unknown')
      do ik=0,m/2
          write(97,*) ik,pacf(ik)
      enddo
      close(97)
c now write the ftacfiso array with interpolated values from the amplitude pacf
      do k=1,m/2
          xk=mod(k-1+m/2,m)-m/2
          do l=1,m
              yk=mod(l-1+m/2,m)-m/2
              wavek=sqrt(xk*xk+yk*yk)
              ik=int(wavek)
              if (ik.lt.m/2-1) then
                  dw=wavek-ik
                  ftacfiso(k,l)=(1.-dw)*pacf(ik)+dw*pacf(ik+1)
              else
                  ftacfiso(k,l)=0
              endif
          enddo
      enddo
c acf now contains the desired power spectrum of the final noise
c this needs to be quadratically differenced from the ACF at each position
c and then convolved with a random noise field.
      return 
      end
