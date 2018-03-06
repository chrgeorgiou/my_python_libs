      real kercoef(mmpoly2,msh-2)
      dimension sky(20000)
      parameter (m=256)
      dimension zz(m,m),out(m,m)
      dimension work1(m*m),work2(m*m)
      dimension naxlen(7)
      complex filt(m/2,m,mmpoly2)
      data naxlen /7*0./ 

c as imxshmap, convolves an image by a spatially variable kernel
c expressed as a shapelet map.  differs in that here the image is split
c up into smaller sections and convolved piece by piece. 


c perhaps best to subtract a bg and reinsert it in the end - not done so far!

      read(*,*) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(n+1)*(n+2)/2
         read(*,*) (kercoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo

      write(0,*) 'Read kernel shapelet map'
      write(0,*) 'Order for shape, spatial variation; beta:',n,nfit,beta

c calculate width of anti-alias border
      nmask=nint(beta*(sqrt(1.*n))+5.)
      if (nmask.gt.m/4) then
         write(0,*) 'Shapelet map too fat for Fourier convolution.'
         stop
      endif

c prepare the filters in Fourier domain
      call shmapft(kercoef,beta,n,nfit,filt,m)

c open image to be convolved
      call openfits('inimage.fits',51,nx,ny)


      mpatch=m-2*nmask
c process the image in sections of mpatch x mpatch, padded to m x m
c outer border of nmask pix wide cannot be used.
      nxsub=int((nx-1-2*nmask)/mpatch)+1
      nysub=int((ny-1-2*nmask)/mpatch)+1
      write(0,*) 'Will process image in ',nxsub,' x ',nysub,' sections.'
      write(0,*) 'Will pad subimages from ',mpatch,mpatch,' to ',m,m
c open output image
c      call newfits('convolved.fits',61,nx,ny)
      naxlen(1)=nx
      naxlen(2)=ny

      call copyfits(51,61,'convolved.fits',5)

c go through the input image, read in a patch, mask a border around it,
c convolve, write out
      do iysub=1,nysub
         if (iysub.eq.1) then
            iy1=1
         else
            iy1=iy1+mpatch
         endif
         iy2=iy1+m-1
         if (iy2.gt.ny) then
            iy2=ny
            iy1=max(1,iy2-m+1)
         endif
         ix1=1
         do ixsub=1,nxsub
            if (ixsub.eq.1) then
               ix1=1
            else
               ix1=ix1+mpatch
            endif
            ix2=ix1+m-1
            if (ix2.gt.nx) then
               ix2=nx
               ix1=max(1,ix2-m+1)
            endif
            mx=ix2-ix1+1
            my=iy2-iy1+1
            istat=igetsec(51,nx,ny,ix1,iy1,mx,my,zz,mx*my)
            if (istat.ne.0.) then
               write(0,*) 'igetsec failed. code = ',istat
               stop
            endif
c convolve with shmap here
            
            call shmapxpix(filt,m,nfit,zz,out,ix1,iy1,xmax,ymax,work1,work2)
c
c convolution done, result in out
c            write(0,*) 'Writing subsection ',ixsub,iysub,ix1,iy1
c pick out the part of out that is unaffected by aliasing (strip off border)
            call stripborder(out,zz,mx-2*nmask,my-2*nmask,m,nmask)
            call putsec(61,zz,mx-2*nmask,my-2*nmask,ix1+nmask,iy1+nmask,naxlen)
         enddo
         write(0,*) 'Finished strip ',iysub,' / ',nysub,
     :       '= pix ',iy1+nmask,' - ',iy1+my-nmask-1
      enddo

c      write(0,*) 'Convolved image in convolved.fits'
c      end


      call closefits(51)

c get target sig from somewhere
      betapsf=0
      open(13,file='ker.in',status='old',err=11)
      read(13,*,err=11) betapsf
      close(13)
   11 continue

      call fitshistory(61,'PSF Gaussianized')
      call fitswritekey(61,'GPSFSIG',betapsf,'SIGMA OF GAUSSIANIZED PSF')

      call closefits(61)
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
c add result into the output pixel array with proper weight
c write out full array including padded parts - could trim this
c              write(0,*) 'Adding in part for x,y powers ',ixpow,iypow
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
c complete the loops over polynomial coefs of the convolution kernel
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


      


c FFT routine already included in workhorse.F file

c      SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign)
c      INTEGER isign,nn1,nn2,nn3
c      COMPLEX data(nn1/2,nn2,nn3),speq(nn2,nn3)
cCU    USES fourn
c      INTEGER i1,i2,i3,j1,j2,j3,nn(3)
c      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
c      COMPLEX c1,c2,h1,h2,w
c      c1=cmplx(0.5,0.0)
c      c2=cmplx(0.0,-0.5*isign)
c      theta=6.28318530717959d0/dble(isign*nn1)
c      wpr=-2.0d0*sin(0.5d0*theta)**2
c      wpi=sin(theta)
c      nn(1)=nn1/2
c      nn(2)=nn2
c      nn(3)=nn3
c      if(isign.eq.1)then
c        call fourn(data,nn,3,isign)
c        do 12 i3=1,nn3
c          do 11 i2=1,nn2
c            speq(i2,i3)=data(1,i2,i3)
c11        continue
c12      continue
c      endif
c      do 15 i3=1,nn3
c        j3=1
c        if (i3.ne.1) j3=nn3-i3+2
c        wr=1.0d0
c        wi=0.0d0
c        do 14 i1=1,nn1/4+1
c          j1=nn1/2-i1+2
c          do 13 i2=1,nn2
c            j2=1
c            if (i2.ne.1) j2=nn2-i2+2
c            if(i1.eq.1)then
c              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
c              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
c              data(1,i2,i3)=h1+h2
c              speq(j2,j3)=conjg(h1-h2)
c            else
c              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
c              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
c              data(i1,i2,i3)=h1+w*h2
c              data(j1,j2,j3)=conjg(h1-w*h2)
c            endif
c13        continue
c          wtemp=wr
c          wr=wr*wpr-wi*wpi+wr
c          wi=wi*wpr+wtemp*wpi+wi
c          w=cmplx(sngl(wr),sngl(wi))
c14      continue
c15    continue
c      if(isign.eq.-1)then
c        call fourn(data,nn,3,isign)
c      endif
c      return
c      END
