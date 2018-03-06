      real kercoef(mmpoly2,msh-2)
      dimension sky(20000)
      parameter (mpatch=210,m=256)
      dimension padded(m,m),out(m,m)
      dimension work1(m*m),work2(m*m)
      dimension zz(mpatch,mpatch),naxlen(7)
      complex filt(m/2,m,mmpoly2)
      data naxlen /7*0./ 

c as imxshmap, convolves an image by a spatially variable kernel
c expressed as a shapelet map.  differs in that here the image is split
c up into smaller sections and convolved piece by piece. 

cc----- STILL NEED TO MEASURE BACKGROUND CORRECTLY

      bg=0

      read(*,*) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(n+1)*(n+2)/2
         read(*,*) (kercoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo

      write(0,*) 'Read kernel shapelet map'
      write(0,*) 'Order for shape, spatial variation; beta:',n,nfit,beta

c prepare the filters in Fourier domain
      call shmapft(kercoef,beta,n,nfit,filt,m)

c open image to be convolved
      call openfits('inimage.fits',51,nx,ny)


cc----TO FIX: READ IN FULL PATCH, MASK OUTSIDE PART WITH BG, PROCESS, PASTE IN 
C APPROPRIATE PART 
C NEED TO ADMINISTRATE WHICH SIDE THE MARGIN IS 


c process the image in sections of mpatch x mpatch, padded to m x m
      nxsub=int((nx-1)/mpatch)+1
      nysub=int((ny-1)/mpatch)+1
      write(0,*) 'Will process image in ',nxsub,' x ',nysub,' sections.'
      write(0,*) 'Will pad subimages from ',mpatch,mpatch,' to ',m,m
c open output image
      call newfits('convolved.fits',61,nx,ny)
      naxlen(1)=nx
      naxlen(2)=ny
      do iysub=1,nysub
         iy1=1+mpatch*(iysub-1)
         iy2=iy1+mpatch-1
         if (iy2.gt.ny) then
            iy2=ny
            iy1=max(1,iy2-mpatch+1)
         endif
         do ixsub=1,nxsub
            ix1=1+mpatch*(ixsub-1)
            ix2=ix1+mpatch-1
            if (ix2.gt.nx) then
               ix2=nx
               ix1=max(1,ix2-mpatch+1)
            endif
            mx=ix2-ix1+1
            my=iy2-iy1+1
            istat=igetsec(51,nx,ny,ix1,iy1,mx,my,zz,mx*my)
            if (istat.ne.0.) then
               write(0,*) 'igetsec failed. code = ',istat
               stop
            endif
c convolve with shmap here
c pad zz array to nearest bigger power of 2 (min margin given by shapelet size)
            call pad(zz,mx,my,padded,bg,m)
            
            call shmapxpix(filt,m,nfit,padded,out,ix1,iy1,xmax,ymax,work1,work2)
c
c convolution done
            write(0,*) 'Writing subsection ',ixsub,iysub,ix1,iy1
            do j=1,my
               do i=1,mx
                  zz(i,j)=out(i,j)+bg
               enddo
            enddo
            call putsec(61,zz,mx,my,ix1,iy1,naxlen)
         enddo
      enddo

c      write(0,*) 'Convolved image in convolved.fits'
c      end


      call closefits(51)
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
         do ixpow=kcpt,0,-1
            iypow=kcpt-ixpow
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
                      x1=2*(ix+ix1-1)-1
                      z2=aimag(ctmp(ix,iy))
                      x2=2*(ix+ix1-1)
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


      subroutine pad(zz,mx,my,padded,bg,m)
      dimension zz(mx,my),padded(m,m)
c pad zz array to padded, after subtracting background level bg
      do j=1,my
         do i=1,mx
            padded(i,j)=zz(i,j)-bg
         enddo
         do i=mx+1,m
            padded(i,j)=0
         enddo
      enddo
      do j=my+1,m
         do i=1,m
            padded(i,j)=0
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


      



      SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign)
      INTEGER isign,nn1,nn2,nn3
      COMPLEX data(nn1/2,nn2,nn3),speq(nn2,nn3)
CU    USES fourn
      INTEGER i1,i2,i3,j1,j2,j3,nn(3)
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      COMPLEX c1,c2,h1,h2,w
      c1=cmplx(0.5,0.0)
      c2=cmplx(0.0,-0.5*isign)
      theta=6.28318530717959d0/dble(isign*nn1)
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      nn(1)=nn1/2
      nn(2)=nn2
      nn(3)=nn3
      if(isign.eq.1)then
        call fourn(data,nn,3,isign)
        do 12 i3=1,nn3
          do 11 i2=1,nn2
            speq(i2,i3)=data(1,i2,i3)
11        continue
12      continue
      endif
      do 15 i3=1,nn3
        j3=1
        if (i3.ne.1) j3=nn3-i3+2
        wr=1.0d0
        wi=0.0d0
        do 14 i1=1,nn1/4+1
          j1=nn1/2-i1+2
          do 13 i2=1,nn2
            j2=1
            if (i2.ne.1) j2=nn2-i2+2
            if(i1.eq.1)then
              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
              data(1,i2,i3)=h1+h2
              speq(j2,j3)=conjg(h1-h2)
            else
              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
              data(i1,i2,i3)=h1+w*h2
              data(j1,j2,j3)=conjg(h1-w*h2)
            endif
13        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          w=cmplx(sngl(wr),sngl(wi))
14      continue
15    continue
      if(isign.eq.-1)then
        call fourn(data,nn,3,isign)
      endif
      return
      END
