c     set up a gaussian source, stddev sqrt(a^2+p^2), and a PSF with
c     stddev p expand in shapelets with scale bs and bp
c
c     construct PSF matrix P that convolves a source with scale q with a
c     PSF of scale bp to yield a result with scale bs
c
c     invert P, apply to source to give deconvolved source = gaussian of
c     scale a, expressed with scale bq
c
c     construct new PSF = gaussian of scale and dispersion q
c
c     construct new PSF matrix that convolves a source of scale q and a
c     PSF of scale q to a result of scale q
c
c     read off the 00 component of the result, which is the
c     gaussian-weighted aperture photometry of the source convolved with
c     gaussian PSF of scale q.

      dimension z(101,101),psf(101,101)
      dimension shobs(msh),shpsf(msh),shpsf2(msh)
      dimension qcoef(0:mm)
      dimension pmat(msh-2,msh-2),pinv(msh-2,msh-2),indx(msh-2),photvec(msh-2)
      parameter (pi=3.1415926535)

c a,p,q are gaussian dispersion of intrinsic gal, PSF, target PSF
      a=4
      p=3
      betaobs=6
      betapsf=4
      q=5
      n=8
      nsh=(n+1)*(n+2)/2
      do i=1,101
         x=i-51
         do j=1,101
            y=j-51
            r2=x*x+y*y
            z(i,j)=exp(-r2/2/(a*a+p*p))/(2*pi*(a*a+p*p))
            psf(i,j)=exp(-r2/2/(p*p))/(2*pi*p*p)
         enddo
      enddo
      sg=0
      sp=0
      do i=1,101
         do j=1,101
            sg=sg+z(i,j)
            sp=sp+psf(i,j)
         enddo
      enddo
      write(0,*) 'Integral, sig of galaxy and PSF:',sg,a,sp,p

      call fitshapelets(z,101,101,101,51.,51.,4*betaobs,0.,shobs,n,betaobs)
      call writeshape('OBS GAL',shobs)
      call fitshapelets(psf,101,101,101,51.,51.,4*betapsf,0.,shpsf,n,betapsf)
      call writeshape('PSF',shpsf)

      do kkk=1,1


c setpsfmat: PSF shape, result scale, input scale, matrix
      call setpsfmat(shpsf,betaobs,q,pmat)

c      do i=1,nsh
c         do j=1,nsh
c            write(0,*) 'P',i,j,pmat(i,j)
c         enddo
c      enddo

      do i=1,nsh
         do j=1,nsh
            pinv(i,j)=0
         enddo
         pinv(i,i)=1
      enddo
      call ludcmp(pmat,nsh,msh-2,indx,detfac)

      k=0
      do nn=0,n
          do j=0,nn
              i=nn-j
              k=k+1
c only calculate the even/even modes
              if (mod(i,2)+mod(j,2).eq.0) then
                  call lubksb(pmat,nsh,msh-2,indx,pinv(1,k))
              endif
          enddo
      enddo
      
c      do i=1,nsh
c         do j=1,nsh
c            write(0,*) 'P^-1',i,j,pinv(i,j)
c         enddo
c      enddo
      shpsf2(1)=q
      shpsf2(2)=n
      shpsf2(3)=1/(2*sqrt(pi)*q)
      do k=4,msh
         shpsf2(k)=0
      enddo
c      call writeshape('NEW PSF',shpsf2)

C
C THE FOLLOWING PSF MATRIX DOES NOT NEED TO BE CALCULATED IN ALL DETAIL
C     IN ANY CASE IT WILL BE VERY SPARSE BETTER WORK OUT DIRECTLY WHAT
C     THE INTEGRAL IS OF A BASIS FUNC TIMES A GAUSSIAN?
c 

c      call setpsfmat(shpsf2,q,q,pmat)
c      do i=1,nsh
c         do j=1,nsh
c            write(0,*) 'PMAT2',i,j,pmat(1,j)
c         enddo
c      enddo

c     coefficients that represent the gaussian-weighted flux in a
c     shapelet of scale q, order k, that is vonvolved with a gaussian of
c     dispersion q (all 1-D):

      qcoef(0)=sqrt(2*q*sqrt(pi)/3)
      do k=2,n,2
          qcoef(k-1)=0
          qcoef(k)=qcoef(k-2) /(3*k) * sqrt(k*(k-1.))
      enddo

c      k=0
c      do nn=0,n
c          do j=0,nn
c              i=nn-j
c              k=k+1
c              if (mod(i,2)+mod(j,2).eq.0) then
c                  write(0,*) i,j,'PMAT2:',pmat(1,k),qcoef(i)*qcoef(j)
c              endif
c          enddo
c      enddo



c      k=0
c      do nn=0,n
c         do j=0,nn
c             i=nn-j
c             k=k+1
c             photvec(k)=0
cc only calculate the even/even modes
c             if (mod(i,2)+mod(j,2).eq.0) then
c                 do l=1,nsh
c                     photvec(k)=photvec(k)+qcoef(i)*qcoef(j)*pinv(l,k)
c                 enddo
c             endif
c         enddo
c      enddo

      phot00=0

      k=0
      do nn=0,n
          do j=0,nn
              i=nn-j
              k=k+1
              if (mod(i,2)+mod(j,2).eq.0) then
                  do l=1,nsh
                      phot00=phot00+qcoef(i)*qcoef(j)*pinv(k,l)*shobs(l+2)
                  enddo
              endif
          enddo
      enddo

c      do k=1,nsh
c         phot00=phot00+photvec(k)*shobs(k+2)
c         write(0,'(i4,2f12.5)') k,photvec(k),shobs(k+2)
c      enddo


      enddo
      write(0,*) 'Exper, theor phot00 values:',kkk,phot00,q*q/(a*a+2*q*q)

      end
