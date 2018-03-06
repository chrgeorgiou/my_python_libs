c read in a psf map (polynomial interpolation coefs for PSF shapelet
c coefs) and construct a spatially variable convolution kernel that
c makes the PSF Gaussian with dispersion beta.

      parameter (pi=3.1415926535)
      dimension psfcoef(mmpoly2,msh-2),pmat(msh-2,msh-2),indx(msh-2),
     :    sh(msh),emat(msh-2,msh-2)
      dimension xgrid((4*mmpoly+1)*(4*mmpoly+1)),
     :    ygrid((4*mmpoly+1)*(4*mmpoly+1))

c first get the PSF map coefficients from stdin

      read(*,*) n,betapsf
      read(*,*) nfit,xmax,ymax
      read(*,*)
      nsh=(n+1)*(n+2)/2
      do k=1,nsh
         read(*,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo

      write(0,*) 'Read psf model'
      write(0,*) 'Order for shape, spatial variation; beta:',n,nfit,betapsf

c calculate ellipticity of the average PSF, using KSB-like 2nd moments
      a00=psfcoef(1,1)
      a10=psfcoef(1,2)
      a01=psfcoef(1,3)
      a20=psfcoef(1,4)
      a11=psfcoef(1,5)
      a02=psfcoef(1,6)
      txx=betapsf**3*sqrt(pi)*(0.5*a00+sqrt(0.5)*a20)
      txy=betapsf**3*sqrt(pi)*0.5*a11
      tyy=betapsf**3*sqrt(pi)*(0.5*a00+sqrt(0.5)*a02)
      e1=(txx-tyy)/(txx+tyy)
      e2=2*txy/(txx+tyy)
      ell=sqrt(e1*e1+e2*e2)
      write(0,*) 'Ellipticity of the average PSF: ',ell

c set output Gaussian radius
      open(13,file='ker.in',status='old',err=1)
      goto 2
    1 write(0,*) '****Specify target PSF size or multiplier in ker.in !'
      stop
    2 read(13,*) fac
      close(13)

      if (fac.lt.0.) then
          betaker=betapsf*abs(fac)/(1-ell)
          write(0,*) 'Output Gaussian radius is PSF beta x ',abs(fac)
          write(0,*) 'Output Gaussian radius will be ',betaker
          open(13,file='ker.in',status='unknown')
          write(13,*) betaker
          close(13)
      else
          betaker=fac
          write(0,*) 'Requested output Gaussian radius is ',betaker
      endif

c check the grid of PSFs to check whether PSF peak is higher than target PSF by factor 1.1 at least. If not increase beta further.

      isubsample=2
      pkmin=1e22
      do ix=1,isubsample*nfit+1
         x=xmax*ix/(isubsample*nfit+2.)
         do iy=1,isubsample*nfit+1
            y=ymax*iy/(isubsample*nfit+2.)
            call interpolpsf(x,y,xmax,ymax,psfcoef,n,betapsf,nfit,sh)
            pk=shapeletpix(0.,0.,sh)
            pkmin=min(pk,pkmin)
         enddo
      enddo
      pkfac=0.90
      if (1./(2*3.1416*betaker**2) .gt. pkfac*pkmin) then
         betaker=sqrt(1./(2*3.1416*pkfac*pkmin))
         write(0,*) '**Peak of target PSF higher than ',nint(100*pkfac),'% of original!'
         write(0,*) 'Lowest peak PSF = ',pkmin,'   target peak ',1./(2*3.1416*betaker**2)
         write(0,*) 'Output Gaussian radius further increased to ',betaker
         open(13,file='ker.in',status='unknown')
         write(13,*) betaker
         close(13)
      endif

c now make a grid of PSF matrices, sampling frequency matched to poly order 
c isubsample specifies how dense the grid is.
      kgrid=0
      isubsample=4
c set up shapelet order for output as n-2
      nshout=n*(n-1)/2
      do ix=1,isubsample*nfit+1
          x=xmax*ix/(isubsample*nfit+2.)
          do iy=1,isubsample*nfit+1
              y=ymax*iy/(isubsample*nfit+2.)
              kgrid=kgrid+1
              xgrid(kgrid)=x
              ygrid(kgrid)=y
              call interpolpsf(x,y,xmax,ymax,psfcoef,n,betapsf,nfit,sh)
              call setpsfmat(sh,betaker,betaker,pmat)
c set up target Gaussian PSF, normalize its flux to 1
              sh(1)=betaker
              sh(2)=n-2
              sh(3)=1
              do k=4,nshout+2
                  sh(k)=0
              enddo
              s=seriessum(sh)
c construct Mij = Pki Pkj matrix, for i,j up to shapelet terms of order n-2 
              do j=1,nshout
                  do i=1,nshout
                      emat(i,j)=0
                      do k=1,nsh
                          emat(i,j)=emat(i,j)+pmat(k,i)*pmat(k,j)
                      enddo
                  enddo
              enddo
c construct Gi = Pik gk where gk is the target shapelet = (1/s,0,0,0,0,0...)
              do i=1,nshout
                  sh(i+2)=pmat(1,i)/s
              enddo
              call ludcmp(emat,nshout,msh-2,indx,detfac)
              call lubksb(emat,nshout,msh-2,indx,sh(3))
              sh(1)=betaker
              sum=seriessum(sh)
c write out these kernel coefs in psf format
              write(*,'(''# target Gaussian PSF sigma: '',g15.5)') betaker
              write(*,'(''#'')')
              write(*,1001) x,y,1.,1e-3,betapsf,1.,1.,0.,0,0
 1001         format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i10,i6)
              write(*,*) 1.,0.,0.,0.,10.,1.e-3
              write(*,*) betaker,n-2
              write(*,*) (sh(k)/sum,k=3,nshout+2)
              write(*,*)
          enddo
      enddo
      end

          
