c read in a psf map (polynomial interpolation coefs for PSF shapelet
c coefs) and construct a spatially variable convolution kernel that
c makes the PSF Gaussian with dispersion beta.

      dimension psfcoef(mmpoly2,msh-2),pmat(msh-2,msh-2),indx(msh-2),
     :    sh(msh)
      dimension xgrid((2*mmpoly+1)*(2*mmpoly+1)),ygrid((2*mmpoly+1)*(2*mmpoly+1))

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

c set output Gaussian radius
      open(13,file='ker.in',status='old',err=1)
      goto 2
    1 write(0,*) '****Specify target PSF size or multiplier in file ker.in****'
      stop
    2 read(13,*) fac
      close(13)

      if (fac.lt.0.) then
          betaout=betapsf*abs(fac)
          write(0,*) 'Output Gaussian radius is PSF beta x ',abs(fac)
          write(0,*) 'Output Gaussian radius will be ',betaout
          open(13,file='ker.in',status='unknown')
          write(13,*) betaout
          close(13)
      else
          betaout=fac
      endif
      betaker=sqrt(betapsf*betaout)

c now make a grid of PSF matrices, sampling frequency matched to poly order 
c isubsample specifies how dense the grid is.
      kgrid=0
      isubsample=4
      do ix=1,isubsample*nfit+1
          x=xmax*ix/(isubsample*nfit+2.)
          do iy=1,isubsample*nfit+1
              y=ymax*iy/(isubsample*nfit+2.)
              kgrid=kgrid+1
              xgrid(kgrid)=x
              ygrid(kgrid)=y
              call interpolpsf(x,y,xmax,ymax,psfcoef,n,betapsf,nfit,sh)
              call setpsfmat(sh,betaout,betaker,pmat)
              call ludcmp(pmat,nsh,msh-2,indx,detfac)
c set up target Gaussian PSF, normalize its flux to 1
              sh(1)=betaout
              sh(2)=n
              sh(3)=1
              do k=4,nsh+2
                  sh(k)=0
              enddo
              s=seriessum(sh)
              sh(3)=sh(3)/s
              call lubksb(pmat,nsh,msh-2,indx,sh(3))
              sh(1)=betaker
              sum=seriessum(sh)
c write out these kernel coefs in psf format
              write(*,'(''# target Gaussian PSF sigma: '',g15.5)') betaout
              write(*,'(''#'')')
              write(*,1001) x,y,1.,1e-3,betapsf,1.,1.,0.,0,0
 1001         format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i10,i6)
              write(*,*) sum,0.,0.,0.,10.,1.e-3
              write(*,*) betaker,n
              write(*,*) (sh(k)/sum,k=3,(n+1)*(n+2)/2+2)
              write(*,*)
          enddo
      enddo
      end

          
