c     verify that the interpolated PSF map has unit integral all over
c     the field: i.e. sums to 1 in the 00 component, and all spatial
c     variation cpts integrate to zero.

      parameter (pi=3.1415926535)
      dimension psfcoef(mmpoly2,msh-2),sh(msh)

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

      do ifit=1,(nfit+1)*(nfit+2)/2
          sh(1)=betapsf
          sh(2)=n
          do k=1,nsh
              sh(k+2)=psfcoef(ifit,k)
          enddo
          write(*,*) ifit,seriessum(sh)
      enddo
      end
