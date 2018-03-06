      real kercoef(mmpoly2,msh-2),sh(msh),shacf(msh)
      dimension naxlen(7)
      data naxlen /7*0./ 
      dimension xgrid((4*mmpoly+1)*(4*mmpoly+1)),
     :    ygrid((4*mmpoly+1)*(4*mmpoly+1))

c construct the ACF of a kernel map; write out as a grid of ACFs that
c can be fitted with fitkermap.  Describe ACF as shapelets with scale
c radius sqrt(2) * original; same shapelet order.

      read(*,*) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(n+1)*(n+2)/2
         read(*,*) (kercoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo
      write(0,*) 'Read kernel shapelet map'
      write(0,*) 'Order for shape, spatial variation; beta:',n,nfit,beta
      nsh=(n+1)*(n+2)/2+2

c now make a grid of kernels, sampling frequency matched to poly order 
c isubsample specifies how dense the grid is.
      kgrid=0
      isubsample=3
c set output scale radius as sqrt(2)*beta
      betaacf=sqrt(2.)*beta
      do ix=1,isubsample*nfit+1
          x=xmax*ix/(isubsample*nfit+2.)
          do iy=1,isubsample*nfit+1
              y=ymax*iy/(isubsample*nfit+2.)
              kgrid=kgrid+1
              xgrid(kgrid)=x
              ygrid(kgrid)=y
              call interpolpsf(x,y,xmax,ymax,kercoef,n,beta,nfit,sh)
              call correlate(sh,sh,betaacf,shacf)
              sum=seriessum(shacf)
c write out these ACF coefs in psf format
              write(*,1001) x,y,1.,1e-3,betaacf,1.,1.,0.,0,0
 1001         format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i10,i6)
              write(*,*) 1.,0.,0.,0.,10.,1.e-3
              write(*,*) betaacf,n
              write(*,*) (shacf(k)/sum,k=3,nsh)
              write(*,*)
          enddo
      enddo
      end
