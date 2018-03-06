c read in psf map and a source catalogue; at each source position
c calculate PSF model and KSB-like polarization

      parameter (pi=3.1415926535)
      real sh(msh),psfcoef(mmpoly2,msh-2)

      open(10,file='psf.map',status='old')
      read(10,*) npsf,betapsf
      read(10,*) nfit,xmax,ymax
      read(10,*)
      do k=1,(npsf+1)*(npsf+2)/2
         read(10,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(10,*)
      enddo
      close(10)
      write(0,*) 'Read psf model from psf.map'
      write(0,*) 'Order for shape, spatial variation; beta:',npsf,nfit,betapsf

      do i=1,mmcat
    1     read(*,*,err=1,end=2) id,x,y
          call interpolpsf(x,y,xmax,ymax,psfcoef,npsf,betapsf,nfit,sh)          
          call centershape(sh,sh,dx,dy)
          a00=sh(3)
c          a10=sh(4)
c          a01=sh(5)
          a20=sh(6)
          a11=sh(7)
          a02=sh(8)
          txx=betapsf**3*sqrt(pi)*(0.5*a00+sqrt(0.5)*a20)
          txy=betapsf**3*sqrt(pi)*0.5*a11
          tyy=betapsf**3*sqrt(pi)*(0.5*a00+sqrt(0.5)*a02)
          e1=(txx-tyy)/(txx+tyy)
          e2=2*txy/(txx+tyy)
          write(*,*) id,e1,e2
      enddo
    2 end


