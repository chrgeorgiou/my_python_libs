      real sh(msh),psfcoef(mmpoly2,msh-2)
c read on a psf map, and at requested x and y write out pixel values.

      read(*,*) x,y
 1    read(*,*,err=1) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(n+1)*(n+2)/2
         read(*,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo
      nplot=max(10,nint(6*beta))
      call interpolpsf(x,y,xmax,ymax,psfcoef,n,beta,nfit,sh)
      do i=-nplot,nplot
          do j=-nplot,nplot
              if (i*i+j*j.lt.nplot*nplot) then
                  write(*,*) i,j,shapeletpix(real(i),real(j),sh)
              endif
          enddo
      enddo
      write(0,*) 'Evaluated PSF map at x,y=',x,y
      end

