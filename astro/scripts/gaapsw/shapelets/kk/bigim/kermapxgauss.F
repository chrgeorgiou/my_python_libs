      real kercoef(mmpoly2,msh-2),sh1(msh),sh2(msh),sh3(msh)

c routine to convolve a shapelet map with a gaussian
c input: sigma to convolve with, output scale radius
      read(*,*) betaconv,betaout

      write(0,*) 'Will convolve kernel map with gauss, sig=',betaconv
      write(0,*) 'Result expressed as shapelet map with beta=',betaout
c read shapelet map
      read(*,*) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(n+1)*(n+2)/2
         read(*,*) (kercoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo

      write(0,*) 'Read kernel shapelet map'
      write(0,*) 'Order for shape, spatial variation; beta:',n,nfit,beta

      sh2(1)=betaconv
      sh2(2)=0
      sh2(3)=0.5/sqrt(3.1415926535)/betaconv
      
      do kfit=1,(nfit+1)*(nfit+2)/2
          sh1(1)=beta
          sh1(2)=n
          do k=1,(n+1)*(n+2)/2
              sh1(k+2)=kercoef(kfit,k)
          enddo
          call convolve(sh1,sh2,betaout,sh3)
          do k=1,(n+1)*(n+2)/2
              kercoef(kfit,k)=sh3(k+2)
          enddo
      enddo
      write(*,*) n,betaout
      write(*,*) nfit,xmax,ymax
      write(*,*) 
      do k=1,(n+1)*(n+2)/2
         write(*,*) (kercoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         write(*,*)
      enddo
      end
