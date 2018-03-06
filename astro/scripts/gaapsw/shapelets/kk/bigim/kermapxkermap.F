      real kercoef(mmpoly2,msh-2),sh1(msh),sh2(msh),sh3(msh)
      real kercoef2(1,msh-2)

c routine to convolve a spatially variable shapelet map with a constant shapelet
c 
c input: variable map concatinated with constant shapelet map
c
c the output map has the same scale radius as the first input map
c
c Can be used to correct a noise autocorrelation function after
c convolution by a shapelet map with a shapelet model of the original
c pixel correlation function.
c

c read variable shapelet map
      read(*,*) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(n+1)*(n+2)/2
         read(*,*) (kercoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo

      write(0,*) 'Read variable kernel shapelet map'
      write(0,*) 'Order for shape, spatial variation; beta:',n,nfit,beta

c read constant shapelet map
      read(*,*) n2,beta2
      read(*,*) nfit2,xmax2,ymax2
      if (nfit2.ne.0) then
          write(0,*) 'Second shapelet map varies spatially.'
          write(0,*) 'Cannot deal with this. Order =',nfit2
          write(0,*) 'Only the average will be used.'
      endif
      read(*,*)
      do k=1,(n2+1)*(n2+2)/2
         read(*,*) kercoef2(1,k),(t,kfit=2,(nfit2+1)*(nfit2+2)/2)
         read(*,*)
      enddo

      write(0,*) 'Read constant kernel shapelet map to convolve with.'
      write(0,*) 'Order for shape, beta:',n2,beta2

c set scale radius for output map
      betaout=beta

      sh2(1)=beta2
      sh2(2)=n2
      do k=1,(n2+1)*(n2+2)/2
          sh2(k+2)=kercoef2(1,k)
      enddo
      
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
