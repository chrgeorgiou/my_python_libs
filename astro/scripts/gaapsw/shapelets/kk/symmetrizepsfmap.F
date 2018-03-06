      
      parameter (mm012=5*(mm/2)-1)
      real psfcoef(mmpoly2,msh-2),roundpsfcoef(mmpoly2,msh-2)
      
c     symmetrize a PSF map: remove all odd coefficients, and the m=2 and
c     m=4 parts.  
c     i.e., a set c_ab --> C_(ab) when a,b even; zero otherwise.
c
c this is equivalent to equalizing x and -x, y and -y directions, as
c well as x and y.

c read input PSF map from std input

      read(*,*) npsf,betapsf
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(npsf+1)*(npsf+2)/2
         read(*,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo
      write(0,*) 'Read psf model from psf.map'
      write(0,*) 'Order for shape, spatial variation; beta:',npsf,nfit,betapsf

c number of poly coefs per shapelet term; number of shapelet terms
      npoly=(nfit+1)*(nfit+2)/2
      nsh=(npsf+1)*(npsf+2)/2

      do kfit=1,npoly
         k=0
         do kk=1,nsh
            roundpsfcoef(kfit,kk)=0
         enddo
         do n=0,npsf
            do m=0,n
               l=n-m
               k=k+1
c only create terms with even total order
               if (mod(n,2).eq.0) then
c find index of term with x,y swapped
                  kswap=k+(l-m)
                  if (mod(m,2).eq.0) then
                     roundpsfcoef(kfit,k)=0.5*(psfcoef(kfit,k)+psfcoef(kfit,kswap))
                  else
                     roundpsfcoef(kfit,k)=0.5*(psfcoef(kfit,k)-psfcoef(kfit,kswap))
                  endif
               endif
            enddo
         enddo
      enddo
      write(*,*) npsf,betapsf
      write(*,*)  nfit,xmax,ymax
      write(*,*)
      do k=1,(npsf+1)*(npsf+2)/2
         write(*,*) (roundpsfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         write(*,*)
      enddo
      write(0,*) 'Wrote out psf map with all odd m, and m=2 filtered out.'
      end
              
