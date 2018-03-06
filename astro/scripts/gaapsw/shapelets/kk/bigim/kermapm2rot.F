c read in a kernel map, and write out the anisotropic part, rotated by
c 90 degrees
c to be used for isotropizing the noise.

      real sh(msh),coef(mmpoly2,msh-2),sh0(msh),cc(0:mm/2)

c get sigma of target gaussian

c read PSF map
      read(*,*) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(n+1)*(n+2)/2
         read(*,*) (coef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo
      write(0,*) 'Read kernel model:'
      write(0,*) 'Order for shape, spatial variation; beta:',n,nfit,beta

      do ish=1,(n+1)*(n+2)/2
         sh(ish+2)=coef(1,ish)
      enddo
      sh(1)=beta
      sh(2)=n
      write(0,*) 'Average PSF component integral is ',seriessum(sh)

c construct circular polar shapelet coefs at order nn=0,2,...,n
      do ish=3,(n+1)*(n+2)/2+2
         sh0(ish)=0
      enddo

      do nn=0,n,2
          k0=3+nn*(nn+1)/2
          m=nn/2
          cc(0)=1
          sum2=cc(0)**2
          do i=1,m
              i2=2*i
              mi2=2*(m-i)
              cc(i)=cc(i-1)*(m-i+1)/i*sqrt(i2*(i2-1)/(mi2+2.)/(mi2+1.))
              sum2=sum2+cc(i)**2
c     write(0,*) 'rnd',nn,i2,cc
          enddo
          sum2=sqrt(sum2)
          do i=0,m
              sh0(k0+2*i)=cc(i)/sum2
          enddo
          write(0,*) nn,(sh0(k0+i),i=0,nn,2)
      enddo

c now go spatial order by spatial order
      do kfit=1,(nfit+1)*(nfit+2)/2
          do ish=1,(n+1)*(n+2)/2
              sh(ish+2)=coef(kfit,ish)
          enddo
c shapelet order by shapelet order
c determine the round component amplitude (even orders only)
          do nn=0,n,2
              sum=0
              k0=3+nn*(nn+1)/2
              do i=0,nn,2
                  sum=sum+sh0(k0+i)*sh(k0+i)
              enddo
c and subtract it
              do i=0,nn,2
                  sh(k0+i)=sh(k0+i)-sum*sh0(k0+i)
              enddo
          enddo
c shapelet by shapelet order, rotate by 90 degrees (all orders)
c this means x-->y, y-->-x
c so swap shapelet cpts ab and ba, and change sign if power of y was odd
          do nn=0,n
              k0=3+nn*(nn+1)/2
c swap ij and ji
              do i=0,nn/2
                  j=nn-i
                  ttt=sh(k0+j)
                  sh(k0+j)=sh(k0+i)
                  sh(k0+i)=ttt
              enddo
c change sign of all odd-y powers
              do i=nn-mod(nn+1,2),0,-2
                  sh(k0+i)=-sh(k0+i)
              enddo
          enddo
c check integral is 0
          write(0,*) 'Cpt ',kfit,': integral is ',seriessum(sh)
c     write back into the kernel map
          do ish=1,(n+1)*(n+2)/2
              coef(kfit,ish)=sh(ish+2)
          enddo
      enddo
      
      write(*,*) n,beta
      write(*,*) nfit,xmax,ymax
      write(*,*)
      do ish=1,(n+1)*(n+2)/2
         write(*,*) (coef(kfit,ish),kfit=1,(nfit+1)*(nfit+2)/2)
         write(*,*)
      enddo
      end


 


