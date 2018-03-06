      real shape1(msh),shape2(msh),z1(30,30),z2(30,30)

      beta=1.8
      nmax=3
      shape1(1)=beta
      shape2(1)=beta
      shape1(2)=nmax
      shape2(2)=nmax
      do i=3,2+(nmax+1)*(nmax+2)/2
          do j=3,2+(nmax+1)*(nmax+2)/2
              do k=3,2+(nmax+1)*(nmax+2)/2
                  shape1(k)=0
                  shape2(k)=0
              enddo
              shape1(i)=1
              shape2(j)=1
              call decode(shape1,z1,30,25)
              call decode(shape2,z2,30,25)
              sum11=0
              sum12=0
              sum22=0
              do ix=1,25
                  do iy=1,25
                      sum11=sum11+z1(ix,iy)**2
                      sum12=sum12+z1(ix,iy)*z2(ix,iy)
                      sum22=sum22+z2(ix,iy)**2
                  enddo
              enddo
              write(0,'(2i4,3f12.5)') i-2,j-2,sum11,sum12,sum22
          enddo
      enddo
      end
