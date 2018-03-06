
c operators that operate on single shape:
c translation, rotation, shear, convergence
c note that translations are in units of beta, 
c                  i.e. 0.1*TRANS1 shifts by 0.1*beta

      subroutine trans1(inshape,outshape)
c generator for translate in x
      real inshape(msh),outshape(msh),tmpshape1(msh),tmpshape2(msh)
      call a1(inshape,tmpshape1)
      call a1d(inshape,tmpshape2)
      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=sqrt(0.5)*(tmpshape2(k)-tmpshape1(k))
      enddo
      return
      end


      subroutine trans2(inshape,outshape)
c generator for translate in y
      real inshape(msh),outshape(msh),tmpshape1(msh),tmpshape2(msh)
      call a2(inshape,tmpshape1)
      call a2d(inshape,tmpshape2)
      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=sqrt(0.5)*(tmpshape2(k)-tmpshape1(k))
      enddo
      return
      end


      subroutine conv(inshape,outshape)
c generator for convergence
      real inshape(msh),outshape(msh),tmpshape1(msh),tmpshape2(msh)

      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=inshape(k)
      enddo
      call a1d(inshape,tmpshape1)
      call a1d(tmpshape1,tmpshape2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=outshape(k)+0.5*tmpshape2(k)
      enddo
      call a2d(inshape,tmpshape1)
      call a2d(tmpshape1,tmpshape2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=outshape(k)+0.5*tmpshape2(k)
      enddo
      call a1(inshape,tmpshape1)
      call a1(tmpshape1,tmpshape2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=outshape(k)-0.5*tmpshape2(k)
      enddo
      call a2(inshape,tmpshape1)
      call a2(tmpshape1,tmpshape2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=outshape(k)-0.5*tmpshape2(k)
      enddo
      return
      end


      subroutine rot(inshape,outshape)
c generator for rotation
      real inshape(msh),outshape(msh),tmpshape1(msh),tmpshape2(msh)

      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      call a2d(inshape,tmpshape1)
      call a1(tmpshape1,tmpshape2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=tmpshape2(k)
      enddo
      call a2(inshape,tmpshape1)
      call a1d(tmpshape1,tmpshape2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=outshape(k)-tmpshape2(k)
      enddo
      return
      end


      subroutine shear1(inshape,outshape)
c generator for shear1
      real inshape(msh),outshape(msh),tmpshape1(msh),tmpshape2(msh)

      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      call a1d(inshape,tmpshape1)
      call a1d(tmpshape1,tmpshape2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=0.5*tmpshape2(k)
      enddo
      call a2d(inshape,tmpshape1)
      call a2d(tmpshape1,tmpshape2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=outshape(k)-0.5*tmpshape2(k)
      enddo
      call a1(inshape,tmpshape1)
      call a1(tmpshape1,tmpshape2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=outshape(k)-0.5*tmpshape2(k)
      enddo
      call a2(inshape,tmpshape1)
      call a2(tmpshape1,tmpshape2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=outshape(k)+0.5*tmpshape2(k)
      enddo
      return
      end




      subroutine shear2(inshape,outshape)
c generator for shear 2
      real inshape(msh),outshape(msh),tmpshape1(msh),tmpshape2(msh)

      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      call a2d(inshape,tmpshape1)
      call a1d(tmpshape1,tmpshape2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=tmpshape2(k)
      enddo
      call a2(inshape,tmpshape1)
      call a1(tmpshape1,tmpshape2)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=outshape(k)-tmpshape2(k)
      enddo
      return
      end


      subroutine shcmult(sh,c)
c multiply coefficients by constant c
      real sh(msh)
      n=sh(2)
      do k=3,(n+1)*(n+2)/2+2
         sh(k)=sh(k)*c
      enddo
      return
      end


      subroutine shear11(inshape,outshape)
      real inshape(msh),outshape(msh),t1(msh),t2(msh),t3(msh),t4(msh)

      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=nint(inshape(1))
      call a1da1d(inshape,t1)
      call a1da1d(t1,t2)
      call a2a2(t1,t3)
      call a2da2d(t1,t4)
      k=2
      do nn=0,n
          do j=0,nn
              i=nn-j
              k=k+1
              outshape(k)=-2*((i+1)**2+(j+1)**2)*inshape(k)+t2(k)+2*t3(k)-2*t4(k)
          enddo
      enddo
      call a1a1(inshape,t1)
      call a1a1(t1,t2)
      call a2a2(t1,t3)
      call a2da2d(t1,t4)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=outshape(k)+t2(k)-2*t3(k)+2*t4(k)
      enddo
      call a2a2(inshape,t1)
      call a2a2(t1,t2)
      call a2da2d(inshape,t3)
      call a2da2d(t3,t4)
      do k=3,2+(n+1)*(n+2)/2
          outshape(k)=(outshape(k)+t2(k)+t4(k)) / 4
      enddo
      return
      end


c Below, the ladder operators

c shape vectors are arranged as
c  beta
c  n
c  a00
c  a10, a01
c  a20, a11, a02
c ...
c below nn is i+j
c and k is a counter starting at 3 for a00

c=====================================================

      subroutine a1(inshape,outshape)

      real inshape(msh),outshape(msh)

      do i=1,msh
          outshape(i)=0
      enddo
      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      k=2
c effect of a1 is to move (nn) to the left, mult by sqrt(i)
      do nn=0,n
          do j=0,nn
              i=nn-j
              k=k+1
              if (i.gt.0) outshape(k-nn)=inshape(k)*sqrt(real(i))
          enddo
      enddo
      return
      end

c========================================================

      subroutine a2(inshape,outshape)

      real inshape(msh),outshape(msh)

      do i=1,msh
          outshape(i)=0
      enddo
      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      k=2
c effect of a2 is to move (nn+1) to the left, mult by sqrt(j)
      do nn=0,n
          do j=0,nn
              i=nn-j
              k=k+1
              if (j.gt.0) outshape(k-nn-1)=inshape(k)*sqrt(real(j))
          enddo
      enddo
      return
      end

c=====================================================

      subroutine a1d(inshape,outshape)

      real inshape(msh),outshape(msh)

      do i=1,msh
          outshape(i)=0
      enddo
      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      k=2
c effect of a1d is to move (nn+1) to the right, mult by sqrt(i+1)
      do nn=0,n-1
          do j=0,nn
              i=nn-j
              k=k+1
              outshape(k+nn+1)=inshape(k)*sqrt(real(i+1))
          enddo
      enddo
      return
      end

c========================================================

      subroutine a2d(inshape,outshape)

      real inshape(msh),outshape(msh)

      do i=1,msh
          outshape(i)=0
      enddo
      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      k=2
c effect of a2d is to move (nn+2) to the right, mult by sqrt(j+1)
      do nn=0,n-1
          do j=0,nn
              i=nn-j
              k=k+1
              outshape(k+nn+2)=inshape(k)*sqrt(real(j+1))
          enddo
      enddo
      return
      end


c=====================================================

      subroutine a1a1(inshape,outshape)

      real inshape(msh),outshape(msh)

      do i=1,msh
          outshape(i)=0
      enddo
      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      k=2
c effect of a1**2 is to move (2nn-1) to the left, mult by sqrt(i(i-1))
      do nn=0,n
          do j=0,nn
              i=nn-j
              k=k+1
              if (i.gt.1) outshape(k-2*nn+1)=inshape(k)*sqrt(real(i*(i-1)))
          enddo
      enddo
      return
      end

c========================================================

      subroutine a2a2(inshape,outshape)

      real inshape(msh),outshape(msh)

      do i=1,msh
          outshape(i)=0
      enddo
      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      k=2
c effect of a2**2 is to move (2nn+1) to the left, mult by sqrt(j(j-1))
      do nn=0,n
          do j=0,nn
              i=nn-j
              k=k+1
              if (j.gt.1) outshape(k-2*nn-1)=inshape(k)*sqrt(real(j*(j-1)))
          enddo
      enddo
      return
      end

c=====================================================

      subroutine a1da1d(inshape,outshape)

      real inshape(msh),outshape(msh)

      do i=1,msh
          outshape(i)=0
      enddo
      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      k=2
c effect of a1d**2 is to move (2nn+3) to the right, mult by sqrt((i+1)(i+2))
      do nn=0,n-2
          do j=0,nn
              i=nn-j
              k=k+1
              outshape(k+2*nn+3)=inshape(k)*sqrt(real((i+1)*(i+2)))
          enddo
      enddo
      return
      end

c========================================================

      subroutine a2da2d(inshape,outshape)

      real inshape(msh),outshape(msh)

      do i=1,msh
          outshape(i)=0
      enddo
      outshape(1)=inshape(1)
      outshape(2)=inshape(2)
      n=inshape(2)
      k=2
c effect of a2d**2 is to move (2nn+5) to the right, mult by sqrt((j+1)(j+2))
      do nn=0,n-2
          do j=0,nn
              i=nn-j
              k=k+1
              outshape(k+2*nn+5)=inshape(k)*sqrt(real((j+1)*(j+2)))
          enddo
      enddo
      return
      end


