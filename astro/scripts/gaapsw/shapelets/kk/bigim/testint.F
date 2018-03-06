      parameter (pi=3.1415926535)
      real imn(0:20,0:20)
      double precision sum
      data imn /441*0./

c calculate the integral of the product of:
c    Hm(x/beta) Hn(y/beta) exp(-r^2/2beta^2)   with    an elliptical Gaussian aperture function, maj,min axis a,b, PA from y to x

      write(*,*) 'm,n,a,b,pa,beta:?'
      read(*,*) m,n,a,b,pa,beta

      parad=pa*pi/180.
      biga2=1./(1./a**2+1./beta**2)
      bigb2=1./(1./b**2+1./beta**2)

c coefficients for recurrence I_mn from I_m-2,n and I_m-1,n-1
      coefmmin2n=4*biga2/beta**2*sin(parad)**2-2*bigb2/beta**2*sin(2*parad)-2
      coefmmin2n=4*biga2/beta**2-2
      coefmmin2n=4*biga2/beta**2*sin(parad)**2+4*bigb2/beta**2*cos(parad)**2-2
      coefmmin1nmin1=2*(biga2-bigb2)/beta**2*sin(2*parad)
c coefficient for recurrence I_0n from I_0,n-2 (obtained from coefmmin2n by swapping biga and bigb, and PA-->90-PA_
      coefmnmin2=4*bigb2/beta**2*cos(parad)**2-2*biga2/beta**2*sin(2*parad)-2
      coefmnmin2=4*bigb2/beta**2-2
      coefmnmin2=4*biga2/beta**2*cos(parad)**2+4*bigb2/beta**2*sin(parad)**2-2

c seed the recurrence. I_00 is simple 2D elliptical Gaussian integral
      imn(0,0)=2*pi*sqrt(biga2*bigb2)
c work up the m=0 axis
      do in=2,2*(n/2)
         imn(0,in)=(in-1)*coefmnmin2*imn(0,in-2)
      enddo
c and fill in the remainder of the table:
c first do each row from left to right, then work up to row above
      do in=0,n
         do im=1,m-(n-in)
            if (im.gt.1 .and. in.gt.0) then
               imn(im,in)=(im-1)*coefmmin2n*imn(im-2,in) + in*coefmmin1nmin1*imn(im-1,in-1)
            elseif (im.gt.1) then
               imn(im,in)=(im-1)*coefmmin2n*imn(im-2,in)
            elseif (in.gt.0) then
               imn(im,in)=in*coefmmin1nmin1*imn(im-1,in-1)
            endif
            write(*,*) im,in
         enddo
      enddo

      do in=n,0,-1
         write(*,'(15(i3,i3,g10.3))') (im,in,imn(im,in),im=0,m-(n-in))
      enddo

      write(*,*) imn(m,n)

      nrand=10000000
      do i=1,nrand
         xx=gasdev(idum)*a
         yy=gasdev(idum)*b
         x=xx*sin(parad)-yy*cos(parad)
         y=xx*cos(parad)+yy*sin(parad)
         hx=herm(m,x/beta)
         hy=herm(n,y/beta)
         sum=sum+hx*hy*exp(-(x*x+y*y)/(2*beta**2))
      enddo
      sum=sum/nrand * 2*pi*a*b
      write(*,*) 'Monte carlo integral: ',sum

      end




      function herm(n,x)
      if (n.le.0) then
         herm=1
      elseif (n.eq.1) then
         herm=2*x
      else
         tmin2=1
         tmin1=2*x
         do i=2,n
            t=2*x*tmin1-2*(i-1)*tmin2
            tmin2=tmin1
            tmin1=t
         enddo
         herm=t
      endif
      return
      end

cccccccccccc

c      do i=-10,10
c         x=0.1*i
c         write(*,*) x,1.,2*x,4*x**2-2,8*x**3-12*x,16*x**4-48*x**2+12
c         write(*,*) x,(herm(n,x),n=0,4)
c      enddo


