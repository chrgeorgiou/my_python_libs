      parameter (pi=3.1415926535)

      real sh(msh)
      real zz(801,801)
      real tr(6)
      data tr /-401.,1.,0.,-401.,0.,1./

c     calculate a series of sheared shapelet basis functions; tabulate
c     their shapelet expansion

c      call pgbeg(0,'/xs',4,2)
      n=8
      do ic=0,8,2
         do ig1=-70,70,2
            do ig2=-70,70,2
               if (ig1**2+ig2**2 .gt.70**2) goto 10
               g1=0.01*ig1
               g2=0.01*ig2
               do ix=-400,400
                  x=ix*0.02
                  do iy=-400,400
                     y=iy*0.02
c normalize the map to unit determinant - no magnification, just shear
                     xs=(x*(1-g1)-y*g2)/sqrt(1-g1*g1-g2*g2)
                     ys=(-x*g2+(1+g1)*y)/sqrt(1-g1*g1-g2*g2)
                     zz(401+ix,401+iy)=shround(xs,ys,ic,1.)
                  enddo
               enddo
c plot contours of a subset if desired
c               if (mod(ig,10).eq.0) then
c                   call pgenv(-400.,400.,-400.,400.,1,0)
c                   do icont=-10,10
c                       call pgcont(zz,801,801,1,801,1,801,0.1*icont,1,tr)
c                   enddo
c               endif
               call fitshapelets(zz,801,801,801,401.,401.,400.,1e-20,sh,n,50.)
               write(*,1000) ic,g1,g2,(sh(k+2)/50.,k=1,(n+1)*(n+2)/2)
 1000          format (i4,2f10.6,msh(f15.8))
   10      enddo
         enddo
      enddo


      end







