      real sh(msh),z(101,101),z2(101,101),xplot(101)

      do i=1,101
          xplot(i)=i
      enddo
      call pgbeg(0,'?',3,3)
      rhalf=4.
      do nser=1,4
          xser=nser
c lay down a sersic profile in the z array. Scale radius 2.5 pixels.
          do ig1=0,0
              g1=ig1*0.1
              g2=0
              do ix=1,101
                  x=ix-51.
                  do iy=1,101
                      y=iy-51.
                      xsh= x*(1-g1)-y*g2
                      ysh=-g2*x+(1+g1)*y
                      r=sqrt(xsh*xsh+ysh*ysh)
                      z(ix,iy)=exp(-log(2.)*(r/rhalf)**(1./xser))
                  enddo
              enddo
              beta=rhalf/1.15*1.5
              call fitshapelets(z,101,101,101,51.,51.,4*beta,1e-10,sh,12,beta)
              call writeshapen(' ',sh)
              call plotshape(sh,101,'SHAPELET FIT')
              call decode(sh,z2,101,101)
              call pgenv(51.,101.,-0.1,1.1,0,0)
              call pgline(101,xplot,z(1,51))
              call pgline(101,xplot,z2(1,51))
          enddo
      enddo
      call pgend
      end
