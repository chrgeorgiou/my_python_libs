      real sh1(msh),sh2(msh),sh3(msh),sh4(msh),sh5(msh),sh6(msh),sh7(msh),sh8(msh)


      call pgbeg(0,'/xs',3,3)
      g1=0.3
      n=16
      nsh=(n+1)*(n+2)/2
      sh1(1)=3.
      sh1(2)=n

c shear circular shapelets in turn; calculate to various orders.

      do m=0,n-2,2
c make circular shapelet, of order m

      do k=3,nsh+2
          sh1(k)=0
      enddo
      k0=2+m*(m+1)/2 + 1
      cc=1
      sh1(k0)=cc
      do kk=2,m,2
          cc=cc*(m/2-kk/2+1)/(kk/2)*sqrt((kk-1)*kk /real((m-kk+2)*(m-kk+1)))
          sh1(k0+kk)=cc
      enddo
      sh3(1)=sh1(1)
      sh3(2)=sh1(2)
      sh4(1)=sh1(1)
      sh4(2)=sh1(2)
      sh5(1)=sh1(1)
      sh5(2)=sh1(2)
      sh6(1)=sh1(1)
      sh6(2)=sh1(2)
      sh7(1)=sh1(1)
      sh7(2)=sh1(2)
      sh8(1)=sh1(1)
      sh8(2)=sh1(2)
c     normalize the shape
      sum=seriessum(sh1)
      do kk=0,m,2
          sh1(k0+kk)=sh1(k0+kk)/sum
      enddo
      call plotshape(sh1,40,'ORDER 0')
      call shear1(sh1,sh2)
      do k=3,nsh+2
          sh3(k)=sh1(k)+g1*sh2(k)
      enddo
      call plotshape(sh3,40,'ORDER 1')
      call writeshapen('ORDER 1',sh3)
      call shear1(sh2,sh3)
      do k=3,nsh+2
          sh4(k)=sh1(k)+g1*sh2(k)+g1*g1/2*sh3(k)
      enddo
      call plotshape(sh4,40,'ORDER 2')
      call writeshapen('ORDER 2',sh4)
      call shear1(sh3,sh4)
      do k=3,nsh+2
          sh5(k)=sh1(k)+g1*sh2(k)+g1*g1/2*sh3(k)+g1*g1*g1/6*sh4(k)
      enddo
      call plotshape(sh5,40,'ORDER 3')
      call writeshapen('ORDER 3',sh5)
      call shear1(sh4,sh5)
      do k=3,nsh+2
          sh6(k)=sh1(k)+g1*sh2(k)+g1*g1/2*sh3(k)+g1*g1*g1/6*sh4(k)+g1**4/24*sh5(k)
      enddo
      call plotshape(sh6,40,'ORDER 4')
      call writeshapen('ORDER 4',sh6)
      call shear1(sh5,sh6)
      do k=3,nsh+2
          sh7(k)=sh1(k)+g1*sh2(k)+g1*g1/2*sh3(k)+g1*g1*g1/6*sh4(k)+g1**4/24*sh5(k)+g1**5/120*sh6(k)
      enddo
      call plotshape(sh7,40,'ORDER 5')
      call writeshapen('ORDER 5',sh7)
      call shear1(sh6,sh7)
      do k=3,nsh+2
          sh8(k)=sh1(k)+g1*sh2(k)+g1*g1/2*sh3(k)+g1*g1*g1/6*sh4(k)+g1**4/24*sh5(k)+g1**5/120*sh6(k)+g1**6/720*sh7(k)
      enddo
      call plotshape(sh8,40,'ORDER 6')
      call writeshapen('ORDER 6',sh8)
      call pgadvance
      call pgadvance
      enddo
      call pgend
      end

      
