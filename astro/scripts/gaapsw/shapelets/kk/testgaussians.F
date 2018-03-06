c convolve two asymmetric gaussians to make a symmetric one.

      real z1(41,41),z2(41,41),sh1(msh),sh2(msh),sh3(msh),sh4(msh),sh5(msh)

      call pgbeg(0,'/xs',3,2)
      sig=2.5
      sig2=sig*sig
      n=8
      do i=1,41
          x=i-21
          do j=1,41
              y1=j-23
              y2=j-19
              z1(i,j)=exp(-(x**2/0.8**2+y1**2/1.2**2)/2/sig2)/(2*3.1415926535*sig2)
              z2(i,j)=exp(-(y2**2/0.8**2+x**2/1.2**2)/2/sig2)/(2*3.1415926535*sig2)
          enddo
      enddo
      call encode(z1,41,41,41,21.,21.,20.,0.,sh1,n,sig)
      call writeshapen('shape1',sh1)
      call plotshape(sh1,40,'1')

      call encode(z2,41,41,41,21.,21.,20.,0.,sh2,n,sig)
      call writeshapen('shape2',sh2)
      call plotshape(sh2,40,'2')

      call convolve(sh1,sh2,sqrt(2.)*sig,sh3)
      call writeshapen('convolution',sh3)
      call plotshape(sh3,40,'1 X 2')

      call shear1(sh3,sh4)
      do k=3,2+(n+1)*(n+2)/2
          sh4(k)=sh3(k)+0.05*sh4(k)
      enddo
      call shear1(sh4,sh5)
      do k=3,2+(n+1)*(n+2)/2
          sh4(k)=sh4(k)+0.05*sh5(k)
      enddo
      call writeshapen('sheared',sh4)
      call plotshape(sh4,40,'sheared')

      call shear1(sh4,sh5)
      do k=3,2+(n+1)*(n+2)/2
          sh4(k)=sh4(k)-0.05*sh5(k)
      enddo
      call shear1(sh4,sh5)
      do k=3,2+(n+1)*(n+2)/2
          sh5(k)=sh4(k)-0.05*sh5(k)
      enddo
      call writeshapen('unsheared',sh5)
      call plotshape(sh5,40,'unsheared')

      call pgend

      end
