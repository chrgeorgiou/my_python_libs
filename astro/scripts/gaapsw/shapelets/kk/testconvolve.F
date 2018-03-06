      real shape1(msh),shape2(msh),shape3(msh),shape4(msh),z(50,50),tr(6)

      do i=1,msh
          shape1(i)=0
          shape2(i)=0
          shape3(i)=0
      enddo

      call pgbegin(0,'/xs',2,2)

      s1=3
      s2=4
      beta1=5
      beta2=5
      betaout=5.

      do i=1,50
         x=i-25.5
         do j=1,50
            y=j-25.5
            z(i,j)=exp(-0.5*((x-1)**2+(y+1)**2)/s1**2)/(2*3.1416*s1*s1)
         enddo
      enddo
      call fitshapelets(z,50,50,50,25.5,25.5,20.,0.,shape1,10,beta1)
      call plotshape(shape1,40,'Shape 1')
      call writeshapen('shape1',shape1)

      do i=1,50
         x=i-25.5
         do j=1,50
            y=j-25.5
            z(i,j)=exp(-0.5*((x+1)**2+(y-1)**2)/s2**2)/(2*3.1416*s2*s2)
         enddo
      enddo
      call fitshapelets(z,50,50,50,25.5,25.5,20.,0.,shape2,10,beta2)
      call plotshape(shape2,40,'Shape 2')
      call writeshapen('shape2',shape2)
      write(0,*)

      call convolve(shape1,shape2,betaout,shape3)
      call writeshapen('conv',shape3)
      call plotshape(shape3,40,'Shape 3')

      call shear1(shape3,shape4)
      n=shape3(2)
      do k=3,(n+1)*(n+2)/2
         shape4(k)=shape3(k)+0.1*shape4(k)
      enddo
      call plotshape(shape4,40,'Shape 4')

      call pgend

      end
