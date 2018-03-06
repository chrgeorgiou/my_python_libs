      real z(30,30),shape(msh),shape2(msh)
      real tr(6)

      beta=1.
      beta2=beta
      n=8
      nsh=(n+1)*(n+2)/2
      shape(1)=beta
      shape(2)=n
      k=2
      do nn=0,n
          do j=0,nn
              i=nn-j
              k=k+1
              shape(k)=0
          enddo
      enddo
      shape(3)=1.
      do k=4,(n+1)*(n+2)/2+2
          shape(k)=(2*ran3(idum)-1) * 0.2
      enddo
      call writeshape('shape',shape)
      call decode(shape,z,30,25)


      call pgbeg(0,'/xs',2,2)

      call fitshapelets(z,30,25,25,13.,13.,10.,0.,shape2,n,beta2)
      call plotshape(shape2,25,'Shape 2')
      call writeshape('shape2',shape2)
      do k=3,nsh+2
          shape2(k)=shape2(k)-shape(k)
      enddo
      call writeshape('Diff (fit)',shape2)

      call encode(z,30,25,25,13.,13.,10.,0.,shape2,n,beta2)
      call plotshape(shape2,25,'Shape 2')
      call writeshape('shape2',shape2)
      do k=3,nsh+2
          shape2(k)=shape2(k)-shape(k)
      enddo
      call writeshape('Diff (sum)',shape2)
      call pgend

      end

