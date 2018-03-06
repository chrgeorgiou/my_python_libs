      real sh(msh),sh2(msh),z(40,40)

      call pgbeg(0,'/xs',3,3)

      sh(1)=2.4
      sh(2)=6
      n=sh(2)
      do k=3,msh
          sh(k)=0
      enddo

      sh(3)=1
      sh(7)=2
      call plotshape(sh,22,'test')
      sh(5)=2
      call plotshape(sh,22,'test')
      sh(6)=-1
      call plotshape(sh,22,'test')
      call decode(sh,z,40,40)
      call encode(z,40,40,40,20.5,20.5,20.,0.,sh2,6,3.)
      call plotshape(sh2,22,'test')

      do k=3,msh
          sh(k)=0
      enddo
      do k=3,2+(n+1)*(n+2)/2
          sh(k)=1
          call plotshape(sh,50,'test')
          sh(k)=0
      enddo

      call pgend

      end
