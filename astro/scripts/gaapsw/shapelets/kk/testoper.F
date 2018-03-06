      real shape(msh),shape2(msh),shape3(msh),shape4(msh)
      real z(30,30),tr(6),contlev(20)

      beta=2
      n=8
      nn=(n+1)*(n+2)/2
      shape(1)=beta
      shape(2)=n
      do i=3,2+nn
          shape(i)=0
      enddo
      shape(3)=1
c      shape(4)=0.2
c      shape(10)=0.1
c      shape(20)=0.3
      call pgbeg(0,'/xs',3,3)
      tr(1)=-15.5
      tr(2)=1.
      tr(3)=0
      tr(4)=-15.5
      tr(5)=0
      tr(6)=1
      do i=1,20
          contlev(i)=0.05*i
      enddo

      call writeshape('Shape:',shape)
      call decode(shape,z,30,30)
      call pgenv(-6.,6.,-6.,6.,1,1)
      call pgcons(z,30,30,1,30,1,30,contlev,20,tr)
      
      write(0,*)
      write(0,*) 'TRANS1'
      call trans1(shape,shape2)
      shape3(1)=shape(1)
      shape3(2)=shape(2)
      do i=3,2+nn
          shape3(i)=shape(i)+0.05*shape2(i)
      enddo
      call writeshape('Shape:',shape3)
      call decode(shape3,z,30,30)
      call pgenv(-6.,6.,-6.,6.,1,1)
      call pgcons(z,30,30,1,30,1,30,contlev,20,tr)

      write(0,*)
      write(0,*) 'TRANS2'
      call trans2(shape,shape2)
      shape3(1)=shape(1)
      shape3(2)=shape(2)
      do i=3,2+nn
          shape3(i)=shape(i)+0.05*shape2(i)
      enddo
      call writeshape('Shape:',shape3)
      call decode(shape3,z,30,30)
      call pgenv(-6.,6.,-6.,6.,1,1)
      call pgcons(z,30,30,1,30,1,30,contlev,20,tr)

      write(0,*)
      write(0,*) 'CONV'
      call conv(shape,shape2)
      shape3(1)=shape(1)
      shape3(2)=shape(2)
      do i=3,2+nn
          shape3(i)=shape(i)+0.05*shape2(i)
      enddo
      call writeshape('Shape:',shape3)
      call decode(shape3,z,30,30)
      call pgenv(-6.,6.,-6.,6.,1,1)
      call pgcons(z,30,30,1,30,1,30,contlev,20,tr)

      write(0,*)
      write(0,*) 'ROT'
      call rot(shape,shape2)
      shape3(1)=shape(1)
      shape3(2)=shape(2)
      do i=3,2+nn
          shape3(i)=shape(i)+0.05*shape2(i)
      enddo
      call writeshape('Shape:',shape3)
      call decode(shape3,z,30,30)
      call pgenv(-6.,6.,-6.,6.,1,1)
      call pgcons(z,30,30,1,30,1,30,contlev,20,tr)

      write(0,*)
      write(0,*) 'SHEAR1'
      call shear1(shape,shape2)
      shape3(1)=shape(1)
      shape3(2)=shape(2)
      do i=3,2+nn
          shape3(i)=shape(i)+0.05*shape2(i)
      enddo
      call writeshape('Shape:',shape3)
      call decode(shape3,z,30,30)
      call pgenv(-6.,6.,-6.,6.,1,1)
      call pgcons(z,30,30,1,30,1,30,contlev,20,tr)

      write(0,*)
      write(0,*) 'SHEAR2'
      call shear2(shape,shape2)
      shape3(1)=shape(1)
      shape3(2)=shape(2)
      do i=3,2+nn
          shape3(i)=shape(i)+0.05*shape2(i)
      enddo
      call writeshape('Shape:',shape3)
      call decode(shape3,z,30,30)
      call pgenv(-6.,6.,-6.,6.,1,1)
      call pgcons(z,30,30,1,30,1,30,contlev,20,tr)

      write(0,*)
      write(0,*) 'SHAPE'
      call writeshape('Shape:',shape)
      call decode(shape,z,30,30)
      call pgenv(-6.,6.,-6.,6.,1,1)
      call pgcons(z,30,30,1,30,1,30,contlev,20,tr)

      write(0,*)
      write(0,*) 'SHEAR11'
      call shear11(shape,shape2)
      call writeshape('Shear11:',shape2)
      call decode(shape2,z,30,30)
      call pgenv(-6.,6.,-6.,6.,1,1)
      call pgcons(z,30,30,1,30,1,30,contlev,20,tr)

      write(0,*)
      write(0,*) 'SHEAR11'
      call shear1(shape,shape3)
      call shear1(shape3,shape2)
      call writeshape('Shear11:',shape2)
      call decode(shape2,z,30,30)
      call pgenv(-6.,6.,-6.,6.,1,1)
      call pgcons(z,30,30,1,30,1,30,contlev,20,tr)


      call a1(shape,shape2)
      call a1(shape2,shape3)
      call a1a1(shape,shape4)
      call writeshape('a1 a1:',shape3)
      call writeshape('a1a1:',shape4)

      call a1d(shape,shape2)
      call a1d(shape2,shape3)
      call a1da1d(shape,shape4)
      call writeshape('a1d a1d:',shape3)
      call writeshape('a1da1d:',shape4)

      call a2(shape,shape2)
      call a2(shape2,shape3)
      call a2a2(shape,shape4)
      call writeshape('a2 a2:',shape3)
      call writeshape('a2a2:',shape4)

      call a2d(shape,shape2)
      call a2d(shape2,shape3)
      call a2da2d(shape,shape4)
      call writeshape('a2d a2d:',shape3)
      call writeshape('a2da2d:',shape4)



      call pgend
      end
