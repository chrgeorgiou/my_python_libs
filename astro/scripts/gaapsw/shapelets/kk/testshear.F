      real s1(msh),s2(msh),s3(msh),s4(msh),s5(msh)

      beta=99
      n=4
      s1(1)=beta
      s1(2)=n
      s1(3)=1
      k=3
      do nn=1,n
          do j=0,nn
              i=nn-j
              k=k+1
              if (nn.eq.n) then
                  s1(k)=0
              else 
                  s1(k)=10*i+j
              endif
          enddo
      enddo
      ktot=k
      call writeshape('s1',s1)
      
      call shear1(s1,s2)
      write(0,*) 'shear1:'
      call writeshape('s2',s2)
      call shear2(s1,s2)
      write(0,*) 'shear2:'
      call writeshape('s2',s2)
      call rot(s1,s2)
      write(0,*) 'rot:'
      call writeshape('s2',s2)
      call conv(s1,s2)
      write(0,*) 'conv:'
      call writeshape('s2',s2)
      call trans1(s1,s2)
      write(0,*) 'trans1:'
      call writeshape('s2',s2)
      call trans2(s1,s2)
      write(0,*) 'trans2:'
      call writeshape('s2',s2)

      end
