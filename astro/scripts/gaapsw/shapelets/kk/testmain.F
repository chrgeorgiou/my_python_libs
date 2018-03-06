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
      write(0,1000) (s1(k),k=1,ktot)
      
      call a1(s1,s2)
      write(0,*) 'a1:'
      write(0,1000) (s2(k),k=1,ktot)
      call a1(s1,s2)
      write(0,*) 'a1:'
      write(0,1000) (s2(k),k=1,ktot)
      call a2(s1,s2)
      write(0,*) 'a1d:'
      write(0,1000) (s2(k),k=1,ktot)
      call a1d(s1,s2)
      write(0,*) 'a2:'
      write(0,1000) (s2(k),k=1,ktot)
      call a2d(s1,s2)
      write(0,*) 'a2d:'
      write(0,1000) (s2(k),k=1,ktot)

      call a1d(s1,s2)
      call a1(s2,s3)
      call a1(s1,s2)
      call a1d(s2,s4)
      write(0,*) 'a1 a1d - a1d a1 (=1):'
      write(0,1000) (s3(k)-s4(k),k=1,ktot)

      call a1(s1,s2)
      call a1d(s2,s3)
      write(0,*) 'a1d a1 (=N1):'
      write(0,1000) (s3(k),k=1,ktot)
1000  format (2f6.2,2x,f7.2,2x,2f7.2,2x,3f7.2,2x,4f7.2,2x,5f7.2)

      end
