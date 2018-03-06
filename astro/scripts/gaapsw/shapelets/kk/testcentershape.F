      real sh(msh),sh2(msh)

      n=8
      nsh=(n+1)*(n+2)/2+2
      sh(1)=2.3
      sh(2)=n
      do k=3,nsh
         sh(k)=0
      enddo
      sh(3)=2
      sh(4)=0.1
      sh(5)=-0.2
      sh(6)=0.3
      sh(7)=0.4
      sh(5)=0.5
      call centershape(sh,sh2,dx,dy)
      write(*,*) 'Centershape shifted by dx,dy =',dx,dy,' pixels.'
      call writeshapen('ORIG',sh)
      call writeshapen('AFTER CENTERSHAPE',sh2)
      end
