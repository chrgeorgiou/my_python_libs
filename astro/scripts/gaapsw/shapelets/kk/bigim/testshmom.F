      real sh(msh),z(101,101)
      beta=2.7
      sh(1)=beta
      sh(2)=n
      do k=3,(n+1)*(n+2)/2
         sh(k)=0
      enddo
      sh(3)=1
      call decode (sh,z,101,101)
      write(0,*) 
      write(0,*) '00:'
      call moms(z,beta)

      sh(3)=0
      sh(4)=1
      call decode (sh,z,101,101)
      write(0,*) 
      write(0,*) '10:'
      call moms(z,beta)

      sh(4)=0
      sh(5)=1
      call decode (sh,z,101,101)
      write(0,*) 
      write(0,*) '01:'
      call moms(z,beta)

      sh(5)=0
      sh(6)=1
      call decode (sh,z,101,101)
      write(0,*) 
      write(0,*) '20:'
      call moms(z,beta)

      sh(6)=0
      sh(7)=1
      call decode (sh,z,101,101)
      write(0,*) 
      write(0,*) '11:'
      call moms(z,beta)

      sh(7)=0
      sh(8)=1
      call decode (sh,z,101,101)
      write(0,*) 
      write(0,*) '02:'
      call moms(z,beta)

      end

      subroutine moms(z,beta)
      parameter (pi=3.1415926536)
      real z(101,101)
      sum00=0
      sum10=0
      sum01=0
      sum20=0
      sum11=0
      sum02=0
      do ix=1,101
         do iy=1,101
            x=ix-51
            y=iy-51
            w=exp(-0.5*(x**2+y**2)/beta**2)
            sum00=sum00+w*z(ix,iy)
            sum10=sum10+x*w*z(ix,iy)
            sum01=sum01+y*w*z(ix,iy)
            sum20=sum20+x*x*w*z(ix,iy)
            sum11=sum11+x*y*w*z(ix,iy)
            sum02=sum02+y*y*w*z(ix,iy)
         enddo
      enddo
      write(0,*) sum00/beta/sqrt(pi)
      write(0,*) sum10/beta**2/sqrt(pi),sum01/beta**2/sqrt(pi)
      write(0,*) sum20/beta**3/sqrt(pi),sum11/beta**3/sqrt(pi),sum02/beta**3/sqrt(pi)
      return
      end

            
