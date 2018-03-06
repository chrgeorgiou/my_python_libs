      dimension sh(msh-2)
      character str1*90
      bmin=1e10
      bmax=-bmin
      do i=1,mmcat
    1     read(*,'(a90)',err=1,end=2) str1
          read(*,*) sum,dx,dy,bg,fitradius,sherr
          sherr=sherr*sum
          read(*,*) beta,n
          nsh=(n+1)*(n+2)/2
          read(*,*) (sh(k),k=1,nsh)
          read(*,*)

          write(*,1003) str1,sum,dx,dy,bg,fitradius,sherr,beta,n
 1003     format (a90,g14.7,2f10.5,g15.7,f10.5,g15.6,f10.5,i4)
      enddo
    2 nread=i-1
      write(0,*) 'Read ',nread,' shapelets'
      end

