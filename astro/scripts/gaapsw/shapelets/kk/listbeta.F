      dimension sh(msh-2,mmcat)
      bmin=1e10
      bmax=-bmin
      do i=1,mmcat
    1     read(*,*,err=1,end=2) x,y,f,fe,fw,a,b,pa,id,iflg
          read(*,*) sum,dx,dy,bg,fitradius,sherr
          read(*,*) beta,nn
          nsh=(nn+1)*(nn+2)/2
          read(*,*) (sh(k,i),k=1,nsh)
          read(*,*)
          write(*,*) beta,fw,f,fe,sum,sherr*sum,iflg
      enddo
    2 nread=i-1
      write(0,*) 'Read ',nread,' shapelets'
      end

