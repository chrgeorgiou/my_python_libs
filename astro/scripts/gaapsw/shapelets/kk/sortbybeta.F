      dimension sh(msh-2,mmcat),beta(mmcat),nn(mmcat)
      character str1(mmcat)*91,str2(mmcat)*74
      real indx(mmcat)
      bmin=1e10
      bmax=-bmin
      do i=1,mmcat
    1     read(*,*,err=1,end=2) x,y,f,fe,fw,a,b,th,id,iflg
          write(str1(i),1001) x,y,f,fe,fw,a,b,th,id,iflg
 1001     format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i10,i5)
          read(*,*) sum,dx,dy,bg,fitradius,sherr
          write(str2(i),1002) sum,dx,dy,bg,fitradius,sherr
 1002     format (g14.7,2f10.5,g15.7,f10.5,g15.6)
          read(*,*) beta(i),nn(i)
          nsh=(nn(i)+1)*(nn(i)+2)/2
          read(*,*) (sh(k,i),k=1,nsh)
          read(*,*)
      enddo
    2 nread=i-1
      write(0,*) 'Read ',nread,' shapelets'
      do i=1,nread
          indx(i)=i
      enddo
      call sort2(nread,beta,indx)
      do i=1,nread
          ii=indx(i)
          write(*,'(a91)') str1(ii)
          write(*,'(a74)') str2(ii)
          write(*,*) beta(i),nn(ii)
          nsh=(nn(ii)+1)*(nn(ii)+2)/2
          write(*,*) (sh(k,ii),k=1,nsh)
          write(*,*)
      enddo
      end

