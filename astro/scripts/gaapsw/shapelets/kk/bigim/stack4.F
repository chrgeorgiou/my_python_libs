      character*20 imname
      dimension pix(256,256),out(128,128)
      read(*,'(a20)') imname
      write(*,*) imname
      call get2(imname,nx,ny,pix,256,256)
      do j=1,128
         do i=1,128
            out(i,j)=(pix(i,j)+pix(i+128,j)+pix(i,j+128)+pix(i+128,j+128))/4
         enddo
      enddo
      do j=1,128
         do i=1,128
            out(129-i,129-j)=out(129-i,129-j)+(pix(i,j)+pix(i+128,j)+pix(i,j+128)+pix(i+128,j+128))/4
         enddo
      enddo
      imname(1:1)='s'
      write(*,*) imname
      call put2(imname,128,128,out,128,128)
      end

