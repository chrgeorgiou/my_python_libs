c make an image of the same dimensions as inimage.fits
c pixel values are zero where there were no sources
c where there was a shapelet, insert this with its associated background.
c Result should be something that can be compared with inimage.fits
c The residual should be blank sky where there was no source fitted; or 
c fit error where there was.

      real pix(mmpix,mmpix),sh(msh),bgarray(mmcat)

      call get2('inimage.fits',nx,ny,pix,mmpix,mmpix)
      do i=1,nx
         do j=1,ny
            pix(i,j)=0
         enddo
      enddo

      nread=0
      do i=1,mmcat
    1     read(*,*,err=1,end=99) x,y,f,fe,fw,a,b,th,id,iflg
          read(*,*) sum,dx,dy,bg,fitradius,sherr
          read(*,*) beta,n
          read(*,*) (sh(k),k=3,(n+1)*(n+2)/2+2)
          read(*,*)
          sh(1)=beta
          sh(2)=n
          do k=3,(n+1)*(n+2)/2+2
             sh(k)=sh(k)*sum
          enddo
c          call putinimage(sh,x+dx,y+dy,fitradius,bg,pix,nx,ny)
          call addtoimage(sh,x-dx,y-dy,1.5*fitradius,pix,mmpix,mmpix)
          nread=nread+1
          bgarray(nread)=bg
       enddo
       write(0,*) 'Reconstructimage stopped after ',mmcat,' sources'
 99    continue
c add in median background value
       bg=fmedian(nread,bgarray)
       write(0,*) 'Median background value was ',bg
       do i=1,nx
           do j=1,ny
               pix(i,j)=pix(i,j)+bg
           enddo
       enddo
       write(0,*) 'Reconstructimage inserted ',nread,' sources'
       call put2('outimage.fits',nx,ny,pix,mmpix,mmpix)
       write(0,*) 'Reconstructed image is outimage.fits'
       end
