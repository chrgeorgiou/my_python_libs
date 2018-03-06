      dimension pix(256,256),sh(msh)
      character*21 filename
      integer*8 id

      open(15,file='ringtest.in',status='old')
      read(15,*) g1,g2,scale
      close(15)
      write(0,*) 'Ring test with shears ',g1,g2
      write(0,*) 'Output pixel size is original times ',scale

 5    continue
      read(*,*,end=10) x,y,f,fe,fw,a,b,pa,id,iflg
      read(*,*) shsum,dx,dy,bg,fitradius,sherr
      read(*,*) beta,n
      read(*,*) (sh(k),k=3,(n+1)*(n+2)/2+2)
      sh(1)=beta
      sh(2)=n
      do k=3,(n+1)*(n+2)/2+2
          sh(k)=sh(k)*shsum
      enddo
c shear by g1,g2 in bottom left
      do iy=1,128
          dy=scale*(iy-64.5)
          do ix=1,128
              dx=scale*(ix-64.5)
              xs=(1-g1)*dx-g2*dy
              ys=-g2*dx+(1+g1)*dy
              pix(ix,iy)=shapeletpix(xs,ys,sh)
          enddo
      enddo
c shear by g1,g2 after 45 deg rot in bottom right
      do iy=1,128
          dy=scale*(iy-64.5)
          do ix=129,256
              dx=scale*(ix-192.5)
              xs=(1-g1)*dx-g2*dy
              ys=-g2*dx+(1+g1)*dy
              xr=sqrt(0.5)*(xs-ys)
              yr=sqrt(0.5)*(xs+ys)
              pix(ix,iy)=shapeletpix(xr,yr,sh)
          enddo
      enddo
c shear by g1,g2 after 90 deg rot in top left
      do iy=129,256
          dy=scale*(iy-192.5)
          do ix=1,128
              dx=scale*(ix-64.5)
              xs=(1-g1)*dx-g2*dy
              ys=-g2*dx+(1+g1)*dy
              xr=-ys
              yr=xs
              pix(ix,iy)=shapeletpix(xr,yr,sh)
          enddo
      enddo
c shear by g1,g2 after -45 deg rot in top right
      do iy=129,256
          dy=scale*(iy-192.5)
          do ix=129,256
              dx=scale*(ix-192.5)
              xs=(1-g1)*dx-g2*dy
              ys=-g2*dx+(1+g1)*dy
              xr=sqrt(0.5)*(xs+ys)
              yr=sqrt(0.5)*(-xs+ys)
              pix(ix,iy)=shapeletpix(xr,yr,sh)
          enddo
      enddo

      write(filename,1000) id
 1000 format("rt",i9.9,".fits")
      call put2(filename,256,256,pix,256,256)

      goto 5
 10   continue
      end
