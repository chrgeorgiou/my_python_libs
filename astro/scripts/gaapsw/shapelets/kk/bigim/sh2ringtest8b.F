      dimension pix(768*2,768*2),sh(msh)
      character*21 filename
      integer*8 id

c read in a shapelet expansion from stdin and turn it into an image
c ringtest.fits which has four images of the same source, at eight
c position angles, sheared by the amount specified in ringtest.in


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
c 8 panels with rotated, sheared galaxies
      do iang=0,7
          cosang=cos(iang*3.1415926535/8.)
          sinang=sin(iang*3.1415926535/8.)
          x0=mod(iang,3)*512+1
          y0=(iang/3)*512+1
          xc=255.5+x0
          yc=255.5+y0
          do iy=y0,y0+511
              dy=scale*(iy-yc)
              do ix=x0,x0+511
                  dx=scale*(ix-xc)
                  xs=(1-g1)*dx-g2*dy
                  ys=-g2*dx+(1+g1)*dy
                  xr=cosang*xs-sinang*ys
                  yr=sinang*xs+cosang*ys
                  pix(ix,iy)=shapeletpix(xr,yr,sh)
              enddo
          enddo
      enddo
c add a delta fn at center of 9th panel, so later see PSF in convolved image.
      pix(1279,1279)=shsum

      write(filename,1001) id
 1001 format("pt",i9.9,".fits")
      call put2(filename,768*2,768*2,pix,768*2,768*2)
      
      goto 5
 10   continue
      end
