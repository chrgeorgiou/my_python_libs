      real z(mmpix,mmpix),zerr(mmpix,mmpix)
      parameter (satur=50000.)

c read parameters for PSF stars, and shape of the galaxy
      read(*,*) idum,nx,ny,fs,pa,pb,ppa,ns,fg,gr,g1,g2,ng,bg
      write(0,*) 'Image size to make, bg: ',nx,ny,bg
      write(0,*) 'Stars: flux, rms maj, rms min, PA, number:',
     :    fs,pa,pb,ppa,ns
      write(0,*) 'Galaxies: flux, rms radius, g1, g2, number:',
     :    fg,gr,g1,g2,ng
      do i=1,mmpix
          do j=1,mmpix
              z(i,j)=bg
          enddo
      enddo
      cpa=cos(ppa*3.14159265/180)
      spa=sin(ppa*3.14159265/180)
      covsxx=pa*pa*cpa**2+pb*pb*spa**2
      covsxy=(pa*pa-pb*pb)*cpa*spa
      covsyy=pa*pa*spa**2+pb*pb*cpa**2
      dets=covsxx*covsyy-covsxy**2
      cinvsxx=covsyy/dets
      cinvsxy=-covsxy/dets
      cinvsyy=covsxx/dets
      covgxx=((1+g1)**2+g2**2) * gr**2
      covgxy=2*g2 * gr**2
      covgyy=((1-g1)**2+g2**2) * gr**2
      covxx=covsxx+covgxx
      covxy=covsxy+covgxy
      covyy=covsyy+covgyy
      det=covxx*covyy-covxy**2
      cinvxx= covyy/det
      cinvxy=-covxy/det
      cinvyy= covxx/det
      write(0,*) 'S',covsxx,covsxy,covsyy
      write(0,*) 'Sinv',cinvsxx,cinvsxy,cinvsyy,dets
      write(0,*) 'G',covgxx,covgxy,covgyy
      write(0,*) 'SxG',covxx,covxy,covyy
      write(0,*) 'SxGinv',cinvxx,cinvxy,cinvyy,det
      do istar=1,ns
          xc=ran3(idum)*nx
          yc=ran3(idum)*ny
          ixc=int(xc)
          iyc=int(yc)
          dxc=xc-ixc
          dyc=yc-iyc
          flx=fs*ran3(idum)
          do i=1,1000
              dx=0.1*(i-500)-dxc
              ix=i/10+ixc-50
              if (ix.lt.1 .or. ix.gt.nx) goto 2
              do j=1,1000
                  dy=0.1*(j-500)-dyc
                  e=cinvsxx*dx*dx+2*cinvsxy*dx*dy+cinvsyy*dy*dy
                  if (e.gt.100.) goto 3
                  iy=j/10+iyc-50
                  if (iy.lt.1 .or. iy.gt.ny) goto 3
                  zhi=flx*0.1**2/dets*exp(-0.5*e)
                  z(ix,iy)=z(ix,iy)+zhi
    3         enddo
    2     enddo
          write(*,*) xc,yc,flx,pa,pb,ppa,istar
      enddo
      write(0,*) 'Added ',ns,' stars.'
      do igal=1,ng
          xc=ran3(idum)*nx
          yc=ran3(idum)*ny
          ixc=int(xc)
          iyc=int(yc)
          dxc=xc-ixc
          dyc=yc-iyc
          flx=fg*(ran3(idum)+1.)/2.
          do i=1,1000
              dx=0.1*(i-500)-dxc
              ix=i/10+ixc-50
              if (ix.lt.1 .or. ix.gt.nx) goto 4
              do j=1,1000
                  dy=0.1*(j-500)-dyc
                  e=cinvxx*dx*dx+2*cinvxy*dx*dy+cinvyy*dy*dy
                  if (e.gt.100.) goto 5
                  iy=j/10+iyc-50
                  if (iy.lt.1 .or. iy.gt.ny) goto 5
                  zhi=flx*0.1**2/det*exp(-0.5*e)
                  z(ix,iy)=z(ix,iy)+zhi
    5         enddo
    4     enddo
          write(*,*) xc,yc,flx,gr,gr,ns+igal
      enddo
      write(0,*) 'Added ',ng,' galaxies.'

      do i=1,nx
          do j=1,ny
              zerr(i,j)=sqrt(z(i,j))
              z(i,j)=min(satur,z(i,j)+gasdev(idum)*zerr(i,j))
          enddo
      enddo

      call put2('fake.fits',nx,ny,z,mmpix,mmpix)
      call put2('err.fits',nx,ny,zerr,mmpix,mmpix)
      end
