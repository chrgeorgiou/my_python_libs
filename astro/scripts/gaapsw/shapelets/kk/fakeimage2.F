      real z(mmpix,mmpix),zerr(mmpix,mmpix)
      parameter (satur=50000.)

c makes realization of PSF stars and elliptical, randomly oriented galaxies.
c all are simple 2d gaussians

c read parameters for PSF stars, and shape of the galaxy
      read(*,*) idum,nx,ny,fs,pa,pb,ppa,ns,fg,ga,gb,g1,g2,ng,bg
      write(0,'('' Image size to make, bg: '',2i6,g9.2)') nx,ny,bg
      write(0,'('' Stars: flux, rms maj, rms min, PA, no.:'',
     :    g9.2,2f6.2,f6.0,i6)') fs,pa,pb,ppa,ns
      write(0,'('' Gals:  flux, rms x, rms y, g1, g2, no.:'',
     :    g9.2,4f6.2,i6)') fg,ga,gb,g1,g2,ng
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
      write(0,*) 'S',covsxx,covsxy,covsyy
      write(0,*) 'Sinv',cinvsxx,cinvsxy,cinvsyy,dets
      do istar=1,ns
          xc=ran3(idum)*nx
          yc=ran3(idum)*ny
          ixc=int(xc)
          iyc=int(yc)
          dxc=xc-ixc
          dyc=yc-iyc
          flx=fs*(0.01+0.99*ran3(idum))
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
          gpa=ran3(idum)*3.1415926535*2
          cgpa=cos(gpa)
          sgpa=sin(gpa)
          gxx=ga*ga*cgpa**2 + gb*gb*sgpa**2
          gxy=(ga*ga-gb*gb)*sgpa*cgpa
          gyy=gb*gb*cgpa**2 + ga*ga*sgpa**2
          covgxx=(1+g1)**2 * gxx + 2*(1+g1)*g2 * gxy + g2**2 * gyy
          covgxy=(1+g1)*g2 * gxx + (1-g1**2+g2**2) * gxy + (1-g1)*g2 * gyy
          covgyy=g2**2 * gxx + 2*(1-g1)*g2 * gxy + (1-g1)**2 * gyy
          covxx=covsxx+covgxx
          covxy=covsxy+covgxy
          covyy=covsyy+covgyy
          det=covxx*covyy-covxy**2
          cinvxx= covyy/det
          cinvxy=-covxy/det
          cinvyy= covxx/det
c          write(0,*) 'G',covgxx,covgxy,covgyy
c          write(0,*) 'SxG',covxx,covxy,covyy
c          write(0,*) 'SxGinv',cinvxx,cinvxy,cinvyy,det
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
          write(*,*) xc,yc,flx,ga,gb,gpa*180/3.1415926535,ns+igal
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
