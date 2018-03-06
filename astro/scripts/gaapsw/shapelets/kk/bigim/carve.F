      parameter (mstar=15,mstar2=2*mstar+1,m2=32,nstarmax=10000)
      dimension starpatch(-mstar:mstar,-mstar:mstar+mstar2,nstarmax)
      dimension x(nstarmax),y(nstarmax),bg(nstarmax)

      do istar=1,nstarmax
    1     read(*,*,err=1,end=2) x(istar),y(istar)
      enddo
    2 nstar=istar-1
      write(0,*) 'Read ',nstar,' stars.'

      nim=51
      call openfits('inimage.fits',nim,nx,ny)
      write(0,*) 'Reading in star postage stamps'
      do istar=1,nstar
          ix=nint(x(istar)-mstar)
          iy=nint(y(istar)-mstar)
          xc=x(istar)-ix+1
          yc=y(istar)-iy+1
          isok=igetsec(nim,nx,ny,ix,iy,mstar2,mstar2,
     :        starpatch(-mstar,-mstar,istar),mstar2*mstar2)
          if (isok.ne.0) then
              write(0,*) isok, x(istar),y(istar)
              call ftrprt('stderr',isok)
              goto 1
          endif
      enddo
      call closefits(nim)

      nim=52
      call openfits('weights.fits',nim,nx,ny)
      write(0,*) 'Reading in weight maps for stars'
      do istar=1,nstar
          ix=nint(x(istar)-mstar)
          iy=nint(y(istar)-mstar)
          xc=x(istar)-ix+1
          yc=y(istar)-iy+1
          isok=igetsec(nim,nx,ny,ix,iy,mstar2,mstar2,
     :        starpatch(-mstar,-mstar+mstar2,istar),mstar2*mstar2)
          if (isok.ne.0) then
              write(0,*) isok, x(istar),y(istar)
              call ftrprt('stderr',isok)
              goto 1
          endif
      enddo
      call closefits(nim)

c write out star postage stamps as a single 3D file.
      call put3('stars.fits',mstar2,2*mstar2,nstar,
     :     starpatch,mstar2,2*mstar2,nstar)
      end
