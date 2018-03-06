      dimension z1d(mmpix*mmpix)

      nunit1=49
      call openfits('inimage.fits',nunit1)

      read(*,*) ix1,iy1,nx,ny
      if (igetsec(nunit1,ix1,iy1,nx,ny,z1d,mmpix*mmpix).eq.0) then
         call dive(z1d,nx,ny)
      endif

      call closefits(nunit1)
      end

      subroutine dive(z,nx,ny)
      dimension z(nx,ny)
      write(*,*) z(1,ny),z(nx,ny)
      write(*,*) z(1,1),z(nx,1)
      return
      end
