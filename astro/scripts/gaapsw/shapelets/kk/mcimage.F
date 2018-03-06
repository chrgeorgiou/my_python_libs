      parameter(mc=1000000)
      real z(mmpix,mmpix)
      real xp(mc),yp(mc),xg(mc),yg(mc)

      read(*,*) npix,bg,signoise,betamoffat,fwhmpsf,nsersic
      write(0,*) 'npix,bg,signoise,betamoffat,fwhmpsf,nsersic: ',
     :    npix,bg,signoise,betamoffat,fwhmpsf,nsersic
      call getpsf(xp,yp,betamoffat)
      call getgal(xg,yg,nsersic)
      do j=1,npix
          do i=1,npix
              z(i,j)=bg+gasdev(idum)*signoise
          enddo
      enddo
      write(0,*) 'Initialized background'
      do isource=1,100000
    1     read(*,*,err=1,end=2) x,y,f,re
          do i=1,mc
              xx=xg(i)*re+xp(i)*fwhmpsf+x
              yy=yg(i)*re+yp(i)*fwhmpsf+y
              ix=nint(xx)
              iy=nint(yy)
c              if (isource.eq.2) write(0,*) i,xg(i),yg(i),re,xp(i),yp(i),fwhmpsf,xx,yy,ix,iy
              if (ix.gt.0 .and. ix.le.npix .and. iy.gt.0 .and. iy.le.npix) 
     :            z(ix,iy)=z(ix,iy)+f/mc
          enddo
          if (mod(isource,500).eq.0) write(0,*) 'Done source ',isource
      enddo
    2 continue
      call put2('fake.fits',npix,npix,z,mmpix,mmpix)
      write(0,*) 'Fake image fake.fits created.'
      end

          




      subroutine getgal(xg,yg,nser)
      real xg(1000000),yg(1000000)
      character s
      write(s,'(i1)') nser
      open(10,file='sersic_'//s//'.unf',status='old',form='unformatted')
      read(10) xg,yg
      close(10)
      write(0,*) 'Read sersic_'//s//'.unf'
      return
      end


      subroutine getpsf(xp,yp,betamoffat)
      real xp(1000000),yp(1000000)
      character s3*3
      ellpsf=betamoffat-nint(betamoffat)
      e1psf=ellpsf*cos(1.)
      e2psf=ellpsf*sin(1.)
      write(s3,'(f3.1)') real(nint(betamoffat))
      open(10,file='moffat_'//s3//'.unf',status='old',form='unformatted')
      read(10) xp,yp
      close(10)
      do i=1,1000000
          px= xp(i)*(1-e1psf) - yp(i)*e2psf
          py=-xp(i)*e2psf     + (1+e1psf)*yp(i)
c make xp,yp in units of FWHM (file is in half-width HM)
          xp(i)=px/2
          yp(i)=py/2
      enddo
      write(0,*) 'Read moffat_'//s3//'.unf'
      return
      end

