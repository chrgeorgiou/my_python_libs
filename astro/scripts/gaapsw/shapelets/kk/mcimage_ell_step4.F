      parameter(mc=10000000)
      parameter(pi=3.1415926535)
      real z(mmpix,mmpix)
      real xp(mc),yp(mc),xg(mc),yg(mc)

      read(*,*) npix,bg,signoise,betamoffat,fwhmpsf,nsersic,e1psf,e2psf,g1,g2

      write(0,*) 'Image size, bg, noise/pix:',npix,bg,signoise
      write(0,*) 'PSF beta, FWHM, e1,e2:',betamoffat,fwhmpsf,e1psf,e2psf
      write(0,*) 'Galaxy sersic param, g1, g2: ',nsersic,g1,g2

c overall magnification factor to take out
      det=1/(1-g1**2-g2**2)
      det=sqrt(det)

c make PSF photons
      call getpsf(xp,yp,betamoffat,e1psf,e2psf)
c get galaxy photons. Round galaxy to start with.
      call getgal(xg,yg,nsersic)
c make blank image with background and noise
      open(17,file='ranseed.dat')
      read(17,*) idum
      close(17)
      do j=1,npix
          do i=1,npix
              z(i,j)=bg+gasdev(idum)*signoise
          enddo
      enddo
      write(0,*) 'Initialized background'
c read all sources, add them into image
      do isource=1,100000
    1     read(*,*,err=1,end=2) x,y,f,re,q,pa
          sq=sqrt(q)
          do i=1,mc
              if (re.gt.0.) then
c     rescale galaxy by re, apply intrinsic flattening
                  xe=xg(i)/sq*re
                  ye=yg(i)*sq*re
c     rotate glaxy to intrinsic PA
                  cpa=cos(pa*pi/180)
                  spa=sin(pa*pi/180)
                  xr=xe*cpa-ye*spa
                  yr=xe*spa+ye*cpa
c     shear galaxy
                  xs=(xr*(1+g1)+yr*g2)*det
                  ys=(xr*g2+(1-g1)*yr)*det
              else
                  xs=0
                  ys=0
              endif
c add random PSF displacement, centroid
              xx=xs+xp(i)*fwhmpsf+x
              yy=ys+yp(i)*fwhmpsf+y
c put down in a pixel
              ix=nint(xx)
              iy=nint(yy)
              if (ix.gt.0 .and. ix.le.npix .and. iy.gt.0 .and. iy.le.npix) 
     :            z(ix,iy)=z(ix,iy)+f/mc
          enddo
          if (mod(isource,500).eq.0) write(0,*) 'Done source ',isource
      enddo
    2 continue
      write(0,*) 'Read ',isource-1,' sources.'
      call put2('fake.fits',npix,npix,z,mmpix,mmpix)
      write(0,*) 'Fake image fake.fits created.'
      end

          




      subroutine getgal(xg,yg,nser)
      real xg(10000000),yg(10000000)
      character s
      write(s,'(i1)') nser
      open(10,file='sersic_'//s//'.unf',status='old',form='unformatted')
      read(10)ngal
      if (ngal.ne.10000000) stop 'problem in sersic input file'
      read(10) xg,yg
      close(10)
      write(0,*) 'Read sersic_'//s//'.unf'
      return
      end


      subroutine getpsf(xp,yp,betamoffat,e1psf,e2psf)
      real xp(10000000),yp(10000000)
      character s3*3
      write(s3,'(f3.1)') real(nint(betamoffat))
      open(10,file='moffat_'//s3//'.unf',status='old',form='unformatted')
      read(10)npsf
      if (npsf.ne.1000000) stop 'problem in moffat input file'
      read(10) xp,yp
      close(10)
      do i=1,10000000
          px= xp(i)*(1-e1psf) - yp(i)*e2psf
          py=-xp(i)*e2psf     + (1+e1psf)*yp(i)
c make xp,yp in units of FWHM (file is in half-width HM)
          xp(i)=px/2
          yp(i)=py/2
      enddo
      write(0,*) 'Read moffat_'//s3//'.unf'
      return
      end

