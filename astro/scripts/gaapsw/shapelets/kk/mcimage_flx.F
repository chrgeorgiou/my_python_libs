      parameter(mc=1000000)
      parameter(pi=3.1415926535)
      real z(mmpix,mmpix)
      real xp(mc),yp(mc),xg(mc),yg(mc)
      real kap
      common nonconv

      read(*,*) npix,bg,signoise,betamoffat,fwhmpsf,nsersic,e1psf,e2psf,kap,g1,g2,flexf1,flexf2,flexg1,flexg2

      write(0,*) 'Image size, bg, noise/pix:',npix,bg,signoise
      write(0,*) 'PSF beta, FWHM, e1,e2:',betamoffat,fwhmpsf,e1psf,e2psf
      write(0,*) 'Galaxy sersic param: ',nsersic
      write(0,*) 'Kappa: ',kap
      write(0,*) 'Shears 1,2: ',g1,g2
      write(0,*) 'First flexion 1,2: ',flexf1,flexf2
      write(0,*) 'Second flexion 1,2: ',flexg1,flexg2

      call lenspot(kap,g1,g2,flexf1,flexf2,flexg1,flexg2,psixx,psixy,psiyy,psixxx,psixxy,psixyy,psiyyy)

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
          nonconv=0
          do i=1,mc
              wt=1.
              if (re.gt.0.) then
c     rescale galaxy by re, apply intrinsic flattening
                  xe=xg(i)/sq*re
                  ye=yg(i)*sq*re
c     rotate glaxy to intrinsic PA
                  cpa=cos(pa*pi/180)
                  spa=sin(pa*pi/180)
                  xr=xe*cpa-ye*spa
                  yr=xe*spa+ye*cpa
c     lens galaxy photons
                  call lens(xr,yr,psixx,psixy,psiyy,psixxx,psixxy,psixyy,psiyyy,xs,ys)
                  if (xs.gt.1.e29) then
                      wt=0
                  else
                      d=detlens(xs,ys,psixx,psixy,psiyy,psixxx,psixxy,psixyy,psiyyy)
                      if (d.lt.0.1 .or. d.gt.10.) then
                          wt=0
                          nonconv=nonconv+1
                      else
                          wt=1./d
                      endif
                  endif
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
     :            z(ix,iy)=z(ix,iy)+f*wt/mc
          enddo
          if (mod(isource,500).eq.0) write(0,*) 'Done source ',isource
          write(0,*) nonconv,'/',mc,' non-converged photons!'
      enddo
    2 continue
      write(0,*) 'Read ',isource-1,' sources.'
      call put2('fake.fits',npix,npix,z,mmpix,mmpix)
      write(0,*) 'Fake image fake.fits created.'
      end

          




      subroutine getgal(xg,yg,nser)
      real xg(1000000),yg(1000000)
      character s
      write(s,'(i1)') nser
      open(10,file='sersic_'//s//'.unf',status='old',form='unformatted')
c      read(10)ngal
c      if (ngal.ne.1000000) stop 'problem in sersic input file'
      read(10) xg,yg
      close(10)
      write(0,*) 'Read sersic_'//s//'.unf'
      return
      end


      subroutine getpsf(xp,yp,betamoffat,e1psf,e2psf)
      real xp(1000000),yp(1000000)
      character s3*3
      write(s3,'(f3.1)') real(nint(betamoffat))
      open(10,file='moffat_'//s3//'.unf',status='old',form='unformatted')
c      read(10)npsf
c      if (npsf.ne.1000000) stop 'problem in moffat input file'
      read(10) xp,yp
      close(10)
c fixed the sign error in epsf's here dec 18 2007
      do i=1,1000000
          px= xp(i)*(1+e1psf) + yp(i)*e2psf
          py= xp(i)*e2psf     + (1-e1psf)*yp(i)
c make xp,yp in units of FWHM (file is in half-width HM)
          xp(i)=px/2
          yp(i)=py/2
      enddo
      write(0,*) 'Read moffat_'//s3//'.unf'
      return
      end



      subroutine lenspot(kap,g1,g2,flexf1,flexf2,flexg1,flexg2,psixx,psixy,psiyy,psixxx,psixxy,psixyy,psiyyy)
      real kap

      psixx=(kap+g1)
      psixy=g2
      psiyy=(kap-g1)
      psixxx=1.5*flexf1+0.5*flexg1
      psixxy=0.5*(flexf2+flexg2)
      psixyy=0.5*(flexf1-flexg1)
      psiyyy=1.5*flexf2-0.5*flexg2
      return
      end

      function detlens(x,y,psixx,psixy,psiyy,psixxx,psixxy,psixyy,psiyyy)
      detlens=(1-psixx-x*psixxx-y*psixxy)*(1-psiyy-x*psixyy-y*psiyyy)-(-psixy-x*psixxy-y*psixyy)**2
      return
      end


      subroutine lens(xr,yr,psixx,psixy,psiyy,psixxx,psixxy,psixyy,psiyyy,xs,ys)
      common nonconv
      xs=xr
      ys=yr
      ddx=1
      ddy=1
      dx0=0
      dy0=0
      iter=0
      do while (ddx*ddx+ddy*ddy.gt.1.e-6 .and. iter.lt.10) 
          dx=xs*psixx+ys*psixy+0.5*xs**2*psixxx+xs*ys*psixxy+0.5*ys**2*psixyy
          dy=xs*psixy+ys*psiyy+0.5*xs**2*psixxy+xs*ys*psixyy+0.5*ys**2*psiyyy
          ddx=dx-dx0
          ddy=dy-dy0
c TRY DAMPING CONVERGENCE
          xs=xr+dx
          ys=yr+dy
          dx0=dx
          dy0=dy
          iter=iter+1
      enddo
c flag non-converged rays
      if (iter.ge.10) then
          xs=1e30
          ys=1e30
          nonconv=nonconv+1
      endif
      return
      end
