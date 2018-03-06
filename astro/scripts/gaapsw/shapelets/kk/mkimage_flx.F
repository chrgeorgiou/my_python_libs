      parameter(pi=3.1415926535)
      real z(mmpix,mmpix)
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
          do j=1,npix
              dy=j-y
              do i=1,npix
                  dx=i-x
                  xl=dx-dx*psixx-dy*psixy-0.5*dx**2*psixxx-dx*dy*psixxy-0.5*dy*dy*psixyy
                  yl=dy-dx*psixy-dy*psiyy-0.5*dx**2*psixxy-dx*dy*psixyy-0.5*dy*dy*psiyyy
                  xl=xl/sq
                  yl=yl*sq
                  cpa=cos(pa*pi/180)
                  spa=sin(pa*pi/180)
                  xn= xl*cpa+yl*spa
                  yn=-xl*spa+yl*cpa
                  z(i,j)=z(i,j)+f*gmod(xn,yn,re,nsersic)
              enddo
          enddo
          if (mod(isource,500).eq.0) write(0,*) 'Done source ',isource
      enddo
    2 continue
      write(0,*) 'Read ',isource-1,' sources.'
      call put2('fake.fits',npix,npix,z,mmpix,mmpix)
      write(0,*) 'Fake image fake.fits created.'
      end

          
c DOES NOT SEEM TO AGREE WITH MC SIMULATIONS! CHECK!!
      function gmod(x,y,re,nser)
      gmod=0
      r=sqrt(x*x+y*y)
      if (nser.eq.1) gmod=exp(-1.67*r/re)/(6.28*re**2/1.67**2)
      if (nser.eq.4) gmod=exp(-7.67*((r/re)**0.25-1))/(7.22*3.14*re**2)
      return
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
c CHECK THESE CONVENTIONS?
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
          xs=xr+dx/2
          ys=yr+dy/2
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
