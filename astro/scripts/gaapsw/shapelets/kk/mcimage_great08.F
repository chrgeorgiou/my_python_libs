      parameter(mc=4000000,mmpix=4000)
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
              z(i,j)=0
          enddo
      enddo
      write(0,*) 'Initialized background'
c read all sources, add them into image
      do isource=1,100000
         if (mod(isource,20).eq.0) write(0,*) 'Doing source',isource
         iran=ran3(idum)*mc
    1     read(*,*,err=1,end=2) x,y,f,re,q,pa
          write(0,*) 'Read source ',isource,x,y,f,re,q,pa
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
c start from random point in the PSF array
              ipsf=1+mod(i+iran,mc)
              xx=xs+xp(ipsf)*fwhmpsf+x
              yy=ys+yp(ipsf)*fwhmpsf+y
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
      do j=1,npix
          do i=1,npix
              z(i,j)=z(i,j)+bg+gasdev(idum)*signoise
          enddo
      enddo
      call put2('fake.fits',npix,npix,z,mmpix,mmpix)
      write(0,*) 'Fake image fake.fits created.'
      end

          




      subroutine getgal(xg,yg,nser)
      real xg(4000000),yg(4000000)
      character s
      write(s,'(i1)') nser
      open(10,file='sersic_'//s//'.unf',status='old',form='unformatted')
      read(10)ngal
      if (ngal.ne.4000000) stop 'problem in sersic input file'
      read(10) xg,yg
      close(10)
      write(0,*) 'Read sersic_'//s//'.unf'
      return
      end


      subroutine getpsf(xp,yp,betamoffat,e1psf,e2psf)
      real xp(4000000),yp(4000000)
      character s3*3
      write(s3,'(f3.1)') real(nint(betamoffat))
      open(10,file='moffat_'//s3//'.unf',status='old',form='unformatted')
      read(10)npsf
      if (npsf.ne.4000000) stop 'problem in moffat input file'
      read(10) xp,yp
      close(10)
      do i=1,4000000
          px= xp(i)*(1-e1psf) - yp(i)*e2psf
          py=-xp(i)*e2psf     + (1+e1psf)*yp(i)
          xp(i)=px
          yp(i)=py
      enddo
      write(0,*) 'Read moffat_'//s3//'.unf'
      return
      end

      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END


      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran3(idum)-1.
        v2=2.*ran3(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END



      subroutine put2(image,nx,ny,pix,mx,my)
      real pix(mx,my)
      integer naxlen(7)
      character image*(*)
      logical simple,extend

c delete the file if it already exists

c      call deletefile(image,istat)

      nunit=49
      open(unit=nunit,file=image,status='old',err=1)
      close(nunit,status='delete')
      write(0,*) 'Deleted existing file ',image
    1 continue
      irw=0
      iblocksize=1
      istat=0
      simple=.true.
      extend=.false.
      naxis=2
      naxlen(1)=nx
      naxlen(2)=ny
c open new fits file
      call ftinit(nunit,image,iblocksize, istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Put2: Open new image failed: status=',istat
         write(0,*) image
         stop
      endif
c write required header info
      call ftphpr(nunit,simple,-32,naxis,naxlen,0,1,extend,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Put2: Write header failed: status=',istat
         write(0,*) image
         stop
      endif
c write the data.
      igroup=1
      call ftp2de(nunit,igroup,mx,nx,ny,pix,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Put2: Write data failed: status=',istat
         write(0,*) image
         stop
      endif
c close the image
      call ftclos(nunit, istat)
      end



