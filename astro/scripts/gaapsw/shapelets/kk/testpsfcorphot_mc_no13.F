c     set up a gaussian source, stddev sqrt(a^2+p^2), and a PSF with
c     stddev p expand in shapelets with scale bs and bp
c
c     construct PSF matrix P that convolves a source with scale q with a
c     PSF of scale bp to yield a result with scale bs
c
c     invert P, apply to source to give deconvolved source = gaussian of
c     scale a, expressed with scale bq
c
c     construct new PSF = gaussian of scale and dispersion q
c
c     construct new PSF matrix that convolves a source of scale q and a
c     PSF of scale q to a result of scale q
c
c     read off the 00 component of the result, which is the
c     gaussian-weighted aperture photometry of the source convolved with
c     gaussian PSF of scale q.

      dimension z(101,101),psf(101,101)
      dimension shobs(msh),shpsf(msh)
      parameter (mc=1000000)
      dimension xg(mc),yg(mc),xp(mc),yp(mc)
      parameter (pi=3.1415926535)
      character s,s3*3

c sample sersic model, convolve with gaussian of p pixels
      write(0,*) 'n,Nser, Re, betamoffat,FWHM_PSF?'
      read(*,*) n,nser,re,betamoffat,fwpsf
      write(s,'(i1)') nser
      open(10,file='sersic_'//s,status='old')
      do i=1,1000000
 1       read(10,*,err=1) xg(i),yg(i)
      enddo
      close(10)

      write(s3,'(f3.1)') betamoffat
      open(10,file='moffat_'//s3,status='old')
      do i=1,mc
 2       read(10,*,err=2) xp(i),yp(i)
      enddo
      close(10)

c construct gridded PSF-convolved galaxy
      do i=1,mc
         xx=xg(i)*re+xp(i)*fwpsf/2
         yy=yg(i)*re+yp(i)*fwpsf/2
         if (abs(xx).lt.50. .and. abs(yy).lt.50.)
     :        z(nint(51+xx),nint(51+yy))=z(nint(51+xx),nint(51+yy))+1./mc
      enddo
      do i=1,mc
         xx=xp(i)*fwpsf/2
         yy=yp(i)*fwpsf/2
         if (abs(xx).lt.50. .and. abs(yy).lt.50.)
     :        psf(nint(51+xx),nint(51+yy))=psf(nint(51+xx),nint(51+yy))+1./mc
      enddo
      
      nsh=(n+1)*(n+2)/2
      fitradius=max(10*re,10.)
      call fitbeta(z,101,101,101,51.,51.,fitradius,0.,betaobs)
      fitradius=max(4*betaobs,10.)
      call fitbeta(z,101,101,101,51.,51.,fitradius,0.,betaobs)
c      betaobs=betaobs*1.3
      bgin=0
      call fitshapelets(z,101,101,101,51.,51.,4*betaobs,bgin,shobs,n,betaobs)
      call writeshapen('OBS GAL',shobs)

      fitradius=max(5*fwpsf,10.)
      call fitbeta(psf,101,101,101,51.,51.,fitradius,0.,betapsf)
      fitradius=max(4*betapsf,10.)
      call fitbeta(psf,101,101,101,51.,51.,fitradius,0.,betapsf)
c      betapsf=betapsf*1.3
      bgin=0
      call fitshapelets(psf,101,101,101,51.,51.,4*betapsf,bgin,shpsf,n,betapsf)
      call writeshapen('PSF',shpsf)

      q=betapsf

      do kkk=0,10
         qq=2+kkk
         phot=psfcorphot(shobs,shpsf,qq)
         photth=0
         do i=1,mc
            xx=xg(i)*re+qq*gasdev(idum)
            yy=yg(i)*re+qq*gasdev(idum)
            wt=exp(-(xx*xx+yy*yy)/2/qq**2)
            photth=photth+wt/mc
         enddo
         write(0,*) 'Exper, theor phot00 values:',qq,phot,photth,phot/photth
         write(*,'(i5,i5,3f5.2,2f7.3,f6.1,3g13.5)') n,nser,re,betamoffat,fwpsf,betaobs,betapsf,qq,phot,photth,phot/photth
      enddo

      end
