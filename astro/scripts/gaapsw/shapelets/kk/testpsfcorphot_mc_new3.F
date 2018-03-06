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

      dimension z(101,101),z2(101,101),psf(101,101),psf2(101,101)
      dimension shobs(msh),shpsf(msh)
      parameter (mc=1000000)
      dimension xg(mc),yg(mc),xp(mc),yp(mc)
      parameter (pi=3.1415926535)
      double precision photth,resflx,respsf

      idum=-1
      open(10,file='dummy.seed',status='old',err=222)
      read(10,*) idum
      close(10)
  222 continue

      write(0,*) 'n,Nser, Re, betamoffat,FWHM_PSF,dx,dy?'
      read(*,*) n,nser,re,betamoffat,fwpsf,dx,dy

      call getgal(xg,yg,nser)
      call getpsf(xp,yp,betamoffat)

c construct gridded PSF-convolved galaxy
      do i=1,mc
         xx=xg(i)*re+xp(i)*fwpsf/2+dx
         yy=yg(i)*re+yp(i)*fwpsf/2+dy
         if (abs(xx).lt.50. .and. abs(yy).lt.50.)
     :        z(nint(51+xx),nint(51+yy))=z(nint(51+xx),nint(51+yy))+1./mc
      enddo
c construct gridded PSF
      do i=1,mc
         xx=xp(i)*fwpsf/2+dx
         yy=yp(i)*fwpsf/2+dy
         if (abs(xx).lt.50. .and. abs(yy).lt.50.)
     :        psf(nint(51+xx),nint(51+yy))=psf(nint(51+xx),nint(51+yy))+1./mc
      enddo
c fit convolved source: first determine beta, scale by 1.3, fit expansion     
      nsh=(n+1)*(n+2)/2
      fitradius=max(10*re,10.)
      call fitbeta(z,101,101,101,51.+dx,51.+dy,fitradius,0.,betaobs)
      fitradius=max(4*betaobs,10.)
      call fitbeta(z,101,101,101,51.+dx,51.+dy,fitradius,0.,betaobs)
      betaobs=betaobs*1.3
      call fitbg(z,101,101,101,51.+dx,51.+dy,4*betaobs,bgobs)
cccccccc fix b/g to correct value
c      bgobs=1e-10
cccccccc
      call fitshapelets(z,101,101,101,51.+dx,51.+dy,4*betaobs,bgobs,shobs,n,betaobs)
      write(0,*) 'B/G ',bgobs
      call writeshape('OBS GAL',shobs)

c fit PSF shapelet series in same way
      fitradius=max(5*fwpsf,10.)
      call fitbeta(psf,101,101,101,51.+dx,51.+dy,fitradius,0.,betapsf)
      fitradius=max(4*betapsf,10.)
      call fitbeta(psf,101,101,101,51.+dx,51.+dy,fitradius,0.,betapsf)
      betapsf=betapsf*1.3
      call fitbg(psf,101,101,101,51.+dx,51.+dy,4*betapsf,bgpsf)
cccccccc fix b/g to correct value
c      bgpsf=1e-10
cccccccc
      call fitshapelets(psf,101,101,101,51.+dx,51.+dy,4*betapsf,bgpsf,shpsf,n,betapsf)
      write(0,*) 'B/G ',bgpsf
      call writeshape('PSF',shpsf)

c now go through a range of aperture radii, and calculate the various fluxes
c     here calculate a compensated aperture flux (Fq - F_sqrt(2.5)q /2.5)
C     , ditto for the theoretical one. The compensated flux is
C     equivalent to convolving the raw galaxy with a gaussian of
C     dispersion q, and then calculating the weighted flux with wt fn G(q)-G(2q)/4
      do kkk=0,31
         qq=betaobs*(0.3+2.2*ran3(idum))
         phot2=psfcorphot3(shobs,1.,shpsf,qq,photerr)
     :       - 0.4*psfcorphot3(shobs,1.,shpsf,qq*sqrt(2.5),photerr2)
         qpres2=max(2*qq**2-(betapsf/1.3)**2,0.5)
         qcor2=max(0.,(betaobs/1.3)**2-(betapsf/1.3)**2) + qpres2
         qcor=sqrt(qcor2)
         psfintq=seriessumgauswt(shpsf,qcor)
c calculate theoretical fluxes based on the true light profile of galaxy.
c first the exact result for the GaaP flux of scale qq
         photth=0
         do i=1,mc
            xx=xg(i)*re
            yy=yg(i)*re
            wt=exp(-(xx*xx+yy*yy)/4/qq**2)
            wt=0.5*wt - 0.2*wt**0.4
            photth=photth+wt
         enddo
         photth=photth/mc
c calculate gauss-weighted residuals from source and PSF
ccccccccc  tests of residual pix flux:
ccccccccc  NO residuals
         resflx=0
         respsf=0
ccccccccc using fitted b/g 
c         resflx=pixresid(shobs,51+dx,51+dy,z,101,101,101,bgobs,qpres2)
c         resflx=resflx*qq**2/qpres2
c         respsf=pixresid(shpsf,51+dx,51+dx,psf,101,101,101,bgpsf,qcor2)
c         respsf=respsf*(2*qq**2+max(0.,(betaobs/1.3)**2-(betapsf/1.3)**2))/qcor2
ccccccccc  using TRUE background level (0)
c         resflx=pixresid(shobs,51+dx,51+dy,z,101,101,101,0.,qpres2)
c         resflx=resflx*qq**2/qpres2
c         respsf=pixresid(shpsf,51+dx,51+dx,psf,101,101,101,0.,qcor2)
c         respsf=respsf*(2*qq**2+max(0.,(betaobs/1.3)**2-(betapsf/1.3)**2))/qcor2
ccccccccc         
         write(0,'(8f12.5)') qq,photth,phot2/photth,(phot2+resflx)/photth,(phot2+resflx)/(1+respsf)/photth
         write(*,'(i5,i5,3f5.2,2f7.3,f8.3,13g13.5)') n,nser,re,betamoffat,fwpsf,betaobs,betapsf,qq,photth,phot2,photerr,resflx,respsf,(phot2+resflx)/(1+respsf)
   77    continue
      enddo

      write(19,'(2i4,9f15.7)') n,nser,re,betamoffat,fwpsf,dx,dy,bgpsf,betapsf,bgobs,betaobs

      open(10,file='dummy.seed',status='unknown')
      write(10,*) -int(ran3(idum)*1.e9)
      close(10)
      end



      subroutine getgal(xg,yg,nser)
      real xg(*),yg(*)
      character s
      write(s,'(i1)') nser
      open(10,file='sersic_'//s,status='old')
      do i=1,1000000
 1       read(10,*,err=1) xg(i),yg(i)
      enddo
      close(10)
      return
      end


      subroutine getpsf(xp,yp,betamoffat)
      real xp(*),yp(*)
      character s3*3
      ellpsf=betamoffat-nint(betamoffat)
      e1psf=ellpsf*cos(1.)
      e2psf=ellpsf*sin(1.)
      write(s3,'(f3.1)') real(nint(betamoffat))
      open(10,file='moffat_'//s3,status='old')
      do i=1,1000000
 2       read(10,*,err=2) xp(i),yp(i)
          px= xp(i)*(1-e1psf) - yp(i)*e2psf
          py=-xp(i)*e2psf     + (1+e1psf)*yp(i)
          xp(i)=px
          yp(i)=py
      enddo
      close(10)
      return
      end

