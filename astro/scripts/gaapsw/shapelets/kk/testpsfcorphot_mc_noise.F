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

      dimension z(101,101),z2(101,101),psf(101,101),znonoise(101,101),zerr(101,101)
      dimension shobs(msh),shpsf(msh),shobs2(msh)
      parameter (mc=1000000)
      dimension xg(mc),yg(mc),xp(mc),yp(mc)
      parameter (pi=3.1415926535)
      character s,s3*3
      double precision photth,psfintg,psfbing,galbing,resflx

c sample sersic model, convolve with gaussian of p pixels
      write(0,*) 'n,Nser, Re, betamoffat,FWHM_PSF,dx,dy,noise level?'
      read(*,*) n,nser,re,betamoffat,fwpsf,dx,dy,gnoise
      write(s,'(i1)') nser
      open(10,file='sersic_'//s,status='old')
      do i=1,1000000
 1       read(10,*,err=1) xg(i),yg(i)
      enddo
      close(10)

      do i=1,101
          do j=1,101
              zerr(i,j)=gnoise
          enddo
      enddo

      write(s3,'(f3.1)') betamoffat
      open(10,file='moffat_'//s3,status='old')
      do i=1,mc
 2       read(10,*,err=2) xp(i),yp(i)
      enddo
      close(10)

c construct gridded PSF-convolved galaxy
      do i=1,mc
         xx=xg(i)*re+xp(i)*fwpsf/2+dx
         yy=yg(i)*re+yp(i)*fwpsf/2+dy
         if (abs(xx).lt.50. .and. abs(yy).lt.50.)
     :        znonoise(nint(51+xx),nint(51+yy))=znonoise(nint(51+xx),nint(51+yy))+1./mc
      enddo
      do i=1,mc
         xx=xp(i)*fwpsf/2+dx
         yy=yp(i)*fwpsf/2+dy
         if (abs(xx).lt.50. .and. abs(yy).lt.50.)
     :        psf(nint(51+xx),nint(51+yy))=psf(nint(51+xx),nint(51+yy))+1./mc
      enddo


c fit PSF shapelet series 
      fitradius=max(5*fwpsf,10.)
      call fitbeta(psf,101,101,101,51.+dx,51.+dy,fitradius,0.,betapsf)
      fitradius=max(4*betapsf,10.)
      call fitbeta(psf,101,101,101,51.+dx,51.+dy,fitradius,0.,betapsf)
      betapsf=betapsf*1.3
      bgin=0
      call fitshapelets(psf,101,101,101,51.+dx,51.+dy,4*betapsf,bgin,shpsf,n,betapsf)
c      call writeshapen('PSF',shpsf)
      call centershape(shpsf,shpsf,dddx,dddy)
c      write(0,*) 'Shifted by ',dddx,dddy
c      call writeshapen('PSF',shpsf)


c fit beta for convolved source
      
      fitradius=max(10*re,10.)
      call fitbeta(znonoise,101,101,101,51.+dx,51.+dy,fitradius,0.,betaobs)
      fitradius=max(4*betaobs,10.)
      call fitbeta(znonoise,101,101,101,51.+dx,51.+dy,fitradius,0.,betaobs)
      betaobs=betaobs*1.3



c     make noise realizations
      do inoise=1,100
          do i=1,101
              do j=1,101
                  z(i,j)=znonoise(i,j)+gasdev(idum)*gnoise
              enddo
          enddo
c     fit convolved source: first determine beta, scale by 1.3, fit expansion     
          nsh=(n+1)*(n+2)/2
          bgin=0
          call fitshapeletserr(z,zerr,101,101,101,51.+dx,51.+dy,4*betaobs,bgin,shobs,sherr,n,betaobs)
c          call writeshapen('OBS GAL',shobs)
          call centershape(shobs,shobs2,ddx,ddy)
c          write(0,*) 'Shifted by ',ddx,ddy
c          call writeshapen('OBS GAL',shobs2)
c          write(0,*) 'Noise realization ',inoise,': Shapelet error: ',sherr
          
          phottot=seriessumerr(shobs,sherr,phottoterr)

c now go through a range of aperture radii, and calculate the various fluxes
      do kkk=-2,2
         qq=betaobs*(1+0.2*kkk)
c         phot=psfcorphot(shobs2,shpsf,qq)
         phot2=psfcorphot3(shobs2,sherr,shpsf,qq,photerr)
c         psfint=seriessum(shpsf)
         psfintq=seriessumgauswt(shpsf,qq)
c calculate theoretical fluxes based on the true light profile of galaxy.
c first the exact result for the GaaP flux of scale qq
         photth=0
         do i=1,mc
            xx=xg(i)*re
            yy=yg(i)*re
            wt=exp(-(xx*xx+yy*yy)/4/qq**2)
            photth=photth+wt / 2.
         enddo
         photth=photth/mc
c then the binned PSF gaussian-aperture flux
c         psfintg=0
c         do i=1,mc
c             xx=xp(i)*fwpsf/2+ran3(idum)-0.5
c             yy=yp(i)*fwpsf/2+ran3(idum)-0.5
c             wt=exp(-(xx*xx+yy*yy)/2/qq**2)
c             psfintg=psfintg+wt
c         enddo
c         psfintg=psfintg/mc

c     calculate residual flux from shapelet expansion, add residual
c     without PSF cor but with gauss-weighting to the result
c
         psfintg=1
         call decodectr(shobs,z2,101,101,51+dx,51+dy)
         resflx=0
         psfbing=0
         galbing=0
         do ix=1,101
             x=ix-51-dx
             do iy=1,101
                 y=iy-51-dy
                 eq=exp(-(x*x+y*y)/2/qq**2)
                 resflx=resflx+eq*(z(ix,iy)-z2(ix,iy))
c also calculate the gaussian fluxes of the observed, binned PSF and gal
                 psfbing=psfbing+eq*psf(ix,iy)
                 galbing=galbing+eq*z(ix,iy)
             enddo
         enddo

c         write(0,*) 'Exper, theor phot00 values:',qq,phot,photth,resflx,phot/photth
         write(*,'(i5,i5,3f5.2,2f7.3,f6.1,13g13.5)') n,nser,re,betamoffat,fwpsf,betaobs,betapsf,qq,phot,resflx,photth,phot/photth,(phot+resflx)/photth,psfintq,psfintg,psfbing,psfint,galbing,phot2,photerr,phottoterr
   77    continue
      enddo
      enddo
      end
