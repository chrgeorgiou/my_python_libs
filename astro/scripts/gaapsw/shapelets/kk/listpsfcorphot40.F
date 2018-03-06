c     list-driven GaaP photometry 
c
c     open fits file inimage.fits, and error file inerr.fits
c     then read sets of sextractor sources,beta,q and perform shapelet
C     expansion on position x,y with shapelet scale beta
c     produce list of psfcorphot values for this shapelet at scales q

c     PSF information is gathered from psf.map file
c     PSF star list is taken from psf.sh file (for calculating PSF residuals)

c     big speedup would be to interpolate the inverse PSF coefficients,
c     rather than calculate new PSF and invert it for every source? A
c     bit less accurate though.
c
c     THIS VERSION: NO PIX RESIDUALS, JUST PSF APERTURE CORRECTION
c

      program listpsfcorphot

      parameter (maxpsf=10000)
      dimension shobs(msh),shpsf(msh),shpsfstars(msh,maxpsf)
      dimension psfcoef(mmpoly2,msh-2)
      parameter (maxstamp=800,maxstamp2=maxstamp*maxstamp)
      dimension stamp(maxstamp2),stamperr(maxstamp2)
      dimension xpsf(maxpsf),ypsf(maxpsf),sumpsf(maxpsf),bgpsf(maxpsf),
     :    apcor(maxpsf),wt(maxpsf),bestfit(maxpsf),fsex(maxpsf)
      
c first get the PSF map coefficients

      open(10,file='psf.map',status='old')
      read(10,*) npsf,betapsf
      read(10,*) nfit,xmax,ymax
      read(10,*)
      do k=1,(npsf+1)*(npsf+2)/2
         read(10,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(10,*)
      enddo
      close(10)
      write(0,*) 'Read psf model from psf.map'
      write(0,*) 'Order for shape, spatial variation; beta:',npsf,nfit,betapsf
      n=npsf
      nrejiter=3
      clip=4.

c get the list of PSF star positions from psf.sh
      open(10,file='psf.sh',status='old')
      do i=1,maxpsf
          read(10,*,end=11) xpsf(i),ypsf(i),fsex(i),(ttt,kkk=1,7)
          read(10,*) sumpsf(i),dxpsf,dypsf,bgpsf(i),(ttt,kkk=1,2)
          xpsf(i)=xpsf(i)-dxpsf
          ypsf(i)=ypsf(i)-dypsf
          read(10,*) ttt,ttt
          read(10,*) (ttt,kkk=3,(n+1)*(n+2)/2+2)
          read(10,*)
c          write(0,'(4g15.5)') xpsf(i),ypsf(i),sumpsf(i),bgpsf(i)
      enddo
   11 npsfstars=i-1
      close(10)
      write(0,*) 'Read ',npsfstars,' PSF star positions'

c open the image and error files

      call openfits('inimage.fits',51,nx,ny)
      sig=guess_errors_sec(51,nx*ny)
      do k=1,maxstamp2
         stamperr(k)=sig
      enddo
c get the PSF models at each PSF star location (for later residual calc)
c as well as the aperture flux within radius 8 beta (to normalize PSF properly)

      write(0,*) 'xPSF      yPSF     Raper       BG    '//
     :     'shplt sum   apflx   SEx flx'
      do is=1,npsfstars
          xx=xpsf(is)
          yy=ypsf(is)
c          write(0,*) '*** ',is
          call interpolpsf(xx,yy,xmax,ymax,psfcoef,npsf,betapsf,nfit,shpsfstars(1,is))
          do k=3,(npsf+1)*(npsf+2)/2+2
              shpsfstars(k,is)=shpsfstars(k,is)*sumpsf(is)
          enddo
          npix=32*betapsf+3
          ix1=max(1,min(int(xx-16*betapsf-1),nx-npix+1))
          iy1=max(1,min(int(yy-16*betapsf-1),ny-npix+1))
          xoff=xx-ix1+1
          yoff=yy-iy1+1
          isok=igetsec(51,nx,ny,ix1,iy1,npix,npix,stamp,maxstamp2)
          ap=apflux(xoff,yoff,16*betapsf,bgpsf(is),npix,npix,stamp)
          apcor(is)=ap / sumpsf(is)
          wt(is)=1
          write(0,'(3f10.2,4g14.6)') 
     :        xx,yy,16*betapsf,bgpsf(is),sumpsf(is),ap,fsex(is)
      enddo

c calculate median aperture correction factor
      apcorfac=fmedian(npsfstars,apcor)
      write(0,*) 'Aperture cor coef: ',apcorfac
cccc SET APCORFAC TO 1
      apcorfac=1

c     read through the list of sources, and process each one:
c     make a shapelet expansion
c     compute the GaaP fluxes

      betaprev=-1
      do i=1,mmcat
    1     read(*,*,err=1,end=99) x,y,f,fe,fw,a,b,th,id,iflg,beta,qq
          nread=nread+1
          fitradius=max(10.,4*beta)
          if (x.gt.fitradius+1 .and. x.lt.nx-fitradius-1
     :        .and. y.gt.fitradius+1 .and. y.lt.ny-fitradius-1 
     :        .and. beta.gt.0.) then
c background from wide annulus
            npix=16*fitradius+3
            ix1=nint(x)-(npix-1)/2
            ix1=max(1,min(ix1,nx-npix+1))
            iy1=nint(y)-(npix-1)/2
            iy1=max(1,min(iy1,ny-npix+1))
            isok=igetsec(51,nx,ny,ix1,iy1,npix,npix,stamp,maxstamp2)
            if (isok.eq.0) then
               xoff=x-ix1+1
               yoff=y-iy1+1
c check that the stamp contains non-trivial data: 
c simple check is to make sure peak is not same as any of the corners
               stamppk=stamp(npix*npix/2+1)
               if (stamppk.eq.stamp(1).or.stamppk.eq.stamp(npix)
     :              .or.stamppk.eq.stamp(npix*npix-npix+1)
     :              .or.stamppk.eq.stamp(npix*npix)) then 
                  write(*,1001) x,y,f,fe,fw,a,b,th,id,iflg,
     :                 beta,0.,0.,0.,0.,0.,0.
               else
                  call fitbg(stamp,npix,npix,npix,xoff,yoff,4*fitradius,bg)
c               write(0,*) 'bg:',bg
cccc  ****** FORCE BG TO ZERO
c     bg=1e-6
                  call fitshapeletserr(stamp,stamperr,npix,npix,npix,xoff,yoff,fitradius
     :                 ,bg,shobs,sherr,n,beta)
                  call interpolpsf(x,y,xmax,ymax,psfcoef,npsf,betapsf,nfit,shpsf)
                  fq=psfcorphot3(shobs,sherr,shpsf,qq,fqerr)
                  fq=fq*apcorfac
                  fqresid=0
                  psfresidfac=0
                  write(*,1001) x,y,f,fe,fw,a,b,th,id,iflg,
     :                 beta,qq,fq,fq+fqresid,(fq+fqresid)/(1+psfresidfac),fqerr,bg
               endif
            else
               write(*,1001) x,y,f,fe,fw,a,b,th,id,iflg,
     :              beta,0.,0.,0.,0.,0.,0.
            endif
         else
            write(*,1001) x,y,f,fe,fw,a,b,th,id,iflg,
     :           beta,0.,0.,0.,0.,0.,0.
         endif
          
 99   enddo

 1001 format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i10,i5,2f10.3,4g13.5,g14.6)
      end



          
