c read a sextractor catalogue and a PSF map;
c quantize the betas to use and decompose the sources into shapelets.
c will work much faster if input list is sorted in increasing fwhm!
c

      parameter(maxstamp=800,maxstamp2=maxstamp*maxstamp)
      real stamp(maxstamp2),stamperr(maxstamp2)
      real sh(msh),sh2(msh)

      call openfits('inimage.fits',51,nx,ny)
      err=guess_errors_sec(51,nx*ny)
      do k=1,maxstamp2
         stamperr(k)=err
      enddo

      open(10,file='psf.map',status='old')
      read(10,*) n,betapsf
      close(10)

      nread=0
      nfitted=0
      do i=1,mmcat
    1     read(*,*,err=1,end=99) x,y,f,fe,fw,a,b,th,id,iflg
          nread=nread+1
c set beta from source FWHM - factor 1.2 helps with nongaussian wings
c          beta2=fw/2.3 * 1.3
c MODIFIED: set beta from best-fit gaussian, multiplied by 1.3
          fitradius=max(3.,fw/2.3*4)
          npix=16*fitradius+3
          ix1=nint(x)-(npix-1)/2
          ix1=max(1,min(ix1,nx-npix+1))
          iy1=nint(y)-(npix-1)/2
          iy1=max(1,min(iy1,ny-npix+1))
          isok=igetsec(51,nx,ny,ix1,iy1,npix,npix,stamp,maxstamp2)
          if (isok.eq.0) then
             xoff=x-ix1+1
             yoff=y-iy1+1
             call fitbg(stamp,npix,npix,npix,xoff,yoff,fitradius*2,bg)
c     write(0,*) 'Fitted bg ',bg
             call fitbeta(stamp,npix,npix,npix,xoff,yoff,fitradius,bg,beta2)
c     write(0,*) 'Fitted beta ',beta2
c     write(0,*) 'ID,FWHM,BETAFIT:',id,fw,beta2
             beta2=beta2*1.3
             if (beta2.gt.betapsf) then
c     MODIFIED: TAKE NEW BETA BIN EVERY 2^(1/8) INSTEAD OF 2^(1/4)
                nmag=nint( 8* log(beta2/betapsf) / log(2.) )
                beta=betapsf * 2.**(0.125*nmag)
c     nmag=nint( 4* log(beta2/betapsf) / log(2.) )
c     beta=betapsf * 2.**(0.25*nmag)
c     END MODIFICATIONS
             else
                beta=betapsf
             endif
c     flag un-/barely resolved sources
             if (beta-betapsf.lt.0.05*betapsf) iflg=iflg+1000
c DISABLED:    skip all flagged sources altogether
c             if (iflg.ne.0) goto 1
c             fitradius=max(10.,4*beta)
c MODIFIED fitradius to n-dependent value             
             fitradius=beta*(sqrt(1.*n)+3)
             fitradius=max(10.,fitradius)
             fitradius=min(fitradius,0.25*npix)
c try to keep sources which are within 3 beta of the edge
c             if (xoff.gt.fitradius+1 .and. xoff.lt.npix-fitradius-1
c     :            .and. yoff.gt.fitradius+1 .and. yoff.lt.npix-fitradius-1) then
             if (xoff.gt.3*beta+1 .and. xoff.lt.npix-3*beta
     :            .and. yoff.gt.3*beta+1 .and. yoff.lt.npix-3*beta) then
                bg=0.
                call fitbg(stamp,npix,npix,npix,xoff,yoff,2*fitradius,bg)
                call fitshapeletserr(stamp,stamperr,npix,npix,npix,xoff,yoff,fitradius,bg,
     :           sh,sherr,n,beta)
c                call plotfittedshapelet(stamp,npix,npix,npix,xoff,yoff,fitradius,bg,sh)
                call centershape(sh,sh2,dx,dy)
c FLAG SOURCES WHERE CTR SHIFTED BY MORE THAN SCALE RADIUS, REJECT SHIFT
                if (dx*dx+dy*dy.gt.beta*beta) then
                   iflg=iflg+2000
                   dx=0
                   dy=0
                   do k=3,(n+1)*(n+2)/2+2
                      sh2(k)=sh(k)
                   enddo
                endif
                sum=seriessum(sh2)
                write(*,1001) x,y,f,fe,fw,a,b,th,id,iflg
 1001           format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i10,i5)
                write(*,1002) sum,dx,dy,bg,fitradius,sherr/sum
 1002           format (g14.7,2f10.5,g15.7,f10.5,g15.6)
                write(*,*) beta,n
                write(*,*) (sh2(k)/sum,k=3,(n+1)*(n+2)/2+2)
                write(*,*)
                nfitted=nfitted+1
             endif
          endif
       enddo
      write(*,*) '*** Cat2sh stopped reading after ',mmcat,' sources.'
   99 continue
      write(0,*) 'Cat2sh fitted shapelets to order ',n,' to ',nfitted,
     :     ' / ',nread,' sources'
      end
