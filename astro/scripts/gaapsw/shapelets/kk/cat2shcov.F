c read a sextractor catalogue and a PSF map;
c quantize the betas to use and decompose the sources into shapelets.
c will work much faster if input list is sorted in increasing fwhm!
c

      real pix(mmpix,mmpix),pixerr(mmpix,mmpix)
      real sh(msh),sh2(msh),shcov(msh-2,msh-2),sherr(msh-2)

      call get2err('inimage.fits','err.fits',nx,ny,pix,pixerr,mmpix,mmpix)

      open(10,file='psf.map',status='old')
      read(10,*) n,betapsf
      close(10)

      nread=0
      nfitted=0
      do i=1,mmcat
    1     read(*,*,err=1,end=99) x,y,f,fe,fw,a,b,th,id,iflg
          nread=nread+1
          beta2=fw/2.3
          if (beta2.gt.betapsf) then
              nmag=nint( 4* log(beta2/betapsf) / log(2.) )
              beta=betapsf * 2.**(0.25*nmag)
          else
              beta=betapsf
          endif
c MODIFIED 4 beta --> 5 beta (2-6-05)
          fitradius=max(10.,5*beta)
c END REPLACED
          if (x.gt.fitradius+1 .and. x.lt.nx-fitradius-1
     :        .and. y.gt.fitradius+1 .and. y.lt.ny-fitradius-1) then
              bg=0.
              call fitshapeletscov(pix,pixerr,mmpix,nx,ny,x,y,fitradius,bg,
     :            sh,sherr,shcov,n,beta)
              write(*,1001) x,y,f,fe,fw,a,b,th,id,iflg
 1001         format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i6,i5)
              call centershape(sh,sh2,dx,dy)
c              write(0,*) 'Shifted shape ctr by ',dx,dy
              sum=seriessum(sh2)
              write(*,*) sum,dx,dy,bg,fitradius,sherr(1)/sum
              write(*,*) beta,n
              write(*,*) (sh2(k)/sum,k=3,(n+1)*(n+2)/2+2)
              write(*,*)
              write(*,*) (sherr(k),k=1,(n+1)*(n+2)/2)
              write(*,*)
              write(*,*)
              nfitted=nfitted+1
          endif
      enddo
      write(*,*) '*** Cat2sh stopped reading after ',mmcat,' sources.'
   99 continue
      write(0,*) 'Cat2sh fitted shapelets to order ',n,' to ',nfitted,
     :     ' / ',nread,' sources'
      end
