c     appends fitted beta value * 1.3 to a 10-column sextractor output list

      program listfitbeta

      parameter (maxstamp=800,maxstamp2=maxstamp*maxstamp)
      real stamp(maxstamp2)
      
c first get the PSF map coefficients
      open(10,file='psf.map',status='old')
      read(10,*) npsf,betapsf
      close(10)
      write(0,*) 'Read psf model from psf.map'
      write(0,*) 'beta for psf:',betapsf

c open the fits file


      call openfits('inimage.fits',51,nx,ny)
c     read through the list of sources, and process each one:
c     fit beta
c     make a shapelet expansion
c     calculate the PSF
c     compute the GaaP fluxes

      do i=1,mmcat
    1     read(*,*,err=1,end=99) x,y,f,fe,fw,a,b,th,id,iflg
          nread=nread+1
          fitradius=max(10.,fw/2.3*4)
          if (x.gt.fitradius+1 .and. x.lt.nx-fitradius-1
     :        .and. y.gt.fitradius+1 .and. y.lt.ny-fitradius-1) then
              npix=8*fitradius+3
              ix1=nint(x)-(npix-1)/2
              ix1=max(1,min(ix1,nx-npix+1))
              iy1=nint(y)-(npix-1)/2
              iy1=max(1,min(iy1,ny-npix+1))
              xoff=x-ix1+1
              yoff=y-iy1+1
              isok=igetsec(51,nx,ny,ix1,iy1,npix,npix,stamp,maxstamp2)
              if (isok.eq.0) then
                  call fitbg(stamp,npix,npix,npix,xoff,yoff,fitradius*2,bg)
c     write(0,*) '1',fitradius,bg
                  call fitbeta(stamp,npix,npix,npix,xoff,yoff,fitradius,bg,beta2)
c     quantize beta in powers of 2**1/8 times PSF beta
                  beta2=beta2*1.3
                  if (beta2.gt.betapsf) then
                      nmag=nint( 8* log(beta2/betapsf) / log(2.) )
                      beta=betapsf * 2.**(0.125*nmag)
                  else
                      beta=betapsf
                  endif
              else
                  beta=0
              endif
          else
              beta=0
          endif
c     write out the result
          write(*,1001) x,y,f,fe,fw,a,b,th,id,iflg,beta
 1001     format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i10,i5,f10.3)
          
   99 enddo

      end
