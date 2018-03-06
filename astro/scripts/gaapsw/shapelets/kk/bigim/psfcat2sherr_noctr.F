c read PSF star cat from std input.
c determine a value of beta that fits all, then make shapelet expansion
c adjust centroid to remove 01 and 10 terms
c normalize to integral=1
c and write this out to stdout

c fits image that is read in must be called inimage.fits (eg via logical link)
c here read in postage stamps around each object


      dimension x(mmcat),y(mmcat),f(mmcat),fe(mmcat),fw(mmcat),
     :    a(mmcat),b(mmcat),th(mmcat),id(mmcat),iflg(mmcat),beta(mmcat)
      real sh(msh)
      parameter (maxstamp=200,maxstamp2=maxstamp*maxstamp)
      real stamp(maxstamp2)

c open image
      call openfits('inimage.fits',51,nx,ny)
c get error on pixels from statistics
      sherr=guess_errors_sec(51,nx*ny)

c read all PSF stars
      call readorders(n,nfit)
      do i=1,mmcat
 1       read(*,*,err=1,end=99) 
     :        x(i),y(i),f(i),fe(i),fw(i),a(i),b(i),th(i),id(i),iflg(i)
      enddo
      write(0,*) 'There may be more objects than the ',mmcat,' read.'
 99   nobj=i-1
      write(0,*) 'Read ',nobj,' stars.'

c assign single PSF beta and fitting radius
c MODIFIED:
c find best-fit beta to all PSF objects, take median, multiply by 1.3
c the factor 1.3 is result of simulations that show 
c this helps with moffat profiles, and does no harm for gaussian PSFs.

c FIXED BUG: fitbeta call was with arguments x,y instead of x(i),y(i)

c record background and beta fits in a file. Delete previous version is exists.
      open(23,file='starpars.dat',status='old',err=2)
      close(23,status='delete')
    2 continue
      open(23,file='starpars.dat',status='new')

      do i=1,nobj
          fitradius=max(10.,fw(i)/2.3*4)
          npix=8*fitradius+3
c          ix1=nint(x(i))-(npix-1)/2
c          iy1=nint(y(i))-(npix-1)/2
          ix1=max(1,min(nx-npix+1,nint(x(i))-(npix-1)/2))
          iy1=max(1,min(ny-npix+1,nint(y(i))-(npix-1)/2))
          isok=igetsec(51,nx,ny,ix1,iy1,npix,npix,stamp,maxstamp2)
          if (isok.eq.0) then
             call fitbg(stamp,npix,npix,npix,x(i)-ix1+1,y(i)-iy1+1,fitradius*2,bg)
             call fitbeta(stamp,npix,npix,npix,x(i)-ix1+1,y(i)-iy1+1,fitradius,bg,beta(i))
             fitradius=max(10.,1.3*beta(i)*4)
             fitradius=min(fitradius,0.25*npix)
             call fitbg(stamp,npix,npix,npix,x(i)-ix1+1,y(i)-iy1+1,fitradius*2,bg)
             call fitbeta(stamp,npix,npix,npix,x(i)-ix1+1,y(i)-iy1+1,fitradius,bg,beta(i))
          write(23,*) 'PSF STAR, FWHM, radius, bg, BETA:',
     :        id(i),fw(i),fitradius,bg,beta(i)
          endif
      enddo
      close(23)

      betapsf=1.3*fmedian(nobj,beta)
      write(0,*) 'Adopted betapsf: ',betapsf
c PREVIOUS RECIPE BASED ON FWHM ONLY
c      betapsf=fmedian(nobj,fw)/2.3 * 1.2
c ensure beta is sufficiently big to have nondegenerate basis fns on pix grid
      if (betapsf.lt.0.1*n) betapsf=0.1*n
c replace
c      fitradius=max(10.,min(40.,8*betapsf))
c by
      fitradius=max(10.,4*betapsf)
      do i=1,nobj
          if (x(i).gt.fitradius+1 .and. x(i).lt.nx-fitradius
     :        .and. y(i).gt.fitradius+1 .and. y(i).lt.ny-fitradius) then

             npix=8*fitradius+3
             ix1=max(1,min(nx-npix+1,nint(x(i))-(npix-1)/2))
             iy1=max(1,min(ny-npix+1,nint(y(i))-(npix-1)/2))
c             write(0,*) 'Getting section ',ix1,iy1,npix,npix,'...'
             isok=igetsec(51,nx,ny,ix1,iy1,npix,npix,stamp,maxstamp2)
c             write(0,*) '...Done.'
             call fitbg(stamp,npix,npix,npix,x(i)-ix1+1,y(i)-iy1+1,4*fitradius,bg)

              call fitshapelets(stamp,npix,npix,npix,x(i)-ix1+1,y(i)-iy1+1,2*fitradius,bg,sh,
     :            n,betapsf)
c              call plotfittedshapelet(stamp,npix,npix,npix,x(i)-ix1+1,y(i)-iy1+1,2*fitradius,bg,sh)
c              call centershape(sh,sh2,dx,dy)
c              write(0,*) 'Shifted shape ctr by ',dx,dy
              sum=seriessumerr(sh,sherr,sumerr)
c only write out S/N > 20 sources
              dx=0
              dy=0
              if (sum.ge.sumerr*20) then
                  write(*,1001) x(i),y(i),f(i),fe(i),fw(i),a(i),b(i),th(i),id(i),iflg(i)
 1001             format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i10,i5)
                  write(*,*) sum,dx,dy,bg,fitradius,1.
                  write(*,*) betapsf,n
                  write(*,*) (sh(k)/sum,k=3,(n+1)*(n+2)/2+2)
                  write(*,*)
              endif
          endif
      enddo
      write(0,'('' Psfcat2sh read'',i7,'' stars, used shapelet scale '',
     :     f6.3,'' and order '',i2)') nobj,betapsf,n
      end


