c read PSF star cat from std input.
c determine a value of beta that fits all, then make shapelet expansion
c adjust centroid to remove 01 and 10 terms
c normalize to integral=1
c and write this out to stdout

c fits image that is read in must be called inimage.fits (eg via logical link)

      dimension x(mmcat),y(mmcat),f(mmcat),fe(mmcat),fw(mmcat),
     :    a(mmcat),b(mmcat),th(mmcat),id(mmcat),iflg(mmcat),beta(mmcat)
      real pix(mmpix,mmpix)
      real sh(msh),sh2(msh)

      call get2('inimage.fits',nx,ny,pix,mmpix,mmpix)
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

      do i=1,nobj
          fitradius=fw(i)/2.3*4
          call fitbg(pix,mmpix,nx,ny,x(i),y(i),fitradius*2,bg)
c WAS BUG: in line below x(i),y(i) was x,y so always same beta was used...
          call fitbeta(pix,mmpix,nx,ny,x(i),y(i),fitradius,bg,beta(i))
c          write(0,*) 'PSF STAR, FWHM, radius, bg, BETA:',
c     :        id(i),fw(i),fitradius,bg,beta(i)
      enddo
      betapsf=1.3*fmedian(nobj,beta)
      write(0,*) 'Adopted betapsf: ',betapsf
c PREVIOUS RECIPE BASED ON FWHM ONLY
c      betapsf=fmedian(nobj,fw)/2.3 * 1.2
c ensure beta is sufficiently big to have nondegenerate basis fns on pix grid
      if (betapsf.lt.0.1*n) betapsf=0.1*n
      fitradius=max(10.,min(40.,8*betapsf))
      do i=1,nobj
          if (x(i).gt.fitradius+1 .and. x(i).lt.nx-fitradius
     :        .and. y(i).gt.fitradius+1 .and. y(i).lt.ny-fitradius) then
              call fitbg(pix,mmpix,nx,ny,x(i),y(i),4*fitradius,bg)
              call fitshapelets(pix,mmpix,nx,ny,x(i),y(i),fitradius,bg,sh,
     :            n,betapsf)
              write(*,1001) x(i),y(i),f(i),fe(i),fw(i),a(i),b(i),th(i),id(i),iflg(i)
 1001         format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i6,i5)
              call centershape(sh,sh2,dx,dy)
c              write(0,*) 'Shifted shape ctr by ',dx,dy
              sum=seriessum(sh2)
              write(*,1002) sum,dx,dy,bg,fitradius,1.
 1002         format (g14.7,2f10.5,g15.7,f10.5,g15.6)
              write(*,*) betapsf,n
              write(*,*) (sh2(k)/sum,k=3,(n+1)*(n+2)/2+2)
              write(*,*)
          endif
      enddo
      write(0,'('' Psfcat2sh read'',i4,'' stars, used shapelet scale '',
     :     f6.3,'' and order '',i2)') nobj,betapsf,n
      end


