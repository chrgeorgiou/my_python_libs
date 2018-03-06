c build fits files os multiple-gaussian images

      real z(mmpix,mmpix)

      open(10,file='mcg.par',status='old')
      read(10,*) npix,bg,signoise
      read(10,*) fapsf,sapsf,e1apsf,e2apsf
      read(10,*) fbpsf,sbpsf,e1bpsf,e2bpsf
      close(10)

      write(0,*) 'Image size, bg, noise/pix:',npix,bg,signoise
      write(0,*) 'PSF cpt A: rel.flx, sig, e1,e2',fapsf,sapsf,e1apsf,e2apsf
      write(0,*) 'PSF cpt B: rel.flx, sig, e1,e2',fbpsf,sbpsf,e1bpsf,e2bpsf
      write(0,*) 

c make blank image with background and noise
      do j=1,npix
          do i=1,npix
              z(i,j)=0
          enddo
      enddo
      write(0,*) 'Initialized background'

c read all sources, add them into image
      do isource=1,100000
    1     read(*,*,err=1,end=2) x,y,fa,fb,sa,sb,qa,qb,paa,pab
          ea=(1-qa)/(1+qa)
          eb=(1-qb)/(1+qb)
          e1a=ea*cos(3.1415926535*paa/90.)
          e2a=ea*sin(3.1415926535*paa/90.)
          e1b=eb*cos(3.1415926535*pab/90.)
          e2b=eb*sin(3.1415926535*pab/90.)
c covariance matrices for the 4 combinations of PSF and Gal cpts.
          covaaxx=sapsf**2*(1+e1apsf)+sa**2*(1+e1a)
          covaaxy=sapsf**2*( +e2apsf)+sa**2*( +e2a)
          covaayy=sapsf**2*(1-e1apsf)+sa**2*(1-e1a)
          covabxx=sapsf**2*(1+e1apsf)+sb**2*(1+e1b)
          covabxy=sapsf**2*( +e2apsf)+sb**2*( +e2b)
          covabyy=sapsf**2*(1-e1apsf)+sb**2*(1-e1b)
          covbaxx=sbpsf**2*(1+e1bpsf)+sa**2*(1+e1a)
          covbaxy=sbpsf**2*( +e2bpsf)+sa**2*( +e2a)
          covbayy=sbpsf**2*(1-e1bpsf)+sa**2*(1-e1a)
          covbbxx=sbpsf**2*(1+e1bpsf)+sb**2*(1+e1b)
          covbbxy=sbpsf**2*( +e2bpsf)+sb**2*( +e2b)
          covbbyy=sbpsf**2*(1-e1bpsf)+sb**2*(1-e1b)
c build postage stamps out to 6 sigma
          n=6*sqrt(max(sapsf,sbpsf)**2+max(sa,sb)**2)
          do ix=max(1,nint(x)-n),min(npix,nint(x)+n)
              dx=ix-x
              do iy=max(1,nint(y)-n),min(npix,nint(y)+n)
                  dy=iy-y
                  p=fapsf*fa*gauss2d(dx,dy,covaaxx,covaaxy,covaayy)
     :             +fapsf*fb*gauss2d(dx,dy,covabxx,covabxy,covabyy)
     :             +fbpsf*fa*gauss2d(dx,dy,covbaxx,covbaxy,covbayy)
     :             +fbpsf*fb*gauss2d(dx,dy,covbbxx,covbbxy,covbbyy)
                  z(ix,iy)=z(ix,iy)+p
              enddo
          enddo
          if (mod(isource,50).eq.0) write(0,*) 'Done source ',isource
      enddo

    2 continue

      idum=0
      open(10,file='ranseed.dat',status='old',err=10)
      read(10,*) idum
      close(10)
   10 continue

      do i=1,npix
          do j=1,npix
              z(i,j)=z(i,j)+bg+gasdev(idum)*signoise
          enddo
      enddo

      open(10,file='ranseed.dat',status='old',err=11)
      close(10,status='delete')
   11 open(10,file='ranseed.dat',status='new')
      write(10,*) int(ran3(idum)*1000000)
      close(10)

      write(0,*) 'Read ',isource-1,' sources.'
      call put2('fake.fits',npix,npix,z,mmpix,mmpix)
      write(0,*) 'Fake image fake.fits created.'
      end

      function gauss2d(x,y,cxx,cxy,cyy)
      det=cxx*cyy-cxy**2
      t=(cyy*x*x-2*cxy*x*y+cxx*y*y)/det
      gauss2d=exp(-0.5*t)/(2*3.1415926535*sqrt(det))
      return
      end
