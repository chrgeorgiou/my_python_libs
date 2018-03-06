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

      g1=0
      g2=0
      open(10,file='shear.in',status='old',err=3)
      read(10,*) g1,g2
      close(10)
    3 write(0,*) 'Input shear: ',g1,g2

c make blank image with background and noise
      do j=1,npix
          do i=1,npix
              z(i,j)=0
          enddo
      enddo
      write(0,*) 'Initialized background'

      covapsfxx=sapsf**2*((1+e1apsf)**2+e2apsf**2)
      covapsfxy=sapsf**2*(2*e2apsf)
      covapsfyy=sapsf**2*((1-e1apsf)**2+e2apsf**2)
      covbpsfxx=sbpsf**2*((1+e1bpsf)**2+e2bpsf**2)
      covbpsfxy=sbpsf**2*(2*e2bpsf)
      covbpsfyy=sbpsf**2*((1-e1bpsf)**2+e2bpsf**2)

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

c implement shear on galaxy cpts
          em=(1-g1**2-g2**2)**2
c do not implement magnification part of the distortion
          em=1
          sxx=sa**2*(1+ea**2+2*e1a)
          sxy=sa**2*(2*e2a)
          syy=sa**2*(1+ea**2-2*e1a)
          covaaxx=covapsfxx+((1+g1)**2*sxx+g2**2*syy+2*(1+g1)*g2*sxy)/em
          covaaxy=covapsfxy+((1+g1)*g2*sxx+(1-g1)*g2*syy+(1-g1**2+g2**2)*sxy)/em
          covaayy=covapsfyy+(g2**2*sxx+(1-g1)**2*syy+2*(1-g1)*g2*sxy)/em
          covbaxx=covbpsfxx+((1+g1)**2*sxx+g2**2*syy+2*(1+g1)*g2*sxy)/em
          covbaxy=covbpsfxy+((1+g1)*g2*sxx+(1-g1)*g2*syy+(1-g1**2+g2**2)*sxy)/em
          covbayy=covbpsfyy+(g2**2*sxx+(1-g1)**2*syy+2*(1-g1)*g2*sxy)/em
          sxx=sb**2*(1+eb**2+2*e1b)
          sxy=sb**2*(2*e2b)
          syy=sb**2*(1+eb**2-2*e1b)
          covabxx=covapsfxx+((1+g1)**2*sxx+g2**2*syy+2*(1+g1)*g2*sxy)/em
          covabxy=covapsfxy+((1+g1)*g2*sxx+(1-g1)*g2*syy+(1-g1**2+g2**2)*sxy)/em
          covabyy=covapsfyy+(g2**2*sxx+(1-g1)**2*syy+2*(1-g1)*g2*sxy)/em
          covbbxx=covbpsfxx+((1+g1)**2*sxx+g2**2*syy+2*(1+g1)*g2*sxy)/em
          covbbxy=covbpsfxy+((1+g1)*g2*sxx+(1-g1)*g2*syy+(1-g1**2+g2**2)*sxy)/em
          covbbyy=covbpsfyy+(g2**2*sxx+(1-g1)**2*syy+2*(1-g1)*g2*sxy)/em
c etc... if OK remove next lines and let above ones supercede

c          covaaxx=covapsfxx+sa**2*(1+ea**2+2*e1a)
c          covaaxy=covapsfxy+sa**2*(2*e2a)
c          covaayy=covapsfyy+sa**2*(1+ea**2-2*e1a)
c          covabxx=covapsfxx+sb**2*(1+eb**2+2*e1b)
c          covabxy=covapsfxy+sb**2*(2*e2b)
c          covabyy=covapsfyy+sb**2*(1+eb**2-2*e1b)
c          covbaxx=covbpsfxx+sa**2*(1+ea**2+2*e1a)
c          covbaxy=covbpsfxy+sa**2*(2*e2a)
c          covbayy=covbpsfyy+sa**2*(1+ea**2-2*e1a)
c          covbbxx=covbpsfxx+sb**2*(1+eb**2+2*e1b)
c          covbbxy=covbpsfxy+sb**2*(2*e2b)
c          covbbyy=covbpsfyy+sb**2*(1+eb**2-2*e1b)
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
          if (mod(isource,1000).eq.0) write(0,*) 'Done source ',isource
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
