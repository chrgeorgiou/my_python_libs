c get FWHM from PSF shapelet map
c
c input: 
c name of PSF map file
c x,y coordinates (pixels)
c
c output:
c table of x,y,FWHM
c

      character*132 mapfile
      dimension psfcoef(mmpoly2,msh-2),sh(msh),shround(mm)
      logical readpsfmap

c get name of map file, read it
      read (*,'(a132)') mapfile
      if (.not. readpsfmap(mapfile,psfcoef,nsh,beta,nfit,xmax,ymax)) then
         write(0,*) 'Error reading map ',mapfile
         stop
      endif

c read x,y from stdin, compute FWHM, write out
 1    read(*,*,err=1,end=2) x,y

c     get shapelet model at image position x,y
      call interpolpsf(x,y,xmax,ymax,psfcoef,nsh,beta,nfit,sh)

c     calculate the (n,0) coefficients of the circularised part of this shapelet. 
      call circularizeshapelet(sh,shround)

c     get fwhn
      fwhm = shfwhm(shround,nsh,beta,1e-2*beta)
      write(*,*) x,y,fwhm
      goto 1
 2    end

c-----------------------------------------------------
c     find the FWHM for a circularized shapelet by inverse interpolation. Accuracy rtol.

      function shfwhm(shround,nsh,beta,rtol)
      dimension shround(nsh/2+1)

      r0=0.
      s0=roundshapeletpix(r0,shround,nsh,beta)
      if (s0.le.0.) then
         shfwhm=-1
         return
      endif
      starget=0.5*s0
      r1=beta
      s1=roundshapeletpix(r1,shround,nsh,beta)
c bracket half-max radius
      do while (s1.gt.starget)
         r0=r1
         s0=s1
         r1=r1+beta
         s1=roundshapeletpix(r1,shround,nsh,beta)
      enddo
c bisection search
      iter=0
      do while (r1-r0.gt.rtol .and. iter.lt.20)
         iter=iter+1
         rmid=0.5*(r0+r1)
         smid=roundshapeletpix(rmid,shround,nsh,beta)
         if (smid.gt.starget) then
            r0=rmid
            s0=smid
         else
            r1=rmid
            s1=smid
         endif
c         write(0,*) 'bracket ',r0,r1,s0,s1
      enddo
c return interpolated solution in bracket. turn radius into full width
      shfwhm=2*(r0+(s0-starget)/(s0-s1)*(r1-r0))
      return
      end

c---------------------------------

c read in a PSF map file. Coefficients go into array psf. First index of
c psf() gives spatial variation; second index is shapelet
c coefficient. At position x,y in the image the shapelet model for the
c PSF is
c
c sum_{ik} PSF(i,k) SPACE(i) SHAPE(k)
c
c SPACE(i) are terms of the form X^aY^b in the order
c (ab=00, 10,01, 20,11,02, 30,21,12,03, 40,...) up to a+b=nfit, 
c where X =(2*x/xmax-1), Y=(2*y/ymax-1)
c
c SHAPE(k) are normalized shapelets S_{ab}(dx,dy,beta) in the same order
c as above (ab=00,10,01,20,11,02,30,21,12,03,40,...) up to a+b=nsh.
c
c Format of psfmap file:
c
c  nsh, beta   ! shapelet order, scale radius
c  nfit,xmax,ymax  ! spatial variation polynomial order, scaling for x and y
c  <blank line>
c and then for each shapelet coefficient k, i.e. (nsh+1)*(nsh+2)/2 times:
c  - a line with the (nfit+1)*(nfit+2)/2 polynomial components psf(i,k)  
c  - followed by a blank line

      logical function readpsfmap(filename,psf,nsh,beta,nfit,xmax,ymax)
c read in a shapelet psfmap called 'filename'.
c argument is TRUE on success, FALSE on failure
      character *(*) filename
      real psf(mmpoly2,msh-2)
      readpsfmap=.false.
      open(23,file=filename,status='old',err=2)
 1    read(23,*,err=1) nsh,beta
      read(23,*) nfit,xmax,ymax
      read(23,*) 
      do k=1,(nsh+1)*(nsh+2)/2
         read(23,*,err=2) (psf(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(23,*)
      enddo
      readpsfmap=.true.
 2    return
      end

c---------------------------------------------------------------
c interpolate psf map coefs coef() to spatial position x,y
c use scale factors xm,ym for field position
c n,beta are shapelet order and scale radius
c nfit is spatial polynomial order of coef()
c sh is the resulting shapelet array

      subroutine interpolpsf(x,y,xm,ym,coef,n,beta,nfit,sh)
      real coef(mmpoly2,msh-2),sh(msh)
      sh(1)=beta
      sh(2)=n
      do k=1,(n+1)*(n+2)/2
          sum=0
          kfit=0
          do i=0,nfit
              do m=0,i
                  l=i-m
                  kfit=kfit+1
                  sum=sum+coef(kfit,k)*polymodel(x,y,xm,ym,l,m)
              enddo
          enddo
          sh(k+2)=sum
      enddo
      return
      end

      function polymodel(x,y,xm,ym,l,m)
c simple power of x and y, with range 0..xm remapped to -1..1
      if (l.eq.0) then
          px=1
      else
          px=(2*x/xm-1)**l
      endif
      if (m.eq.0) then
          py=1
      else
          py=(2*y/ym-1)**m
      endif
      polymodel=px * py
      return
      end

c---------------------------------------------------
c     circularize shapelet sh. Result is the coeffs of Sn0, n=0,2,4,...
c
      subroutine circularizeshapelet(sh,shround)
      real sh(msh),shround(mm+1)

      beta=sh(1)
      nsh=sh(2)
      do nn=0,nsh,2
         m=nn/2
         k0=3+nn*(nn+1)/2
         cc=1
         sum2=cc**2
         sumcoef=cc*sh(k0)
         do i=1,m
            i2=2*i
            mi2=2*(m-i)
            cc=cc*(m-i+1)/i*sqrt(i2*(i2-1)/(mi2+2.)/(mi2+1.))
            sum2=sum2+cc**2
            sumcoef=sumcoef+cc*sh(k0+i2)
         enddo
         shround(m+1)=sumcoef/sum2
      enddo
      return
      end

c---------------------------------------------------------
c evaluate shapelet model sh at position x,y with respect to the center
c of the parent gaussian

      function shapeletpix(x,y,sh)
      real sh(msh),xtab(mm+1),ytab(mm+1)
      beta=sh(1)
      n=sh(2)
      k=2
      xnorm=x/beta
      call htab(xnorm,n,xtab)
      ex=exp(-0.5*xnorm**2)
      ynorm=y/beta
      call htab(ynorm,n,ytab)
      ey=exp(-0.5*ynorm**2)
      zz=0
      fac=ex*ey/sqrt(3.1415926535)/beta
      do nn=0,n
         do j=0,nn
            i=nn-j
            k=k+1
            if (j.ne.0) then
               fac=fac/sqrt(j/real(i+1))
            else 
               if (i.ne.0) fac=fac/sqrt(real(i))
            endif
            zz=zz+fac*sh(k)*xtab(i+1)*ytab(j+1)
         enddo
         fac=fac/sqrt(2.)
      enddo
      shapeletpix=zz
      return
      end

c as shapeletpix but now for radially symmetric shapelet.
c shround is the shapelet coeff of (n,0) for even n. The coeffs of Hn-r(x)Hr(y) are then simply 
c (n-r,r) = n/2Cr/2 for even n,r.
      function roundshapeletpix(r,shround,n,beta)
      real shround(mm),rtab(mm+1),ytab(mm+1)
      rnorm=r/beta
      call htab(rnorm,n,rtab)
      call htab(0.,n,ytab)
      e=exp(-0.5*rnorm**2)
      zz=0
c shapelet Sab normalized by 1/sqrt(2^(a+b) a!b! pi). fac is the coeff for the a=n,b=0 case.
      fac=e/sqrt(3.1415926535)/beta
      do nn=0,n,2
         m=nn/2
         k0=3+nn*(nn+1)/2
         cc=1
         zz=zz+cc*fac*shround(m+1)*rtab(nn+1)*ytab(1)
         do i=1,m
            cc=cc*(m-i+1)/i
            i2=2*i
            zz=zz+cc*fac*shround(m+1)*rtab(nn+1-i2)*ytab(i2+1)
c            write(0,*) 'rnd',nn,i2,cc,cc*fac*shround(m+1)*rtab(nn+1-i2)*ytab(i2+1)
         enddo
c update fac to next even order
         fac=fac/sqrt(4.*(nn+1)*(nn+2))
      enddo
      roundshapeletpix=zz
      return
      end

      subroutine htab(x,nn,t)
c make table t of Hermite polynomial evaluated at x, for order 0..n
      real t(nn+1)
      t(1)=1
      if (nn.gt.0) t(2)=2*x
c Hn is stored in t(n+1)
c recurrence Hn+1 = 2x Hn - 2n Hn-1
      do n=1,nn-1
          i=n+1
          t(n+2)=2*x*t(n+1) - 2*n*t(n)
      enddo
c      write(0,*) 'HTAB: ',x,n
c      write(0,*) (t(i),i=1,n+1)
      return
      end



