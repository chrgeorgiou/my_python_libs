c     given fits image inimage.fits and a PSF catalogue, construct a shapelet Gaussianization kernel for each star, and write this out as a shapelet catalogue (for fitting by fitkermap).



      dimension x(mmcat),y(mmcat),f(mmcat),fe(mmcat),fw(mmcat),
     :    a(mmcat),b(mmcat),th(mmcat),id(mmcat),iflg(mmcat),beta(mmcat)
      parameter (maxstamp=200,maxstamp2=maxstamp*maxstamp)
      real stamp(maxstamp2)
      real kercoef(msh)


c open image
      call openfits('inimage.fits',51,nx,ny)

c     read in the target fitting radius.
c     If input negative or zero, determine beta automatically
      read(*,*) betain

c     read all PSF stars
      call readorders(n,nfit)
      do i=1,mmcat
 1       read(*,*,err=1,end=99) 
     :        x(i),y(i),f(i),fe(i),fw(i),a(i),b(i),th(i),id(i),iflg(i)
c skip objects that fall outside the fits image boundaries
         if (x(i).lt.1. .or. x(i).gt.nx*1. .or. y(i).lt.1. .or. y(i).gt.ny*1.) goto 1
      enddo
      write(0,*) 'There may be more objects than the ',mmcat,' read.'
 99   nobj=i-1
      write(0,*) 'Read ',nobj,' stars.'

c     if input beta negative or zero determine optimal value from the image

      if (betain.le.0.) then

c     assign single PSF beta and fitting radius
c     find best-fit beta to all PSF objects, take median, multiply by 1.3
c     the factor 1.3 is result of simulations that show 
c     this helps with moffat profiles, and does no harm for gaussian PSFs.
         
         do i=1,nobj
            fitradius=max(10.,fw(i)/2.3*4)
            npix=min(maxstamp,nint(4*fitradius+3))
            ix1=max(1,min(nx-npix+1,nint(x(i))-(npix-1)/2))
            iy1=max(1,min(ny-npix+1,nint(y(i))-(npix-1)/2))
            isok=igetsec(51,nx,ny,ix1,iy1,npix,npix,stamp,maxstamp2)
            if (isok.eq.0) then
               call fitbg(stamp,npix,npix,npix,x(i)-ix1+1,y(i)-iy1+1,fitradius*2,bg)
               call fitbeta(stamp,npix,npix,npix,x(i)-ix1+1,y(i)-iy1+1,fitradius,bg,beta(i))
               fitradius=max(10.,1.3*beta(i)*4)
               fitradius=min(fitradius,0.25*npix-1)
               call fitbg(stamp,npix,npix,npix,x(i)-ix1+1,y(i)-iy1+1,fitradius*2,bg)
               call fitbeta(stamp,npix,npix,npix,x(i)-ix1+1,y(i)-iy1+1,fitradius,bg,beta(i))
            endif
         enddo
         betapsf=1.3*fmedian(nobj,beta)
      else
         betapsf=betain
      endif
      if (betapsf.lt.0.1*n) betapsf=0.1*n
      write(0,*) 'Adopted betapsf: ',betapsf

      fitradius=max(10.,betapsf*(4+sqrt(1.*n)))
      fitradius=min(maxstamp*0.25-1.,fitradius)
      write(0,*) 'Using fitting radius of',fitradius

c now go through the objects and fit a shapelet kernel that Gaussianizes the PSF
      do i=1,nobj
          if (x(i).gt.fitradius+1 .and. x(i).lt.nx-fitradius
     :        .and. y(i).gt.fitradius+1 .and. y(i).lt.ny-fitradius) then

             npix=min(maxstamp,nint(4*fitradius+3))
             ix1=max(1,min(nx-npix+1,nint(x(i))-(npix-1)/2))
             iy1=max(1,min(ny-npix+1,nint(y(i))-(npix-1)/2))
             isok=igetsec(51,nx,ny,ix1,iy1,npix,npix,stamp,maxstamp2)
             call fitbg(stamp,npix,npix,npix,x(i)-ix1+1,y(i)-iy1+1,2*fitradius,bg)
c     using catalogue centre, make a kernel that produces a gaussian at that position
             call fitkernel(stamp,npix,npix,x(i)-ix1+1,y(i)-iy1+1,fitradius,bg,betapsf,n,kercoef)
c     find kernel centroid and move target PSF by that amount
             call centershape(kercoef,kercoef,dx,dy)
c             write(0,*) i,'Initial kernel shift was',dx,dy 
c     repeat kernel fit with updated centre for target gaussian
             call fitkernel(stamp,npix,npix,x(i)+dx-ix1+1,y(i)+dy-iy1+1,fitradius,bg,betapsf,n,kercoef)
             call centershape(kercoef,kercoef,dx2,dy2)
c             write(0,*) i,'Second kernel shift was',dx2,dy2
             dx=dx+dx2
             dy=dy+dy2
             call fitkernel(stamp,npix,npix,x(i)+dx-ix1+1,y(i)+dy-iy1+1,fitradius,bg,betapsf,n,kercoef)
             call centershape(kercoef,kercoef,dx2,dy2)
c             write(0,*) i,'Third kernel shift was',dx2,dy2 
             dx=dx+dx2
             dy=dy+dy2
c write out the kernel shapelet for this star             
             write(*,1001) x(i),y(i),f(i),fe(i),fw(i),a(i),b(i),th(i),id(i),iflg(i)
 1001        format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i10,i5)
             write(*,*) 1.,dx,dy,bg,fitradius,1.
             write(*,*) betapsf,n
             write(*,1002) (kercoef(k),k=3,(n+1)*(n+2)/2+2)
 1002        format(msh(g15.6))
             write(*,*)
          endif
      enddo
      write(0,'('' Psfcat2sh read'',i7,'' stars, used shapelet scale '',
     :     f6.3,'' and order '',i2)') nobj,betapsf,n
      end

c-------------------------------------------
      
      subroutine fitkernel(stamp,m,n,xc,yc,rfit,bg,beta,nsh,sh)
      real sh(msh),stamp(m,m)
      real z(64,64),model(64,64,msh-2),amat(msh-2,msh-2),cvec(msh-2),xtab(mm+1),ytab(mm+1)
      integer indx(msh-2)
      complex cz(32,64),cmodel(32,64,msh-2),fac,speq(64)
      equivalence (model,cmodel),(z,cz)
      double precision sum
      parameter (ntaper=8,twopi=6.28318530717958)
      logical verbose
      
c     from postage stamp (nxn data in mxm array) fit a gaussianization kernel to star with ctr at xc,yc.
c     kernel is a sum of gaussians - fit of kernel coefs is directly made to the pixels of the star.

c     first make the model functions to superpose linearly: observed PSF
c     convolved with each shapelet mask postage stamp with circular
c     hamming filter, FT, and then multiply by shapelets in turn and
c     inverse-FFT

      verbose=.false.
      
      do j=1,64
         do i=1,64
            z(i,j)=0
         enddo
      enddo

      ixc=nint(xc)
      iyc=nint(yc)
      xc64=xc-ixc+32
      yc64=yc-iyc+32
c copy postage stamp into z, shifting so that ixc,iyc lands on 32,32      
      do j=max(1,1-32+iyc),min(n,64-32+iyc)
         do i=max(1,1-32+ixc),min(n,64-32+ixc)
            z(i+32-ixc,j+32-iyc)=stamp(i,j)-bg
         enddo
      enddo
c     hamming filter outside circle of radius (32.-ntaper)
      do j=1,64
         do i=1,64
            r=sqrt((32.-i)**2+(32.-j)**2)
            if (r.ge.32.) then
               z(i,j)=0
            elseif (r.ge.32.-ntaper) then
               taperfac=(1.-cos(3.1415926535*(32.-r)/ntaper))/2.
               z(i,j)=z(i,j)*taperfac
            endif
         enddo
      enddo
c     FT of postage stamp
      call rlft3(z,speq,64,64,1,1)
c     now multiply by all k shapelets' FTs in turn and inverse FT. results in model(x,y,k)
      
      npix=64
c     construct shapelet function's FT at frequency fxnorm, fynorm
c     Hn(x/beta)exp(-0.5*x**2/beta**2) --> i^n Hn(f*beta) exp(-0.5*f**2
c     *beta**2)

      fmax2=(sqrt(nsh*1.)+3)**2

c make shapelets convolved with stamp in Fourier space i.e. FT of z multiplied by FT of shapelet k
      do iy=1,npix
         fynorm=(mod(iy-1+npix/2,npix)-npix/2)*twopi*beta/npix
         call htab(fynorm,n,ytab)
         do ix=1,npix/2
            fxnorm=(mod(ix-1+npix/2,npix)-npix/2)*twopi*beta/npix
            fnorm2=fxnorm*fxnorm+fynorm*fynorm
            if (fnorm2.gt.fmax2) then
               do k=1,(nsh+1)*(nsh+2)/2
                  cmodel(ix,iy,k)=0
               enddo
            else
               call htab(fxnorm,n,xtab)
               fac=exp(-0.5*fnorm2)/sqrt(3.1415926535)*beta*twopi/(npix*npix/2) * cz(ix,iy)
               k=0
               do nn=0,nsh
                  do j=0,nn
                     i=nn-j
                     k=k+1
                     if (j.ne.0) then
                        fac=fac/sqrt(j/real(i+1))
                     else 
                        if (i.ne.0) fac=fac/sqrt(real(i))
                     endif
                     cmodel(ix,iy,k)=xtab(i+1)*ytab(j+1)*fac*(0.,1.)**nn
                  enddo
                  fac=fac/sqrt(2.)
               enddo
            endif
         enddo
      enddo

c     now inverse-FT all models to give real-space convolutions of stamps with shapelet basis.
      
      nmod=(nsh+1)*(nsh+2)/2
      do k=1,nmod
         do i=1,64
            speq(i)=0
         enddo
         call rlft3(model(1,1,k),speq,64,64,1,-1)
      enddo

      if (verbose) then
         write(0,*) 'Writing out',nmod,'models to display with showim'
         do k=1,nmod
c     write(0,*) 'k=',k,nmod
            write(k+100,*) 'MODEL',k
            write(k+100,*) 64,64
            do j=1,64
               do i=1,64
                  write(k+100,*) i,j,model(i,j,k)
               enddo
            enddo
            close(k+100)
         enddo
      endif
      
c     build the model x model matrix. Only fit pixels out to rfit + beta from the centre.
      rout2=(rfit+beta)**2
      do k=1,nmod
         do l=k,nmod
            sum=0
            do iy=1,64
               do ix=1,64
                  if ((ix-xc64)**2+(iy-yc64)**2 .le. rout2) then
                     sum=sum+model(ix,iy,k)*model(ix,iy,l)
                  endif
               enddo
            enddo
            amat(k,l)=sum
            amat(l,k)=sum
         enddo
c     build the model x target vector (latter is the Gaussian we want to convolve to)
         sum=0
         do iy=1,64
            do ix=1,64
               r2=(ix-xc64)**2+(iy-yc64)**2
               if (r2 .le. rout2) then
                  sum=sum+model(ix,iy,k)*exp(-r2/(2*beta**2))/(twopi*beta**2)
               endif
            enddo
         enddo
         cvec(k)=sum
      enddo
c     now solve for the kernel coeffs.
      call ludcmp(amat,nmod,msh-2,indx,det)
      call lubksb(amat,nmod,msh-2,indx,cvec)

      if (verbose) then
         write(0,*) (cvec(k),k=1,nmod)
c     write model and residual to fort.100 and fort.99
         write(99,*) 'TARGET'
         write(99,*) 64,64
         write(100,*) 'BEST MODEL - TARGET'
         write(100,*) 64,64
         do iy=1,64
            do ix=1,64
               sum=0
               do k=1,nmod
                  sum=sum+model(ix,iy,k)*cvec(k)
               enddo
               write(100,*) ix,iy,real(sum)
               r2=(ix-xc64)**2+(iy-yc64)**2
               if (r2 .le. rout2) then
                  write(99,*) ix,iy,sum-exp(-r2/(2*beta**2))/(twopi*beta**2)
               else
                  write(99,*) ix,iy,sum
               endif
            enddo
         enddo
         close(99)
         close(100)
      endif
      
c write the shapelet coefs in the sh vector to be returned.
      sh(1)=beta
      sh(2)=nsh
      do k=1,nmod
         sh(k+2)=cvec(k)
      enddo
      sum=seriessum(sh)
      do k=1,nmod
         sh(k+2)=sh(k+2)/sum
      enddo
      return
      end
      
