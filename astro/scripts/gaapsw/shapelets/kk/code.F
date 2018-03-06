
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





c routines to encode and decode images into shapelets

      subroutine decode(shape,z,mz,nzin)

      real shape(msh),z(mz,mz),xtab(mm+1),ytab(mm+1)

      nz=abs(nzin)
      beta=shape(1)
      n=shape(2)
      k=2
c      write(0,*) 'Shape:',(shape(i),i=1,2+(n+1)*(n+2)/2)
c make image in z(1:nz,1:nz)
      ctr=0.5*(nz+1)
      s00=0
      s10=0
      s01=0
      s20=0
      s11=0
      s02=0
      do ix=1,nz
          xnorm=(ix-ctr)/beta
          call htab(xnorm,n,xtab)
          ex=exp(-0.5*xnorm**2)
          do iy=1,nz
              ynorm=(iy-ctr)/beta
              call htab(ynorm,n,ytab)
              ey=exp(-0.5*ynorm**2)
              z(ix,iy)=0.
              fac=ex*ey/sqrt(3.1415926535)/beta
c              write(0,*) ix,iy,xnorm,ynorm,fac
              k=2
              do nn=0,n
                  do j=0,nn
                      i=nn-j
                      k=k+1
                      if (j.ne.0) then
                          fac=fac/sqrt(j/real(i+1))
                      else 
                          if (i.ne.0) fac=fac/sqrt(real(i))
                      endif
                      z(ix,iy)=z(ix,iy)+fac*shape(k)*xtab(i+1)*ytab(j+1)
                  enddo
                  fac=fac/sqrt(2.)
              enddo
              zz=z(ix,iy)
              s00=s00+zz
              s10=s10+xnorm*zz*beta
              s01=s01+ynorm*zz*beta
              s20=s20+xnorm**2*zz*beta**2
              s11=s11+xnorm*ynorm*zz*beta**2
              s02=s02+ynorm**2*zz*beta**2
          enddo
      enddo
      if (nzin.lt.0) then
          write(0,*) 'Decoded shape into matrix of order ',nz
          write(0,'(''Moments: '',6f10.2)') s00,s10,s01,s20,s11,s02
          write(0,'(''Normalized:        '',5f10.2)') s10/s00,s01/s00,s20/s00,s11/s00,s02/s00
          write(0,*) 'Series sum flux=',seriessum(shape)
          write(0,*)
      endif
      return
      end




      subroutine decodectr(shape,z,mz,nzin,xc,yc)

c     write out the shape as a pixel array in z(mz,mz) using the first
c     |nzin| pixels of each dimension, and pix pos xc,yc as center 

      real shape(msh),z(mz,mz),xtab(mm+1),ytab(mm+1)

      nz=abs(nzin)
      beta=shape(1)
      n=shape(2)
      k=2
c      write(0,*) 'Shape:',(shape(i),i=1,2+(n+1)*(n+2)/2)
c make image in z(1:nz,1:nz)
      s00=0
      s10=0
      s01=0
      s20=0
      s11=0
      s02=0
      do ix=1,nz
          xnorm=(ix-xc)/beta
          call htab(xnorm,n,xtab)
          ex=exp(-0.5*xnorm**2)
          do iy=1,nz
              ynorm=(iy-yc)/beta
              call htab(ynorm,n,ytab)
              ey=exp(-0.5*ynorm**2)
              z(ix,iy)=0.
              fac=ex*ey/sqrt(3.1415926535)/beta
c              write(0,*) ix,iy,xnorm,ynorm,fac
              k=2
              do nn=0,n
                  do j=0,nn
                      i=nn-j
                      k=k+1
                      if (j.ne.0) then
                          fac=fac/sqrt(j/real(i+1))
                      else 
                          if (i.ne.0) fac=fac/sqrt(real(i))
                      endif
                      z(ix,iy)=z(ix,iy)+fac*shape(k)*xtab(i+1)*ytab(j+1)
                  enddo
                  fac=fac/sqrt(2.)
              enddo
              zz=z(ix,iy)
              s00=s00+zz
              s10=s10+xnorm*zz*beta
              s01=s01+ynorm*zz*beta
              s20=s20+xnorm**2*zz*beta**2
              s11=s11+xnorm*ynorm*zz*beta**2
              s02=s02+ynorm**2*zz*beta**2
          enddo
      enddo
      if (nzin.lt.0) then
          write(0,*) 'Decoded shape into matrix of order ',nz
          write(0,'(''Moments: '',6f10.2)') s00,s10,s01,s20,s11,s02
          write(0,'(''Normalized:        '',5f10.2)') s10/s00,s01/s00,s20/s00,s11/s00,s02/s00
          write(0,*) 'Series sum flux=',seriessum(shape)
          write(0,*)
      endif
      return
      end





      function seriessum(shape)
      real shape(msh)
      double precision sumtot

      beta=shape(1)
      n=shape(2)
      sumtot=0
      faci=1
      do i=0,n,2
c          if (i.gt.0) faci=faci*i*(i-1)/(i/2.)/(i/2.)/4.
          if (i.gt.0) faci=faci*(i-1)/i
          facj=1
          do j=0,n-i,2
c              if (j.gt.0) facj=facj*j*(j-1)/(j/2.)/(j/2.)/4.
              if (j.gt.0) facj=facj*(j-1)/j
              sumtot=sumtot+sqrt(faci*facj)*shape(3+(i+j)*(i+j+1)/2+j)
          enddo
      enddo
      seriessum=sumtot*sqrt(3.1415926535)*2*beta
      return
      end

      function seriessumerr(shape,sherr,sumerr)
c     evaluate the integral of the shapelet series, and its
c     error. sherr=error on each shapelet coef
      real shape(msh)
      double precision sumtot,vartot

      beta=shape(1)
      n=shape(2)
      sumtot=0
      vartot=0
      faci=1
      do i=0,n,2
c          if (i.gt.0) faci=faci*i*(i-1)/(i/2.)/(i/2.)/4.
          if (i.gt.0) faci=faci*(i-1)/i
          facj=1
          do j=0,n-i,2
c              if (j.gt.0) facj=facj*j*(j-1)/(j/2.)/(j/2.)/4.
              if (j.gt.0) facj=facj*(j-1)/j
              sumtot=sumtot+sqrt(faci*facj)*shape(3+(i+j)*(i+j+1)/2+j)
              vartot=vartot+faci*facj
          enddo
      enddo
      seriessumerr=sumtot*sqrt(3.1415926535)*2*beta
      sumerr=sqrt(vartot*3.1415926535)*2*beta*sherr
      return
      end


      function seriessumgauswt(shape,q)
      real shape(msh)
      double precision sumtot

c computes the sum of the shape times a 2-D gaussian of dispersion q

      if (q.le.0.) then
          seriessumgauswt=0
          return
      endif

      beta=shape(1)
      n=shape(2)
      sumtot=0
      faci=1
      b2=beta*beta
      q2=q*q
      t=(q2-b2)/(q2+b2)
      do i=0,n,2
          if (i.gt.0) faci=faci*sqrt((i-1.)/i)*t
          facj=1
          do j=0,n-i,2
              if (j.gt.0) facj=facj*sqrt((j-1.)/j)*t
              sumtot=sumtot+faci*facj*shape(3+(i+j)*(i+j+1)/2+j)
          enddo
      enddo
      seriessumgauswt=sumtot*sqrt(3.1415926535)*2*beta/(1+b2/q2)
      return
      end


      subroutine shwtmoments(sh,w0,wx,wy,wxx,wxy,wyy)
c compute weighted moments of shapelet sh using Gaussian wt fn of size beta
c because of orthogonality relations this depends on coefs up to 2nd order only
      real shape(msh)
      parameter (pi=3.1415926535)
      beta=sh(1)
      s00=sh(3)
      s10=sh(4)
      s01=sh(5)
      s20=sh(6)
      s11=sh(7)
      s02=sh(8)
      w0=beta*sqrt(pi)*s00
      wx=beta**2*sqrt(pi)*sqrt(0.5)*s10
      wy=beta**2*sqrt(pi)*sqrt(0.5)*s01
      wxx=beta**3*sqrt(pi)*(0.5*s00+sqrt(0.5)*s20)
      wxy=beta**3*sqrt(pi)*0.5*s11
      wyy=beta**3*sqrt(pi)*(0.5*s00+sqrt(0.5)*s02)
      return
      end


      subroutine putinimage(shape,xc,yc,radius,bg,pix,nx,ny)
      real shape(msh),pix(mmpix,mmpix),xtab(mm+1),ytab(mm+1)

      ix1=max(1,int(xc-radius))
      ix2=min(nx,int(xc+radius+1))
      iy1=max(1,int(yc-radius))
      iy2=min(ny,int(yc+radius+1))
      beta=shape(1)
      n=shape(2)
      do ix=ix1,ix2
         x=ix-xc
         xnorm=x/beta
         ex=exp(-0.5*xnorm**2)
         call htab(xnorm,n,xtab)
         do iy=iy1,iy2
            y=iy-yc
            if (x*x+y*y .le. radius*radius) then
               ynorm=y/beta
               call htab(ynorm,n,ytab)
               ey=exp(-0.5*ynorm**2)
               sum=0.
               fac=ex*ey/sqrt(3.1415926535)/beta
c              write(0,*) ix,iy,xnorm,ynorm,fac
               k=2
               do nn=0,n
                  do j=0,nn
                     i=nn-j
                     k=k+1
                     if (j.ne.0) then
                        fac=fac/sqrt(j/real(i+1))
                     else 
                        if (i.ne.0) fac=fac/sqrt(real(i))
                     endif
                     sum=sum+fac*shape(k)*xtab(i+1)*ytab(j+1)
                  enddo
                  fac=fac/sqrt(2.)
               enddo
               pix(ix,iy)=bg+sum
            endif
         enddo
      enddo
      return
      end




      subroutine addtoimage(shape,xc,yc,radius,pix,nx,ny)
      real shape(msh),pix(nx,ny),xtab(mm+1),ytab(mm+1)

      ix1=max(1,int(xc-radius))
      ix2=min(nx,int(xc+radius+1))
      iy1=max(1,int(yc-radius))
      iy2=min(ny,int(yc+radius+1))
      beta=shape(1)
      n=shape(2)
      do ix=ix1,ix2
         x=ix-xc
         xnorm=x/beta
         ex=exp(-0.5*xnorm**2)
         call htab(xnorm,n,xtab)
         do iy=iy1,iy2
            y=iy-yc
            if (x*x+y*y .le. radius*radius) then
               ynorm=y/beta
               call htab(ynorm,n,ytab)
               ey=exp(-0.5*ynorm**2)
               sum=0.
               fac=ex*ey/sqrt(3.1415926535)/beta
c              write(0,*) ix,iy,xnorm,ynorm,fac
               k=2
               do nn=0,n
                  do j=0,nn
                     i=nn-j
                     k=k+1
                     if (j.ne.0) then
                        fac=fac/sqrt(j/real(i+1))
                     else 
                        if (i.ne.0) fac=fac/sqrt(real(i))
                     endif
                     sum=sum+fac*shape(k)*xtab(i+1)*ytab(j+1)
                  enddo
                  fac=fac/sqrt(2.)
               enddo
               pix(ix,iy)=pix(ix,iy)+sum
            endif
         enddo
      enddo
      return
      end




      subroutine encode(z,mz,nx,ny,xc,yc,robj,bgin,shape,n,betain)
      real z(mz,mz),shape(msh),xtab(mm+1),ytab(mm+1)
c encode an array z into shapelet coefs of scale beta, to order n
c do this by direct summation (no fit).
      do i=1,msh
          shape(i)=0
      enddo

      if (bgin.eq.0.) then
          call fitbg(z,mz,nx,ny,xc,yc,robj,bg)
      else
          bg=bgin
      endif
      if (betain.le.0.) then
          call fitbeta(z,mz,nx,ny,xc,yc,robj,bg,beta)
          beta=max(beta,-betain)
      else
          beta=betain
      endif

      shape(1)=beta
      shape(2)=n
      robjnorm2=(robj/beta)**2
      ix1=max(1,int(xc-robj))
      ix2=min(nx,int(xc+robj)+1)
      iy1=max(1,int(yc-robj))
      iy2=min(ny,int(yc+robj)+1)
      do ix=ix1,ix2
          xnorm=(ix-xc)/beta
          call htab(xnorm,n,xtab)
          do iy=iy1,iy2
              ynorm=(iy-yc)/beta
              rnorm2=xnorm*xnorm+ynorm*ynorm
              if (rnorm2.gt.robjnorm2) goto 10
              call htab(ynorm,n,ytab)
              fac=exp(-0.5*rnorm2)/sqrt(3.1415926535)/beta
              zz=z(ix,iy)-bg
c              write(0,*) ix,iy,xnorm,ynorm,fac
              k=2
              do nn=0,n
                  do j=0,nn
                      i=nn-j
                      k=k+1
                      if (j.ne.0) then
                          fac=fac/sqrt(j/real(i+1))
                      else 
                          if (i.ne.0) fac=fac/sqrt(real(i))
                      endif
                      shape(k)=shape(k)+zz*fac*xtab(i+1)*ytab(j+1)
                  enddo
                  fac=fac/sqrt(2.)
              enddo
   10     enddo
      enddo
      write(0,*) 'Encoded matrix into shapelets of scale, order ',beta,n
      return
      end


      subroutine fitshapelets(z,mz,nx,ny,xc,yc,robj,bgin,shape,nin,betain)
      real z(mz,mz),shape(msh),xtab(mm+1),ytab(mm+1)
c encode an array z into shapelet coefs of scale beta, to order n
c do this by fitting. If beta or bg are equal to zero, their value
c is fitted. In case of beta, the fitted beta will be at least |betain|
      logical makematrix
      real sksl(msh-2,msh-2),bvec(msh-2),model(msh-2)
      integer indx(msh-2)
      common /sfitshapemat/ sksl,indx,betafit,dxpix,dypix,rfit,nfit

      
      n=abs(nin)
c for nin<0, report inputs
      if (nin.lt.0) then
         write(0,*) 'FITSHAPELETS:'
         write(0,*) 'mz,nx,ny,xc,yc,robj,bgin,nin,betain:'
         write(0,*) mz,nx,ny,xc,yc,robj,bgin,nin,betain
      endif

c for bg exactly equal to zero, fit a value
      if (bgin.eq.0.) then
          call fitbg(z,mz,nx,ny,xc,yc,2*robj,bg)
          if (nin.lt.0) write(0,*) 'Fitted background ',bg
      else
          bg=bgin
      endif
c for negative betain, fit a value. Force result as least as large as |betain|
      if (betain.le.0.) then
          call fitbeta(z,mz,nx,ny,xc,yc,robj,bg,beta)
          beta=max(beta,-betain)
          if (nin.lt.0) write(0,*) 'Fitted beta=',beta
      else
          beta=betain
      endif

      dx=xc-nint(xc)
      dy=yc-nint(yc)
      makematrix=(betafit.ne.beta .or. rfit.ne.robj .or. nfit.ne.n 
     :    .or. abs(dx-dxpix).gt.1e-3 .or. abs(dy-dypix).gt.1e-3)
      if (makematrix) then
          do k=1,msh-2
              do l=1,msh-2
                  sksl(k,l)=0
              enddo
          enddo
          if (nin.lt.0) write(0,*) 'Building matrix again'
      endif
      do k=1,msh-2
          bvec(k)=0
          shape(k+2)=0
      enddo
      nsh=(n+1)*(n+2)/2
      shape(1)=beta
      shape(2)=n
      robjnorm2=(robj/beta)**2
      ix1=max(1,int(xc-robj))
      ix2=min(nx,int(xc+robj)+1)
      iy1=max(1,int(yc-robj))
      iy2=min(ny,int(yc+robj)+1)
      if (nin.lt.0) write(*,*) 'Fitting x, y indices',ix1,' to ',ix2,'; ',iy1,' to ',iy2
      do ix=ix1,ix2
          xnorm=(ix-xc)/beta
          if (xnorm*xnorm.gt.robjnorm2) goto 11
          call htab(xnorm,n,xtab)
          do iy=iy1,iy2
              ynorm=(iy-yc)/beta
              rnorm2=xnorm*xnorm+ynorm*ynorm
              if (rnorm2.gt.robjnorm2) goto 10
              call htab(ynorm,n,ytab)
              fac=exp(-0.5*rnorm2)/sqrt(3.1415926535)/beta
              k=0
              do nn=0,n
                  do j=0,nn
                      i=nn-j
                      k=k+1
                      if (j.ne.0) then
                          fac=fac/sqrt(j/real(i+1))
                      else 
                          if (i.ne.0) fac=fac/sqrt(real(i))
                      endif
                      model(k)=fac*xtab(i+1)*ytab(j+1)
                  enddo
                  fac=fac/sqrt(2.)
              enddo
              do k=1,nsh
                  if (makematrix) then
                      do l=1,nsh
                          sksl(k,l)=sksl(k,l)+model(k)*model(l)
                      enddo
                  endif
                  bvec(k)=bvec(k)+(z(ix,iy)-bg)*model(k)
              enddo
   10     enddo
   11 enddo
      if (makematrix) then
c          do l=nsh-2,1,-1
c              write(0,*) (nint(10*sksl(k,l)/sqrt(sksl(k,k)*sksl(l,l))),k=1,nsh-2)
c          enddo
          call ludcmp(sksl,nsh,msh-2,indx,det)
          betafit=beta
          rfit=robj
          nfit=n
          dxpix=dx
          dypix=dy
          if (nin.lt.0) 
     :        write(0,*) 'New mat for beta,r,n,dx,dy=',
     :        betafit,rfit,nfit,dxpix,dypix
      endif
      call lubksb(sksl,nsh,msh-2,indx,bvec)
      do k=1,nsh
          shape(k+2)=bvec(k)
      enddo
      if (nin.lt.0) 
     :    write(0,*) 'Fitted shapelets of scale, order ',beta,n
      return
      end


      subroutine fitshapeletserr(z,zerr,mz,nx,ny,xc,yc,robj,bg,
     :    shape,sherr,nin,betain)
      real z(mz,mz),zerr(mz,mz),shape(msh),xtab(mm+1),ytab(mm+1)
c encode an array z into shapelet coefs of scale beta, to order n
c do this by fitting. If beta or bg are equal to zero, their value
c is fitted. In case of beta, the fitted beta will be at least |betain|
c Propagate the measurement errors given in the zerr array.
c covariances between coefs are ignored.
c the fitting matrix is always rebuilt

      real sksl(msh-2,msh-2),bvec(msh-2),model(msh-2)
      integer indx(msh-2)

      
      n=abs(nin)
c for bg exactly equal to zero, fit a value
      if (bg.eq.0.) then
          call fitbg(z,mz,nx,ny,xc,yc,2*robj,bg)
          if (nin.lt.0) write(0,*) 'Fitted background ',bg
      endif
c for negative betain, fit a value. Force result as least as large as |betain|
      if (betain.le.0.) then
          call fitbeta(z,mz,nx,ny,xc,yc,robj,bg,beta)
          beta=max(beta,-betain)
      else
          beta=betain
      endif

      dx=xc-nint(xc)
      dy=yc-nint(yc)
      do k=1,msh-2
          do l=1,msh-2
              sksl(k,l)=0
          enddo
      enddo
      do k=1,msh-2
          bvec(k)=0
          shape(k+2)=0
      enddo
      nsh=(n+1)*(n+2)/2
      shape(1)=beta
      shape(2)=n
      robjnorm2=(robj/beta)**2
      ix1=max(1,int(xc-robj))
      ix2=min(nx,int(xc+robj)+1)
      iy1=max(1,int(yc-robj))
      iy2=min(ny,int(yc+robj)+1)
c     at each pixel that is to be fitted, construct the value model(k)
c     of each shapelet at that pixel, and build up the matrix sksl and
c     vector bvec needed for the least-squares fitting.
      do ix=ix1,ix2
          xnorm=(ix-xc)/beta
          if (xnorm*xnorm.gt.robjnorm2) goto 11
          call htab(xnorm,n,xtab)
          do iy=iy1,iy2
              if (zerr(ix,iy).le.0.)  goto 10
              weight=1/zerr(ix,iy)**2
              ynorm=(iy-yc)/beta
              rnorm2=xnorm*xnorm+ynorm*ynorm
              if (rnorm2.gt.robjnorm2) goto 10
              call htab(ynorm,n,ytab)
              fac=exp(-0.5*rnorm2)/sqrt(3.1415926535)/beta
              k=0
              do nn=0,n
                  do j=0,nn
                      i=nn-j
                      k=k+1
                      if (j.ne.0) then
                          fac=fac/sqrt(j/real(i+1))
                      else 
                          if (i.ne.0) fac=fac/sqrt(real(i))
                      endif
                      model(k)=fac*xtab(i+1)*ytab(j+1)
                  enddo
                  fac=fac/sqrt(2.)
              enddo
              do k=1,nsh
                  do l=1,nsh
                      sksl(k,l)=sksl(k,l)+model(k)*model(l)*weight
                  enddo
                  bvec(k)=bvec(k)+(z(ix,iy)-bg)*model(k)*weight
              enddo
   10     enddo
   11 enddo
c solve sksl * shape = bvec
      call ludcmp(sksl,nsh,msh-2,indx,det)
      call lubksb(sksl,nsh,msh-2,indx,bvec)
      do k=1,nsh
          shape(k+2)=bvec(k)
      enddo
c errors on all coefs should be practically identical.
c use variance of 00 component as indication
      bvec(1)=1
      do k=2,nsh
          bvec(k)=0
      enddo
      call lubksb(sksl,nsh,msh-2,indx,bvec)
      if (bvec(1).gt.0.) then
          sherr=sqrt(bvec(1))
      else
          sherr=1.e30
      endif
      if (nin.lt.0) 
     :    write(0,*) 'Fitted shapelets of scale, order ',beta,n
      return
      end


      subroutine plotfittedshapelet(z,mz,nx,ny,xc,yc,robj,bg,
     :    shape)
      real z(mz,mz),shape(msh),xtab(mm+1),ytab(mm+1)

c write out table of model and data for object at xc,yc
c 

      real model(msh-2)
      
      beta=shape(1)
      n=shape(2)
c for bg exactly equal to zero, fit a value
      if (bg.eq.0.) then
          call fitbg(z,mz,nx,ny,xc,yc,2*robj,bg)
      endif

      dx=xc-nint(xc)
      dy=yc-nint(yc)
      nsh=(n+1)*(n+2)/2+1
      robjnorm2=(robj/beta)**2
      ix1=max(1,int(xc-robj))
      ix2=min(nx,int(xc+robj)+1)
      iy1=max(1,int(yc-robj))
      iy2=min(ny,int(yc+robj)+1)
c     at each pixel that is to be fitted, construct the value model(k)
c     of each shapelet at that pixel, and build up the shapelet object.
      do ix=ix1,ix2
          xnorm=(ix-xc)/beta
          if (xnorm*xnorm.gt.robjnorm2) goto 11
          call htab(xnorm,n,xtab)
          do iy=iy1,iy2
              ynorm=(iy-yc)/beta
              rnorm2=xnorm*xnorm+ynorm*ynorm
              if (rnorm2.gt.robjnorm2) goto 10
              call htab(ynorm,n,ytab)
              fac=exp(-0.5*rnorm2)/sqrt(3.1415926535)/beta
              k=0
              fitted=0
              do nn=0,n
                  do j=0,nn
                      i=nn-j
                      k=k+1
                      if (j.ne.0) then
                          fac=fac/sqrt(j/real(i+1))
                      else 
                          if (i.ne.0) fac=fac/sqrt(real(i))
                      endif
                      model(k)=fac*xtab(i+1)*ytab(j+1)
                      fitted=fitted+model(k)*shape(k+2)
                  enddo
                  fac=fac/sqrt(2.)
              enddo
              write(15,*) ix,iy,sqrt(rnorm2)*beta,z(ix,iy)-bg,fitted,z(ix,iy)-bg-fitted
   10     enddo
   11 enddo
      return
      end


      subroutine fitshapeletscov(z,zerr,mz,nx,ny,xc,yc,robj,bg,
     :    shape,sherr,shcov,nin,betain)
      real z(mz,mz),zerr(mz,mz),shape(msh),xtab(mm+1),ytab(mm+1)
      real sherr(msh-2),shcov(msh-2,msh-2)
c encode an array z into shapelet coefs of scale beta, to order n
c do this by fitting. If beta or bg are equal to zero, their value
c is fitted. In case of beta, the fitted beta will be at least |betain|
c Propagate the measurement errors given in the zerr array.
c fully propagate photon errors into covariance matrix and error vector

      real sksl(msh-2,msh-2),bvec(msh-2),model(msh-2)
      integer indx(msh-2)

      
      n=abs(nin)
c for bg exactly equal to zero, fit a value
      if (bg.eq.0.) then
          call fitbg(z,mz,nx,ny,xc,yc,robj,bg)
          if (nin.lt.0) write(0,*) 'Fitted background ',bg
      endif
c for negative betain, fit a value. Force result as least as large as |betain|
      if (betain.le.0.) then
          call fitbeta(z,mz,nx,ny,xc,yc,robj,bg,beta)
          beta=max(beta,-betain)
      else
          beta=betain
      endif

      dx=xc-nint(xc)
      dy=yc-nint(yc)
      do k=1,msh-2
          do l=1,msh-2
              sksl(k,l)=0
          enddo
      enddo
      do k=1,msh-2
          bvec(k)=0
          shape(k+2)=0
      enddo
      nsh=(n+1)*(n+2)/2
      shape(1)=beta
      shape(2)=n
      robjnorm2=(robj/beta)**2
      ix1=max(1,int(xc-robj))
      ix2=min(nx,int(xc+robj)+1)
      iy1=max(1,int(yc-robj))
      iy2=min(ny,int(yc+robj)+1)
c     at each pixel that is to be fitted, construct the value model(k)
c     of each shapelet at that pixel, and build up the matrix sksl and
c     vector bvec needed for the least-squares fitting.
      do ix=ix1,ix2
          xnorm=(ix-xc)/beta
          if (xnorm*xnorm.gt.robjnorm2) goto 11
          call htab(xnorm,n,xtab)
          do iy=iy1,iy2
              if (zerr(ix,iy).le.0.)  goto 10
              weight=1/zerr(ix,iy)**2
              ynorm=(iy-yc)/beta
              rnorm2=xnorm*xnorm+ynorm*ynorm
              if (rnorm2.gt.robjnorm2) goto 10
              call htab(ynorm,n,ytab)
              fac=exp(-0.5*rnorm2)/sqrt(3.1415926535)/beta
              k=0
              do nn=0,n
                  do j=0,nn
                      i=nn-j
                      k=k+1
                      if (j.ne.0) then
                          fac=fac/sqrt(j/real(i+1))
                      else 
                          if (i.ne.0) fac=fac/sqrt(real(i))
                      endif
                      model(k)=fac*xtab(i+1)*ytab(j+1)
                  enddo
                  fac=fac/sqrt(2.)
              enddo
              do k=1,nsh
                  do l=1,nsh
                      sksl(k,l)=sksl(k,l)+model(k)*model(l)*weight
                  enddo
                  bvec(k)=bvec(k)+(z(ix,iy)-bg)*model(k)*weight
              enddo
   10     enddo
   11 enddo
c solve sksl * shape = bvec
      call ludcmp(sksl,nsh,msh-2,indx,det)
      call lubksb(sksl,nsh,msh-2,indx,bvec)
      do k=1,nsh
          shape(k+2)=bvec(k)
      enddo
c recover the full error array and cov matrix---that is just sksl **-1)      

      do k=1,nsh
          do l=1,nsh
              bvec(l)=0
          enddo
          bvec(k)=1
          call lubksb(sksl,nsh,msh-2,indx,bvec)
          sherr(k)=sqrt(bvec(k))
          do l=1,nsh
              shcov(k,l)=bvec(l)
          enddo
      enddo
      if (nin.lt.0) 
     :    write(0,*) 'Fitted shapelets of scale, order ',beta,n
      return
      end


      subroutine to012(sh,sh012,norder)
c     convert a shapelet sh into a polar shapelet containing the m=0,1,2
c     coefficients for radial orders up to norder. Used specifically with
c     sh2gerrshift012 and fitshearshift_err_012.

      parameter (mm012=5*(mm/2)-1)
      real sh(msh),sh012(mm012)

c structure of sh012: first the m=0 terms of order 0,2,4,...N-2;
c then the cos(phi) terms of order 1,3,5,...,N-3 followed by the sine terms
c then the cos(2phi) terms of order 2,4,...N-4 followed by the sine terms
      
      n=nint(sh(2))
      if (norder.gt.n) then
         write(0,*) 'IN TO012: REQUESTED ORDER EXCEEDS SHAPELET ORDER'
         stop
      endif

c first construct the axisymmetric components (leading order r^2m)
      do nn=0,norder,2
         m=nn/2
         k0=3+nn*(nn+1)/2
         cc=1
         sum2=cc**2
         sumcoef=cc*sh(k0)
c         write(0,*) 'rnd',nn,0,cc
         do i=1,m
            i2=2*i
            mi2=2*(m-i)
            cc=cc*(m-i+1)/i*sqrt(i2*(i2-1)/(mi2+2.)/(mi2+1.))
            sum2=sum2+cc**2
            sumcoef=sumcoef+cc*sh(k0+i2)
c            write(0,*) 'rnd',nn,i2,cc
         enddo
         sh012(m+1)=sumcoef/sqrt(sum2)
      enddo
c now the cos-phi cpts
      do nn=1,norder-1,2
         m=(nn-1)/2
         k0=3+nn*(nn+1)/2
         cc=1
         sum2=cc**2
         sumcoef=cc*sh(k0)
c         write(0,*) 'cos-th',nn,0,cc
         do i=1,m
            i2=2*i
            mi2=2*(m-i)
            cc=cc*(m-i+1)/i*sqrt(i2*(i2-1)/(mi2+3.)/(mi2+2.))
            sum2=sum2+cc**2
            sumcoef=sumcoef+cc*sh(k0+i2)
c            write(0,*) 'cos-th',nn,i2,cc
         enddo
         sh012((norder/2)+1+m+1)=sumcoef/sqrt(sum2)
      enddo
c now the sin-phi cpts
      do nn=1,norder-1,2
         m=(nn-1)/2
         k0=3+nn*(nn+1)/2
         cc=1
         sum2=cc**2
         sumcoef=cc*sh(k0+1)
c         write(0,*) 'sin-th',nn,1,cc
         do i=1,m
            i2=2*i
            mi2=2*(m-i)
            cc=cc*(m-i+1)/i*sqrt((i2+1)*i2/(mi2+2.)/(mi2+1.))
            sum2=sum2+cc**2
            sumcoef=sumcoef+cc*sh(k0+i2+1)
c            write(0,*) 'sin-th',nn,i2+1,cc
         enddo
         sh012(2*(norder/2)+1+m+1)=sumcoef/sqrt(sum2)
      enddo
c now the cos-2phi cpts
      do nn=2,norder-2,2
         m=(nn-2)/2
         k0=3+nn*(nn+1)/2
         cc=1
         sum2=cc**2
         sumcoef=cc*sh(k0)
c         write(0,*) 'cos-2th',nn,0,cc
         cplus=cc
         fac=1
         do i=1,m+1
            i2=2*i
            mi2=2*(m-i)
            cmin=cplus
            cplus=cplus*(m-i+1)/i
            fac=fac*(i2-1)*i2/(mi2+4.)/(mi2+3.)
            cc=(cplus-cmin)*sqrt(fac)
            sum2=sum2+cc**2
            sumcoef=sumcoef+cc*sh(k0+i2)
c            write(0,*) 'cos-2th',nn,i2,cc
         enddo
         sh012(3*(norder/2)+1+m+1)=sumcoef/sqrt(sum2)
      enddo
c now the sin-2phi cpts
      do nn=2,norder-2,2
         m=(nn-2)/2
         k0=3+nn*(nn+1)/2
         cc=1
         sum2=cc**2
         sumcoef=cc*sh(k0+1)
c         write(0,*) 'sin-2th',nn,1,cc
         do i=1,m
            i2=2*i
            mi2=2*(m-i)
            cc=cc*(m-i+1)/i*sqrt((i2+1)*i2/(mi2+3.)/(mi2+2.))
            sum2=sum2+cc**2
            sumcoef=sumcoef+cc*sh(k0+i2+1)
c            write(0,*) 'sin-2th',nn,i2+1,cc
         enddo
         sh012(4*(norder/2)+m+1)=sumcoef/sqrt(sum2)
      enddo
c      write(0,*) 'Last coef:',4*(norder/2)+(norder-4)/2+1
c      write(0,*) (sh012(k),k=1,5*(norder/2)-1)

      return
      end



      subroutine fitbg(z,mz,nx,ny,xc,yc,robj,bg)
      parameter (kmax=1000)
      real z(mz,mz),sky(kmax)
c find sky value from an annulus between robj and 2*robj
      robj2=robj*robj
      rout2=2*robj2
      rout=1.41*robj
      if (3.14*(rout2-robj2).gt.0.7*kmax) then
          rout2=0.7*kmax/3.14+robj2
          rout2=max(rout2,(robj+1)**2)
          rout=sqrt(rout2)
      endif
      ix1=max(1, int(xc-rout))
      ix2=min(nx,int(xc+rout)+1)
      iy1=max(1, int(yc-rout))
      iy2=min(ny,int(yc+rout)+1)
      robj2=robj*robj
      k=0
      do ix=ix1,ix2
          do iy=iy1,iy2
              r2=(ix-xc)**2+(iy-yc)**2
              if (r2.gt.robj2 .and. r2.lt.rout2) then
                  if (k.eq.kmax) goto 1
                  k=k+1
                  sky(k)=z(ix,iy)
              endif
          enddo
      enddo
    1 bg=fmedian(k,sky)
      return
      end


      subroutine fitbeta(z,mz,nx,ny,xc,yc,robj,bg,beta)
      real z(mz,mz),coef(2),covar(2,2),alpha(2,2)
      parameter (nmaxfit=10000)
      real r2(nmaxfit),sig(nmaxfit),tofit(nmaxfit)
      external fgaus
c fit (flx-bg) vs radius for all pixels inside radius robj
c k is the counter of the number of pixels to be fitted
      k=0
      zmax=-1e30
      ix1=max(1,int(xc-robj))
      ix2=min(nx,int(xc+robj)+1)
      iy1=max(1,int(yc-robj))
      iy2=min(ny,int(yc+robj)+1)
c check if need to subsample
      nfitapprox=(iy2-iy1+1)*(ix2-ix1+1)
      nsub=1
      if (nfitapprox.gt.nmaxfit) nsub=1+int(sqrt(real(nfitapprox)/nmaxfit))
      do i=ix1,ix2,nsub
          x=i-xc
          do j=iy1,iy2,nsub
              y=j-yc
              rr=x*x+y*y
c NEXT LINE WAS rr<robj**2/4
              if (rr.le.robj*robj) then
                  k=k+1
                  r2(k)=rr
c downweight large radii - choose weight to ensure equal weight per ring.
                  sig(k)=sqrt(max(1.,rr))
                  tofit(k)=z(i,j)-bg
                  zmax=max(zmax,tofit(k))
c WAS                  zmax=max(zmax,z(i,j))
              endif
          enddo
      enddo
c check dispersion of tofit array to check there is something to fit there
      call sort2(k,tofit,r2)
c want difference between 10th and 90th percentile - otherwise too many junk values
      if (k.le.10) then
         k10=1
         k90=k
      else
         k10=k/10+1
         k90=k-k10+1
      endif
      if (tofit(k10).eq.tofit(k90)) then
c break out - no sensible Gaussian can be fitted to constant array
          beta=-1
          return
      endif
      if (k.gt.nmaxfit) then
          write(0,*) '*** PROBLEM IN FITBETA --- SUBSAMPLING FAILED ***'
          stop
      endif
c      call pgenv(-0.1*robj,robj*1.1,-0.1*zmax,1.1*zmax,0,0)
c      do i=1,k
c          call pgpt(1,sqrt(r2(i)),tofit(i),5)
c      enddo
      alam=-1
      coef(1)=zmax
c WAS      coef(1)=zmax*3
      coef(2)=0.01
c      write(0,*) alam,chisq,coef
      do iter=1,10
          call mrqmin(r2,tofit,sig,k,coef,2,covar,alpha,2,chisq,fgaus,alam)
c          write(0,*) alam,chisq,coef
c escape if you're fitting a point spike (sigma < 0.2 pix)
          if (coef(2).gt.12.5) goto 3
      enddo
    3 if (coef(2).gt.0.) then
          beta=(2*coef(2))**(-0.5)
c          write(0,*) 'Fitted a beta value of ',beta
      else
          beta=1.
c          write(0,*) '*** could not find a good beta. Assumed beta=1'
      endif
      return
      end

      subroutine centershape(sh,sh2,dx,dy)
c shift the shape so that the 10 and 01 components are zero
      real sh(msh),sh2(msh),shx(msh),shy(msh)
      beta=sh(1)
      n=sh(2)
      if (n.eq.0) then
          do k=1,(n+1)*(n+2)/2+2
              sh2(k)=sh(k)
          enddo
          dx=0
          dy=0
          return
      endif
      call trans1(sh,shx)
      call trans2(sh,shy)
      det=shx(4)*shy(5)-shx(5)*shy(4)
      dx=-( sh(4)*shy(5)-sh(5)*shy(4))/det
      dy=-(-sh(4)*shx(5)+sh(5)*shx(4))/det
      sh2(1)=sh(1)
      sh2(2)=sh(2)
      do k=3,(n+1)*(n+2)/2+2
         sh2(k)=sh(k)+dx*shx(k)+dy*shy(k)
      enddo
      dx=dx*beta
      dy=dy*beta
      return
      end


      subroutine fgaus(r2,coef,y,dyda,ma)
      real coef(2),dyda(2)
      ee=exp(-r2*coef(2))
      y=ee*coef(1)
      dyda(1)=ee
      dyda(2)=-r2*y
      return
      end



c write out the shape coefs in a nice table
      subroutine writeshape(text,shape)
      real shape(msh)
      character text*(*)
 1000 format (a15,f6.2,f5.0/f8.4/2f8.4/3f8.4/4f8.4/5f8.4/6f8.4/7f8.4/8f8.4/9f8.4/10f8.4/11f8.4/12f8.4/13f8.4/14f8.4/15f8.4/16f8.4/17f8.4)
      n=shape(2)
      write(0,1000) text,(shape(i),i=1,2+(n+1)*(n+2)/2)
      write(0,*)
      return
      end

c write out the normalized shape coefs in a nice table
      subroutine writeshapen(text,shape)
      real shape(msh)
      character text*(*)
 1000 format (a15,f6.2,f5.0,'  (Normalized to largest coef)'/
     :    f8.4/2f8.4/3f8.4/4f8.4/5f8.4/6f8.4/7f8.4/8f8.4/9f8.4/10f8.4/11f8.4/12f8.4/13f8.4/14f8.4/15f8.4/16f8.4/17f8.4)
      n=shape(2)
      nsh=(n+1)*(n+2)/2
      t=0
      do k=3,2+nsh
          t=max(abs(shape(k)),t)
      enddo
c      t=seriessum(shape)/shape(1)
      write(0,1000) text,shape(1),shape(2),(shape(i)/t,i=3,2+nsh)
      write(0,*)
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



