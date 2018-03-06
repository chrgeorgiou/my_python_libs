
      function psfcorphot(shobs,shpsf,q)
c
c     takes the observed shapelet shobs, and the PSF in shpsf.
c     calculates the gaussian-weighted aperture for the intrinsic source
c     after smearing with a gaussian of scale radius q, and weighting
c     with a gaussian of same size.
c
c     this value is returned as psfcorphot
c

      dimension shobs(msh),shpsf(msh)
      dimension qcoef(0:mm)
      dimension pmat(msh-2,msh-2),pinv(msh-2,msh-2),indx(msh-2)
      parameter (pi=3.1415926535)

c first calculate the PSF matrix from the shapelet PSF supplied.
c use q as scale for shapelet to operate on, and for result.
c setpsfmat: PSF shape, result scale, input scale, matrix

      if (q.lt.0.) then
          psfcorphot=0
          return
      endif

      betaobs=shobs(1)
      n=nint(shobs(2))
      nsh=(n+1)*(n+2)/2
      call setpsfmat(shpsf,betaobs,q,pmat)

c     now invert the matrix; result pinv need all rows, but not all
c     columns; unfortunately this means that full matrix must be
c     inverted. If could work with transpose of pmat then fewer (about
c     half) columns need to be calculated.
c     
c     Result is the deconvolved source, expressed as a shapelet series
c     of scale q

      do i=1,nsh
         do j=1,nsh
            pinv(i,j)=0
         enddo
         pinv(i,i)=1
      enddo
      call ludcmp(pmat,nsh,msh-2,indx,detfac)

      do k=1,nsh
          call lubksb(pmat,nsh,msh-2,indx,pinv(1,k))
      enddo
      
c     tabulate the elementary weighted fluxes: qcoef(n) is the flux of a
c     1-D shapelet of order n, scale q, that is convolved with a
c     normalized gaussian of scale q, and for a gaussian-weighted flux
c     (also scale q).These can then be summed over all coefs of shobs.

      qcoef(0)=sqrt(2*q*sqrt(pi)/3)
      do k=2,n,2
          qcoef(k-1)=0
          qcoef(k)=qcoef(k-2) /(3*k) * sqrt(k*(k-1.))
      enddo

c     finally calculate the total weighted flux, cpt per cpt. Apply the
c     qcoef weights in each dimension to the deconvolved shape
c      PINV . SHOBS

c     could track the error here as well by accumulating coefs of each
c     shobs(l) and taking the rms

      phot00=0

      k=0
      do nn=0,n
          do j=0,nn
              i=nn-j
              k=k+1
              if (mod(i,2)+mod(j,2).eq.0) then
                  do l=1,nsh
                      phot00=phot00+qcoef(i)*qcoef(j)*pinv(k,l)*shobs(l+2)
                  enddo
              endif
          enddo
      enddo

      psfcorphot=phot00

      end



      function psfcorphot2(shobs,shpsf,q)
c
c     takes the observed shapelet shobs, and the PSF in shpsf.
c     calculates the gaussian-weighted aperture for the intrinsic source
c     after smearing with a gaussian of scale radius q, and weighting
c     with a gaussian of same size.
c
c     this value is returned as psfcorphot2
c
c     THIS VERSION SAME AS PSFCORPHOT, BUT SCALE FOR INTRINSIC SOURCE IS
c     TAKEN TO BE betaobs, AND NOT Q.
c     THIS MODIFIES THE SERIES TO BE SUMMED.
c

      dimension shobs(msh),shpsf(msh)
      dimension qcoef(0:mm)
      dimension pmat(msh-2,msh-2),pinv(msh-2,msh-2),indx(msh-2)
      parameter (pi=3.1415926535)

c first calculate the PSF matrix from the shapelet PSF supplied.
c use betaobs as scale for shapelet to operate on, and for result.
c setpsfmat: PSF shape, result scale, input scale, matrix


      if (q.lt.0.) then
          psfcorphot2=0
          return
      endif

      betaobs=shobs(1)
      n=nint(shobs(2))
      nsh=(n+1)*(n+2)/2
      call setpsfmat(shpsf,betaobs,betaobs,pmat)

c     now invert the matrix; result pinv need all rows, but not all
c     columns; unfortunately this means that full matrix must be
c     inverted. If could work with transpose of pmat then fewer (about
c     half) columns need to be calculated.
c     
c     Result is the deconvolved source, expressed as a shapelet series
c     of scale betaobs

      do i=1,nsh
         do j=1,nsh
            pinv(i,j)=0
         enddo
         pinv(i,i)=1
      enddo
      call ludcmp(pmat,nsh,msh-2,indx,detfac)

      do k=1,nsh
          call lubksb(pmat,nsh,msh-2,indx,pinv(1,k))
      enddo
      
c     tabulate the elementary weighted fluxes: qcoef(n) is the flux of a
c     1-D shapelet of order n, scale betaobs, that is convolved with a
c     normalized gaussian of scale q, and for a gaussian-weighted flux
c     (also scale q).These can then be summed over all coefs of shobs.

      q2=q*q
      b2=betaobs*betaobs
      t=(2*q2-b2)/(2*q2+b2)
      qcoef(0)=sqrt(betaobs*sqrt(pi)/(1+0.5*b2/q2))
      do k=2,n,2
          qcoef(k-1)=0
          qcoef(k)=qcoef(k-2) * sqrt((k-1.)/k) * t
      enddo

c     finally calculate the total weighted flux, cpt per cpt. Apply the
c     qcoef weights in each dimension to the deconvolved shape
c      PINV . SHOBS

c     could track the error here as well by accumulating coefs of each
c     shobs(l) and taking the rms

      phot00=0

      k=0
      do nn=0,n
          do j=0,nn
              i=nn-j
              k=k+1
              if (mod(i,2)+mod(j,2).eq.0) then
                  do l=1,nsh
                      phot00=phot00+qcoef(i)*qcoef(j)*pinv(k,l)*shobs(l+2)
                  enddo
              endif
          enddo
      enddo

      psfcorphot2=phot00

      end


      function psfcorphot3(shobs,sherr,shpsf,q,photerr)
c
c     takes the observed shapelet shobs, and the PSF in shpsf.
c     calculates the gaussian-weighted aperture for the intrinsic source
c     after smearing with a gaussian of scale radius q, and weighting
c     with a gaussian of same size.
c
c     this value is returned as psfcorphot3
c
c     THIS VERSION SAME AS PSFCORPHOT, BUT SCALE FOR INTRINSIC SOURCE IS
c     TAKEN TO BE betaobs, AND NOT Q.
c     THIS MODIFIES THE SERIES TO BE SUMMED.
c

      dimension shobs(msh),shpsf(msh),c(msh-2)
      dimension qcoef(0:mm)
      dimension pmat(msh-2,msh-2),pinv(msh-2,msh-2),indx(msh-2)
      double precision phot00,photvar
      parameter (pi=3.1415926535)

c first calculate the PSF matrix from the shapelet PSF supplied.
c use betaobs as scale for shapelet to operate on, and for result.
c setpsfmat: PSF shape, result scale, input scale, matrix


      if (q.lt.0.) then
          psfcorphot3=0
          photerr=1
          return
      endif

      betaobs=shobs(1)
      n=nint(shobs(2))
      nsh=(n+1)*(n+2)/2
      call setpsfmat(shpsf,betaobs,betaobs,pmat)

c     now invert the matrix; result pinv need all rows, but not all
c     columns; unfortunately this means that full matrix must be
c     inverted. If could work with transpose of pmat then fewer (about
c     half) columns need to be calculated.
c     
c     Result is the deconvolved source, expressed as a shapelet series
c     of scale betaobs

      do i=1,nsh
         do j=1,nsh
            pinv(i,j)=0
         enddo
         pinv(i,i)=1
      enddo
      call ludcmp(pmat,nsh,msh-2,indx,detfac)

      do k=1,nsh
          call lubksb(pmat,nsh,msh-2,indx,pinv(1,k))
      enddo
      
c     tabulate the elementary weighted fluxes: qcoef(n) is the flux of a
c     1-D shapelet of order n, scale betaobs, that is convolved with a
c     normalized gaussian of scale q, and for a gaussian-weighted flux
c     (also scale q).These can then be summed over all coefs of shobs.

      q2=q*q
      b2=betaobs*betaobs
      t=(2*q2-b2)/(2*q2+b2)
      qcoef(0)=sqrt(betaobs*sqrt(pi)/(1+0.5*b2/q2))
      do k=2,n,2
          qcoef(k-1)=0
          qcoef(k)=qcoef(k-2) * sqrt((k-1.)/k) * t
      enddo

c     finally calculate the total weighted flux, cpt per cpt. Apply the
c     qcoef weights in each dimension to the deconvolved shape
c      PINV . SHOBS

c     track the error here as well by accumulating coefs of each
c     shobs(l) and taking the rms

      do l=1,nsh
          c(l)=0
      enddo

      k=0
      do nn=0,n
          do j=0,nn
              i=nn-j
              k=k+1
              if (mod(i,2)+mod(j,2).eq.0) then
                  do l=1,nsh
                      c(l)=c(l)+qcoef(i)*qcoef(j)*pinv(k,l)
                  enddo
              endif
          enddo
      enddo
      phot00=0
      photvar=0
      do l=1,nsh
          phot00=phot00+c(l)*shobs(l+2)
          photvar=photvar+c(l)*c(l)
      enddo

      psfcorphot3=phot00
      photerr=sherr*sqrt(photvar)

      end




      function pixresid(sh,x,y,pix,mpix,nx,ny,bg,qq2)
c     gauss-weighted residual flux in image pix after subtraction of
c     shape sh centered on position x. Variance of gaussian wt fn is qq2.
      parameter (rbeta=6.)
      dimension pix(mpix,mpix),sh(msh),xtab(mm+1),ytab(mm+1)
      double precision zz,sum
      beta=sh(1)
      n=sh(2)
      rmax=max(10.,rbeta*max(beta,sqrt(qq2)))
      ix1=max(1,int(x-rmax))
      ix2=min(nx,int(x+rmax+1))
      iy1=max(1,int(y-rmax))
      iy2=min(ny,int(y+rmax+1))
      sum=0
      do ix=ix1,ix2
          xnorm=(ix-x)/beta
c excise as background differences can still be important :
c          if (abs(xnorm).gt.rbeta) goto 11
          call htab(xnorm,n,xtab)
          ex=exp(-0.5*xnorm**2)
          do iy=iy1,iy2
              ynorm=(iy-y)/beta
c excise as background differences can still be important :
c              if (xnorm*xnorm+ynorm*ynorm.gt.rbeta*rbeta) goto 10
              call htab(ynorm,n,ytab)
              ey=exp(-0.5*ynorm**2)
              fac=ex*ey/sqrt(3.1415926535)/beta
              k=2
              zz=0
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
c              write(10,'(2i6,6f15.5)') ix,iy,pix(ix,iy),zz,zz+bg,exp(-0.5*((ix-x)**2+(iy-y)**2)/qq2),sqrt(qq2),(pix(ix,iy)-zz-bg)*exp(-0.5*((ix-x)**2+(iy-y)**2)/qq2)
              zz=pix(ix,iy)-zz-bg
              sum=sum+zz*exp(-0.5*((ix-x)**2+(iy-y)**2)/qq2)
   10         continue
          enddo
   11     continue
      enddo
c      write(0,*) sum,rmax,beta
c      stop
      pixresid=sum
      return
      end


      function apflux(x,y,r,bg,nx,ny,z)
      dimension z(nx,ny)
      double precision sum
      sum=0
      do ix=max(1,int(x-r)),min(nx,int(x+r)+1)
          do iy=max(1,int(y-r)),min(ny,int(y+r)+1)
              r2=(ix-x)**2+(iy-y)**2
              if (r2.le.r*r) sum=sum+z(ix,iy)-bg
          enddo
      enddo
      apflux=sum
      return
      end
