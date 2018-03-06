      real kercoef(mmpoly2,msh-2)
      parameter (m=512,maxnx=30000)
      dimension zz(m,m),out(m,m),zbufin(m*maxnx),zbufout(m*maxnx)
      dimension work1(m*m),work2(m*m)
      dimension zclip(100),iclip(100),jclip(100)
      dimension naxlen(7)
      complex filt(m/2,m,mmpoly2)
      data naxlen /7*0/ 

c as imxshmap, convolves an image by a spatially variable kernel
c expressed as a shapelet map.  differs in that here the image is split
c up into smaller sections and convolved piece by piece. 
c
c this version (2) reads and writes entire strips of with m, to curtail disk I/O

      read(*,*) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(n+1)*(n+2)/2
         read(*,*) (kercoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo

      write(0,*) 'Read kernel shapelet map'
      write(0,*) 'Order for shape, spatial variation; beta:',n,nfit,beta

c set width of taper edge
      ntap=nint(4*beta)

c calculate width of anti-alias border
      nmask=nint(beta*(5.+sqrt(1.*n)))+ntap
      if (nmask.gt.m/4) then
         write(0,*) 'Shapelet map too fat for Fourier convolution.'
         stop
      endif

c prepare the filters in Fourier domain
      call shmapft(kercoef,beta,n,nfit,filt,m)

c open image to be convolved
      call openfits('inimage.fits',51,nx,ny)
      if (nx.gt.maxnx) then
         write(0,*) ' *** IMAGE X dimension exceeds ',maxnx,': too big. *** '
         stop
      endif
      naxlen(1)=nx
      naxlen(2)=ny

c max size that can be fitted into the m x m array that will be FT'd,
c accounting for borders
      mpatch=m-2*nmask

c process the image in sections of mpatch x mpatch, padded to m x m
c outer border of nmask pix wide cannot be used.
      nxsub=int((nx-1)/mpatch)+1
      nysub=int((ny-1)/mpatch)+1
      write(0,*) 'Will process image in ',nxsub,' x ',nysub,' sections.'
      write(0,*) 'Will pad subimages from ',mpatch,mpatch,' to ',m,m

c open output image
      call copyfits_nodata('inimage.fits',61,'convolved.fits')
c get target sig from somewhere, write to header
      betapsf=beta
      call fitshistory(61,'PSF Gaussianized')
      call fitswritekey(61,'GPSFSIG',betapsf,'SIGMA OF GAUSSIANIZED PSF')

c background to subtract before FT processing (and then add back)
      zerr=guess_errors_sec(51,nx*ny)
      bg=get_bg_sec(51,nx*ny)

c go through the input image, read in a patch, mask a border around it,
c convolve, write out
      do iysub=1,nysub
         if (iysub.eq.1) then
            iy1=1-nmask
         else
            iy1=iy1+mpatch
         endif
         iy2=iy1+m-1
         if (iy2.gt.ny+nmask) then
            iy2=ny+nmask
            iy1=max(1,iy2-m+1)
         endif
c get and b/g subtract strip between iy1 and iy2; pad with zeroes to width m if needed
         istat=igetsecpad(51,nx,ny,1,iy1,nx,m,zbufin,nx*m,bg)
c         write(0,*) 'Read lines ',iy1,iy2,' of ',ny
         do ixsub=1,nxsub
            if (ixsub.eq.1) then
               ix1=1-nmask
            else
               ix1=ix1+mpatch
            endif
            ix2=ix1+m-1
            if (ix2.gt.nx+nmask) then
               ix2=nx+nmask
               ix1=max(1,ix2-m+1)
            endif
            mx=ix2-ix1+1
            my=iy2-iy1+1
c            write(0,*) 'Getting image section ',ix1,ix2,iy1,iy2
            call frombuf(nx,m,zbufin,ix1,1,mx,my,zz)
c            istat=igetsecpad(51,nx,ny,ix1,iy1,mx,my,zz,mx*my,bg)
            if (istat.ne.0.) then
               write(0,*) 'igetsecpad failed. code = ',istat
               stop
            endif
c convolve with shmap here
c taper the array to reduce aliasing
            call taperedge(zz,mx,my,ntap)    
c clip to remove extreme 1000sig outliers
            call clipimage(zz,mx,my,zerr,nclip,zclip,iclip,jclip)
            call shmapxpix(filt,m,nfit,zz,out,ix1,iy1,xmax,ymax,work1,work2)
c restore the clipped outliers by adding in the kernel by hand
            if (nclip.gt.0) then
                call addker(out,mx,my,ix1,iy1,nclip,zclip,iclip,jclip,n,beta,nfit,kercoef,xmax,ymax)
                write(0,*) 'Kept ',nclip,'outliers outside FFT.',
     :              ' (in block with columns ',ix1,'-',ix2,')'
            endif
c
c convolution done, result in out
c            write(0,*) 'Writing subsection ',ixsub,iysub,ix1,iy1
c pick out the part of out that is unaffected by aliasing (strip off border)
            call stripborder(out,zz,mx-2*nmask,my-2*nmask,m,nmask)
            call subbg(zz,(mx-2*nmask)*(my-2*nmask),-bg)
c            write(0,*) 'Write to section starting at ',ix1,iy1
            call tobuf(nx,ny,zbufout,ix1+nmask,1,mx-2*nmask,my-2*nmask,zz)
c            call putsec(61,zz,mx-2*nmask,my-2*nmask,ix1+nmask,iy1+nmask,naxlen)
         enddo
c         write(0,*) 'Write out lines ',iy1+nmask,iy1-nmask+my-1,' of ',nx,ny
         call putsec(61,zbufout,nx,my-2*nmask,1,iy1+nmask,naxlen)
         write(0,*) 'Finished strip ',iysub,' / ',nysub,
     :       ' = row ',iy1+nmask,' - ',iy1+my-nmask-1
      enddo

c      write(0,*) 'Convolved image in convolved.fits'
c      end


      call closefits(51)

c get target sig from somewhere
      betapsf=beta
c      open(13,file='ker.in',status='old',err=11)
c      read(13,*,err=11) betapsf
c      close(13)
c   11 continue

      call closefits(61)
      end


      subroutine frombuf(nx,ny,zbuf,ix1,iy1,mx,my,zz)
      dimension zz(mx,my),zbuf(nx,ny)
c copy mx*my array starting at ix1,iy1 from nx*ny array zbuf
c if index starts below 1, pad with zero rows
      do iy=iy1,0
         iyout=iy-iy1+1
         do ix=1,mx
            zz(ix,iy)=0
         enddo
      enddo         
c rows that correspond to lines 1-ny
      do iy=max(1,iy1),min(iy1+my-1,ny)
         iyout=iy-iy1+1
         do ix=ix1,0
            ixout=ix-ix1+1
            zz(ixout,iyout)=0
         enddo
         do ix=max(1,ix1),min(ix1+mx-1,nx)
            ixout=ix-ix1+1
            zz(ixout,iyout)=zbuf(ix,iy)
         enddo
         do ix=ixout+1,mx
            zz(ix,iyout)=0
         enddo
      enddo
c rows that exceed ny
      do iy=iyout+1,my
         do ix=1,mx
            zz(ix,iy)=0
         enddo
      enddo
      return
      end


      subroutine tobuf(nx,ny,zbuf,ix1,iy1,mx,my,zz)
      dimension zz(mx,my),zbuf(nx,ny)
c copy mx*my array zz to nx*ny array zbuf, starting at ix1,iy1. do not copy across the edges.
      do iy=iy1,min(iy1+my-1,ny)
         iyin=iy-iy1+1
         do ix=ix1,min(ix1+mx-1,nx)
            ixin=ix-ix1+1
            zbuf(ix,iy)=zz(ixin,iyin)
         enddo
      enddo
      return
      end


      subroutine shmapxpix(filt,m,nfit,cpix,cout,ix1,iy1,xm,ym,ftim,ctmp)

c convolve pixel array pix with shapelet map shmap
c Use the first npix x npix pixels of the mpix x mpix pixel array pix
c pix can be called as a real npix x xpix array - inside this routine it is
c treated as a complex npix/2 x npix array cpix
c     ftim and ctmp are temporary arrays that need to be dimensioned as
c     real npix*npix for the FT routine

c in this implementation evaluate kernel at the four courners of the
c postage stamp and use bilinear interpolation.
c should be a lot faster than evaluating all polynomial terms of the
c kernel variation and adding those.


      parameter(twopi=2*3.1415926535)
      complex filt(m/2,m,mmpoly2)
      complex speq(m),cpix(m/2,m),cout(m/2,m)
      complex ftim(m/2,m),ctmp(m/2,m)

c Make Fourier transform of pixel array cpix (nondestructively!), result in ftim
      do iy=1,m
         do ix=1,m/2
            ftim(ix,iy)=cpix(ix,iy)
         enddo
      enddo
      call rlft3(ftim,speq,m,m,1,1)

      if (nfit.eq.0) then
c constant kernel - just filter with the one kernel filt(ix,iy,1)
         do iy=1,m
            do ix=1,m/2
               cout(ix,iy)=filt(ix,iy,1)*ftim(ix,iy)
            enddo
            speq(iy)=0
         enddo
         call rlft3(cout,speq,m,m,1,-1)
      else
c variable kernel - approximate with bilinear interpolation between the corners of the postage stamp.
c loop over the four corners (note that linear case could by done slightly more efficiently)
         do iy=1,m
            do ix=1,m/2
               cout(ix,iy)=0.
            enddo
         enddo
         do jxcorner=0,1
            do jycorner=0,1
               ixcorner=ix1+jxcorner*(m-1)
               iycorner=iy1+jycorner*(m-1)
               do iy=1,m
                  do ix=1,m/2
                     ctmp(ix,iy)=0
                  enddo
                  speq(iy)=0
               enddo
c build FT of kernel at this corner from FT of polynomial components in filt(ix,iy,icpt) array
               icpt=0
               do kcpt=0,nfit
                  do iypow=0,kcpt
                     ixpow=kcpt-iypow
                     icpt=icpt+1
                     if (ixpow.eq.0) then
                        xfac=1
                     else
                        xfac=(2*(ixcorner-1.)/xm-1)**ixpow
                     endif
                     if (iypow.eq.0) then
                        yfac=1
                     else
                        yfac=(2*(iycorner-1.)/ym-1)**iypow
                     endif
                     do iy=1,m
                        do ix=1,m/2
                           ctmp(ix,iy)=ctmp(ix,iy)+filt(ix,iy,icpt)*xfac*yfac
                        enddo
                        speq(iy)=0
                     enddo
                  enddo
               enddo
c convolve pixel array with filter
               do iy=1,m
                  do ix=1,m/2
                     ctmp(ix,iy)=ctmp(ix,iy)*ftim(ix,iy)
                  enddo
               enddo
               call rlft3(ctmp,speq,m,m,1,-1)
c add result into output using bilinear interp weight.
c t and u run from 0 to 1 in x and y over the patch.
c in x the array is stored as complex numbers, so need to pull out pair of pixels as real and imag pt
               do iy=1,m
                  u=(iy-1.)/(m-1.)
                  if (jycorner.eq.0.) u=1.-u
                  do ix=1,m/2
                     x1=2*ix-1
                     x2=x1+1
                     t1=(2*ix-2.)/(m-1.)
                     t2=(2*ix-1.)/(m-1.)
                     if (jxcorner.eq.0) then
                        t1=1.-t1
                        t2=1.-t2
                     endif
                     z1=real(ctmp(ix,iy))
                     z2=aimag(ctmp(ix,iy))
                     cout(ix,iy)=cout(ix,iy)+u*complex(t1*z1,t2*z2)
                  enddo
               enddo
c     return to add in remaining corners
            enddo
         enddo
      endif
      return
      end


      subroutine clipimage(zz,mx,my,zerr,nclip,zclip,iclip,jclip)
      parameter (nclipmax=50)
      dimension zz(mx,my),zclip(nclipmax),iclip(nclipmax),jclip(nclipmax)
c find the up to 50 worst 1E6-sigma outliers in z and replace them with
c zero before the Fourier analysis.
      nclip=0
      zclip(1)=0.
      thresh=1e6*zerr
      do j=1,my
          do i=1,mx
              if (abs(zz(i,j)).gt.thresh) then
                  nclip=min(nclipmax,nclip+1)
                  zc=zz(i,j)
                  do ic=1,nclip
                      if (abs(zc).gt.abs(zclip(ic))) then
                          do k=nclipmax,ic+1,-1
                              zclip(k)=zclip(k-1)
                              iclip(k)=iclip(k-1)
                              jclip(k)=jclip(k-1)
                          enddo
                          zclip(ic)=zc
                          iclip(ic)=i
                          jclip(ic)=j
                          goto 5
                      endif
                  enddo
    5         endif
          enddo
      enddo
      do ic=1,nclip
          zz(iclip(ic),jclip(ic))=0.
      enddo
      return
      end

      subroutine addker(zz,mx,my,ix1,iy1,nclip,zclip,iclip,jclip,n,beta,nfit,kercoef,xmax,ymax)
      dimension zz(mx,my),zclip(nclip),iclip(nclip),jclip(nclip)
      real kercoef(mmpoly2,msh-2),shker(msh)
      do ic=1,nclip
          xc=iclip(ic)
          yc=jclip(ic)
          zc=zclip(ic)
          x=ix1+xc-1
          y=iy1+yc-1
          call interpolpsf(x,y,xmax,ymax,kercoef,n,beta,nfit,shker)
          r=(5+sqrt(1.*n))*beta
          i1=max(1,nint(xc-r))
          i2=min(mx,nint(xc+r))
          j1=max(1,nint(yc-r))
          j2=min(my,nint(yc+r))
          do j=j1,j2
              do i=i1,i2
                  if ((i-xc)**2+(j-yc)**2 .le. r*r) then
                      zz(i,j)=zz(i,j)+zc*shapeletpix(i-xc,j-yc,shker)
                  endif
              enddo
          enddo
      enddo
      return
      end


      subroutine subbg(zz,m,bg)
c subtract the constant value bg from the zz array
      dimension zz(m)
      do i=1,m
          zz(i)=zz(i)-bg
      enddo
      return
      end

      subroutine taperedge(zz,mx,my,ntap)
c taper a strip of width ntap down to zero smoothly
      dimension zz(mx,my)
      do iy=1,ntap
         filt=(iy-1.)/ntap
         filt=(1-cos(3.1415926535*filt))/2.
         do ix=1,mx
            zz(ix,iy)=zz(ix,iy)*filt
         enddo
         do ix=1,mx
            zz(ix,my-iy+1)=zz(ix,my-iy+1)*filt
         enddo
      enddo
      do ix=1,ntap
         filt=(ix-1.)/ntap
         filt=(1-cos(3.1415926535*filt))/2.
         do iy=1,my
            zz(ix,iy)=zz(ix,iy)*filt
         enddo
         do iy=1,my
            zz(mx-ix+1,iy)=zz(mx-ix+1,iy)*filt
         enddo
      enddo
      return
      end

      subroutine stripborder(full,stripped,mx,my,m,nmask)
c write the mx,my sized block of the full array starting at
c nmask+1,nmask+1 to the stripped array
      dimension stripped(mx,my),full(m,m)
      do j=nmask+1,min(m,my+nmask)
         do i=nmask+1,min(m,mx+nmask)
            stripped(i-nmask,j-nmask)=full(i,j)
         enddo
      enddo
      return
      end


      subroutine shmapft(shmap,beta,n,nfit,filt,npix)
      parameter(twopi=2*3.1415926535)
      complex filt(npix/2,npix,mmpoly2),tfilt
      real xtab(n+1),ytab(n+1),shmap(mmpoly2,msh-2)

c     construct the FT of a shapelet map
c     one filter per spatial polynomial term
c     output is a complex array that can be multiplied with output of
c     rlft3 FT of an m x m array   

      icpt=0
      do kcpt=0,nfit
         do ixpow=kcpt,0,-1
            iypow=kcpt-ixpow
            icpt=icpt+1
c     construct shapelet function's FT at frequency fxnorm, fynorm
c     Hn(x/beta)exp(-0.5*x**2/beta**2) --> i^n Hn(f*beta) exp(-0.5*f**2
c     *beta**2)
            fmax2=(sqrt(n*1.)+3)**2
            do iy=1,npix
               fynorm=(mod(iy-1+npix/2,npix)-npix/2)*twopi*beta/npix
               call htab(fynorm,n,ytab)
               do ix=1,npix/2
                  fxnorm=(mod(ix-1+npix/2,npix)-npix/2)*twopi*beta/npix
                  call htab(fxnorm,n,xtab)
                  fnorm2=fxnorm*fxnorm+fynorm*fynorm
                  tfilt=0
                  if (fnorm2.gt.fmax2) goto 7 
                  fac=exp(-0.5*fnorm2)/sqrt(3.1415926535)*beta*twopi/(npix*npix/2)
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
                        tfilt=tfilt+(0.,1.)**nn*shmap(icpt,k)*fac*xtab(i+1)*ytab(j+1)
                     enddo
                     fac=fac/sqrt(2.)
                  enddo
c     filter FT of pixel array
 7                filt(ix,iy,icpt)=tfilt
               enddo
            enddo
c            write(0,*) 'Made FT of PSF map poly cpt ',icpt,ixpow,iypow
         enddo
      enddo
      write(0,*) 'Made FT''s of PSF map poly cpts 1 to',icpt
      return
      end


      

c     fit the ndat data points (x,y,tofit) with polynomials in x and y
c     up to total order norder. Weights w are given as input; 
c     for clipped data w is set to zero at output. 
c     xm and ym are the factors by which to rescale x and y, 
c     niter the number of sigma-clipping iterations (at clip sigmas).  
c     output is the list of coefficients coefs, in order
c     00, 10, 01, 20, 11, 02, ...
c     and a list of best-fit values for the data points.


      subroutine fitpoly(ndat,xm,ym,x,y,tofit,w,norder,niter,clip,coefs,bestfit)
      real x(ndat),y(ndat),tofit(ndat),bestfit(ndat),w(ndat)
      real coefs(mmpoly2),mat(mmpoly2,mmpoly2),model(mmpoly2)
      double precision matd(mmpoly2,mmpoly2),coefsd(mmpoly2)
      real resid(ndat)
      logical more
      integer indx(mmpoly2)
      double precision sum

      nmod=(norder+1)*(norder+2)/2

      do iter=1,niter
          do kk=1,nmod
              do ll=1,nmod
                  matd(kk,ll)=0
              enddo
              coefsd(kk)=0
          enddo
          do i=1,ndat
              kk=0
              do k=0,norder
                  do m=0,k
                      l=k-m
                      kk=kk+1
                      model(kk)=polymodel(x(i),y(i),xm,ym,l,m)
                  enddo
              enddo
              do kk=1,nmod
                  do ll=1,nmod
                      matd(kk,ll)=matd(kk,ll)+model(kk)*model(ll)*w(i)
                  enddo
                  coefsd(kk)=coefsd(kk)+model(kk)*tofit(i)*w(i)
              enddo
          enddo
          do kk=1,nmod
              do ll=1,nmod
                  mat(kk,ll)=matd(kk,ll)
              enddo
              coefs(kk)=coefsd(kk)
          enddo
          call ludcmp(mat,nmod,mmpoly2,indx,det)
          call lubksb(mat,nmod,mmpoly2,indx,coefs)
          more=.false.
          do i=1,ndat
              kk=0
              do k=0,norder
                  do m=0,k
                      l=k-m
                      kk=kk+1
                      model(kk)=polymodel(x(i),y(i),xm,ym,l,m)
                  enddo
              enddo
              sum=0
              do kk=1,nmod
                  sum=sum+model(kk)*coefs(kk)
              enddo
              bestfit(i)=sum
              resid(i)=tofit(i)-bestfit(i)
          enddo
c calculate median and mean abs dev of the residuals
          residmed=fmedian(ndat,resid)
          sum0=0
          sum1=0
          do i=1,ndat
              sum0=sum0+w(i)
              sum1=sum1+abs(resid(i)-residmed)*w(i)
          enddo
          avgdev=sum1/sum0
          do i=1,ndat
              if (abs(resid(i)-residmed).gt.clip*avgdev .and. w(i).gt.0.) then 
                  w(i)=0
                  more=.true.
              endif
          enddo
          if (.not.more) goto 20
      enddo
   20 continue
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



      SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign)
      INTEGER isign,nn1,nn2,nn3
      COMPLEX data(nn1/2,nn2,nn3),speq(nn2,nn3)
CU    USES fourn
      INTEGER i1,i2,i3,j1,j2,j3,nn(3)
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      COMPLEX c1,c2,h1,h2,w
      c1=cmplx(0.5,0.0)
      c2=cmplx(0.0,-0.5*isign)
      theta=6.28318530717959d0/dble(isign*nn1)
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      nn(1)=nn1/2
      nn(2)=nn2
      nn(3)=nn3
      if(isign.eq.1)then
        call fourn(data,nn,3,isign)
        do 12 i3=1,nn3
          do 11 i2=1,nn2
            speq(i2,i3)=data(1,i2,i3)
11        continue
12      continue
      endif
      do 15 i3=1,nn3
        j3=1
        if (i3.ne.1) j3=nn3-i3+2
        wr=1.0d0
        wi=0.0d0
        do 14 i1=1,nn1/4+1
          j1=nn1/2-i1+2
          do 13 i2=1,nn2
            j2=1
            if (i2.ne.1) j2=nn2-i2+2
            if(i1.eq.1)then
              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
              data(1,i2,i3)=h1+h2
              speq(j2,j3)=conjg(h1-h2)
            else
              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
              data(i1,i2,i3)=h1+w*h2
              data(j1,j2,j3)=conjg(h1-w*h2)
            endif
13        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          w=cmplx(sngl(wr),sngl(wi))
14      continue
15    continue
      if(isign.eq.-1)then
        call fourn(data,nn,3,isign)
      endif
      return
      END


      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      REAL data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
     *k2,n,nprev,nrem,ntot
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END



      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END


      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END



      function fmedian(n,x)
      real x(n),xx(10*mmcat)
      do i=1,n
          xx(i)=x(i)
      enddo
      call sort(n,xx)
      if ((n/2)*2.eq.n) then
          fmedian=( xx(n/2)+xx(n/2+1) ) /2
      else
          fmedian=xx(n/2+1)
      endif
      return
      end


      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END


      subroutine openfits(imagename,nunit,nx,ny)
      character*(*) imagename
      dimension naxlen(7)
      irw=0
      iblocksize=0
      istat=0
      call ftnopn(nunit,imagename,irw,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Open failed: status=',istat,', ',imagename
         stop
      endif
      call ftgidm(nunit,naxis,istat)
      call ftrprt('stderr',istat)
      call ftgisz(nunit,naxis,naxlen,istat)
      call ftrprt('stderr',istat)
      if (naxis.ne.2) then
         write(0,*) 'Not a 2D image: ',imagename
         stop
      endif
      nx=naxlen(1)
      ny=naxlen(2)
      write(0,*) 'Opened ',nx,' by ',ny,' image ',imagename
      return
      end

      subroutine closefits(nunit)
      istat=0
      call ftclos(nunit,istat)
      call ftrprt('stderr',istat)
      return
      end
      

      subroutine copyfits_nodata(image1,nunit2,image2)
      character image2*(*),image1*(*)
c copy fits file (file image1) to a new file, image2 and open as unit 2.
c Headers of image1, not its data, are copied.
c delete existing image2 if necessary
      open(unit=nunit2,file=image2,status='old',err=1)
      close(nunit2,status='delete')
      write(0,*) 'Deleted existing file ',image2
    1 continue
      istat=0
      call fttplt(nunit2,image2,image1,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Creating new image failed: status=',istat
         write(0,*) image2
         stop
      endif
      return
      end

      subroutine fitshistory(nunit,comment)
      character comment*(*)
      istat=0
      call ftphis(nunit,comment,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Write history failed: status=',istat
         write(0,*) 'image unit ',nunit
         stop
      endif
      return
      end


      subroutine fitswritekey(nunit,keyword,value,comment)
      character*(*) keyword,comment
      istat=0
      call ftpkye(nunit,keyword,value,-6,comment,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'FITS: Write keyword failed: status=',istat
         write(0,*) 'image unit ',nunit
         write(0,*) 'keyword ',keyword
         write(0,*) 'Value ',value
         write(0,*) 'Comment ',comment
         stop
      endif
      return
      end


      function guess_errors_sec(nunit,ntot)
      parameter (nprime=10099)
      real sky(nprime)
      logical anyf
c take prime number to avoid many samples form same col or row.

c sample image at ~ 100x100 points, look at the median and lower quartile
c assign a sigma based on that

c subsample from the ntot pixel 1-d data stream
      k=0
      istat=0
      if (ntot.le.nprime) then
         i1=1
         istep=1
         nstep=ntot
      else
         istep=ntot/nprime 
         i1=(ntot-istep*(nprime-1))/2 + 1
         nstep=nprime
      endif
      do i=i1,i1+(nstep-1)*istep,istep
         k=k+1
         call ftgpve(nunit,0,i,1,0.,sky(k),anyf,istat)
c do not include values exactly equal to zero
c (should do something better here with masks/weights...)
         if (sky(k).eq.0.) k=k-1
c         call ftrprt('stderr',istat)
c         write(0,*) i,k,sky(k)
      enddo
c find median and lower quartile
      call sort(k,sky)
      call get_med_quartile(k,sky,bg,qdev)
c for a Gaussian, quartile is at 0.675 sigma
      sig=qdev/0.675
      write(0,2000) k,bg,sig
 2000 format (' B/g, noise estimate from sample of ',i5,' pixels:',2g12.4)
c use first estimate to remove values very close to the median (empty
c regions of the image)
      zaplo=bg-0.001*sig
      zaphi=bg+0.001*sig
      ilo=1
c make ilo first item above lower limit
      do while (sky(ilo).lt.zaplo)
          ilo=ilo+1
      enddo
      nkill=0
c make ihi last item below lower limit
      do while (sky(ilo+nkill).lt.zaphi)
          nkill=nkill+1
      enddo
c excise nkill items from the list, starting at ilo
      if (nkill.gt.0) then
          do i=ilo,k-nkill
              sky(i)=sky(i+nkill)
          enddo
      endif
c repeat median and quartile determination:
      k=k-nkill
      call get_med_quartile(k,sky,bg,qdev)
      sig=qdev/0.675
      write(0,2000) k,bg,sig

      guess_errors_sec=sig

      open(27,file='bgnoise.dat',status='unknown')
      write(27,*) bg,sig
      close(27)
c      call pgbeg(0,'noise-est.ps/ps',1,1)
c      plotmin=sky(nint(0.05*k)+1)
c      plotmax=sky(nint(0.95*k))
c      plotwidth=plotmax-plotmin
c      plotmin=plotmin-0.1*plotwidth
c      plotmax=plotmax+0.1*plotwidth
c      call pghist(k, sky, plotmin,plotmax, 100, 2)
c      call pgsls(3)
c      call pgmove(bg,0.)
c      call pgdraw(bg,real(k))
c      call pgmove(bg-sig,0.)
c      call pgdraw(bg-sig,real(k))
c      call pgsls(1)
c      call pglab('PIXEL VALUES','#','BG AND NOISE ESIMATE')
c      call pgend
      return
      end


      function guess_bg_sec(nunit,ntot)
      parameter (nprime=10099)
      real sky(nprime)
      logical anyf
c estimate median value of an image by subsampling 1D pixel stream
c take prime number to avoid many samples form same col or row.

c sample image at ~ 100x100 points, look at the median and lower quartile
c assign a bg based on that

c subsample from the ntot pixel 1-d data stream
      k=0
      istat=0
      if (ntot.le.nprime) then
         i1=1
         istep=1
         nstep=ntot
      else
         istep=ntot/nprime 
         i1=(ntot-istep*(nprime-1))/2 + 1
         nstep=nprime
      endif
      do i=i1,i1+(nstep-1)*istep,istep
         k=k+1
         call ftgpve(nunit,0,i,1,0.,sky(k),anyf,istat)
c do not include values exactly equal to zero
         if (sky(k).eq.0.) k=k-1
      enddo

      guess_bg_sec=fmedian(k,sky)
      return
      end

      function igetsecpad(nunit,nxtot,nytot,ix1,iy1,nx,ny,zz,maxx,bg)

c     from previously opened fits file with unit nunit, size nxtot x nytot,
c     extract nx by ny pixel rectangle starting at ix1,iy1.
c     Background subtract all pixels with value bg
c     Result gets put in 1D zz array (which can be
c     declared as a 2d nx by ny array elsewhere).
c     
c     differs from igetsec: if the extracted rectangle overflows the
c     border of the fits file the patch is padded with zeroes, rather
c     than generating an error
c
c     igetsec takes a diagnostic value on return:
c     0 if all OK
c     -999 if nx*ny too large for compiled size of zz array
c     istat code from cfitsio otherwise.

      dimension zz(maxx),naxlen(7)
      integer fpix(7),lpix(7),incpix(7)
      data fpix,lpix,incpix /21*1/
      data naxlen /7*0/
      logical anyf

      istat=0
      igetsecpad=0
      naxlen(1)=nxtot
      naxlen(2)=nytot
      if (nx*ny.gt.maxx) then
         write(0,*) 'iGETSECPAD: patch dimensions',nx,ny,' too large.'
         igetsecpad=-999
         return
      endif
c extract the part of the rectangle that fits inside the fits image
c write it in straight fortran order, no gaps
      fpix(1)=max(1,ix1)
      fpix(2)=max(1,iy1)
      lpix(1)=min(ix1+nx-1,naxlen(1))
      lpix(2)=min(iy1+ny-1,naxlen(2))
ccc      write(0,*) 'Getting section',fpix(1),fpix(2),lpix(1),lpix(2)
      call ftgsve(nunit,0,2,naxlen,fpix,lpix,incpix,0,zz,anyf,istat)
      call ftrprt('stderr',istat)
      igetsecpad=istat
c check whether the rectangle needs to be padded
      if (ix1.lt.1 .or. ix1+nx-1.gt.naxlen(1) .or.
     :     iy1.lt.1 .or. iy1+ny-1.gt.naxlen(2)) then
ccc         write(0,*) 'iGETSECPAD: Image section exceeds boundary',ix1,nx,nxtot,iy1,ny,nytot
ccc         write(0,*) 'Padding with zeros to size ',nx,' x ',ny
c last line/col to zero-pad at beginning
         ixpad1=max(0,1-ix1)
         iypad1=max(0,1-iy1)
c first line/col to zero-pad at the end
         ixpad2=min(nx+1,naxlen(1)-ix1+2)
         iypad2=min(ny+1,naxlen(2)-iy1+2)
ccc         write(0,*) 'Padding borders in x,y:',ixpad1,ixpad2,iypad1,iypad2
c now copy the zz array from high to low index, inserting zeros as we go.
c also do the background subtraction here
c keep two running indices with location where pixels should be written and where they come from
c izin should always be <= izout, so can work from high to low
         izout=nx*ny
         izin=(lpix(1)-fpix(1)+1)*(lpix(2)-fpix(2)+1)
ccc         write(0,*) 'start at izout,izin=',izout,izin
c top lines that need to be zeroed out entirely
         do iy=ny,iypad2,-1
            do ix=nx,1,-1
               zz(izout)=0
               izout=izout-1
            enddo
         enddo
c middle lines that may need to be zero-padded in x
         do iy=iypad2-1,iypad1+1,-1
            do ix=nx,ixpad2,-1
               zz(izout)=0
               izout=izout-1
            enddo
            do ix=ixpad2-1,ixpad1+1,-1
               zz(izout)=zz(izin)-bg
               izout=izout-1
               izin=izin-1
            enddo
            do ix=ixpad1,1,-1
               zz(izout)=0
               izout=izout-1
            enddo
         enddo
c bottom lines that need to be zeroed out entirely
         do iy=iypad1,1,-1
            do ix=nx,1,-1
               zz(izout)=0
               izout=izout-1
            enddo
         enddo
ccc        write(0,*) 'izout, izin:',izout,izin
         if (izout.ne.0 .or. izin.ne.0) then
            write(0,*) 'iGETSECPAD: Something went wrong, ',
     :           'not all pixels correctly copied'
            igetsecpad=-1
         endif
      else
c no need to pad; just du background subtraction
ccc         write(0,*) 'Background subtracting: nx,ny,bg=',nx,ny,bg
         do i=1,nx*ny
            zz(i)=zz(i)-bg
         enddo
      endif
c      write(0,*) fpix,lpix,incpix
c      write(0,*) 'Read',nx,' x ',ny,' image section.'
      return
      end


      function get_bg_sec(nunit,ntot)
      parameter (nprime=10099)
      real sky(nprime)
      logical anyf
c take prime number to avoid many samples form same col or row.

c sample image at 100x100 points, look at the median and lower quartile
c assign a uniform sigma based on that

c subsample from the ntot pixel 1-d data stream
      k=0
      istat=0
      if (ntot.le.nprime) then
         i1=1
         istep=1
         nstep=ntot
      else
         istep=ntot/nprime 
         i1=(ntot-istep*(nprime-1))/2 + 1
         nstep=nprime
      endif
      do i=i1,i1+(nstep-1)*istep,istep
         k=k+1
         call ftgpve(nunit,0,i,1,0.,sky(k),anyf,istat)
c do not include values exactly equal to zero
c (should do something better here with masks/weights...)
         if (sky(k).eq.0.) k=k-1
c         call ftrprt('stderr',istat)
c         write(0,*) i,k,sky(k)
      enddo
      call sort(k,sky)
      bg=0.5*(sky(k/2+1)+sky((k+1)/2))
      write(0,*) 'Estimated bg level from sample of ',k,
     :    ' pixels: ',bg
      get_bg_sec=bg
      return
      end


      subroutine get_med_quartile(n,y,ym,qdev)
c from SORTED array y(n), iteratively find median and lower quartile width, in a way that is robust to excess of positive values
c first find true median M and lower quartile Q1
c then find median of all values between Q1 and M + (M-Q1)
c then recalculate Q1 as median of all values below M
c return M and qdev=M-Q1
      dimension y(n)
      ym=(y(n/2+1) + y((n+1)/2)) /2
      iq1=nint(0.25*n+0.75)
      q1=y(iq1)
      q2=2*ym-q1
      iq2=iq1+1
      do while (y(iq2).lt.q2 .and. iq2.lt.n)
         iq2=iq2+1
      enddo
      ym=0.5*(y((iq1+iq2)/2+1)+y((iq1+iq2+1)/2))
      iq1=nint(0.25*(iq1+iq2+0.5))
      q1=y(iq1)
      qdev=ym-q1
      return
      end


      subroutine putsec(nunit,pix,nx,ny,ix1,iy1,naxlen)
c write nx,ny array pix into open fits image with unit nunit.
c bottom left pixel goes to ix1,iy1
c array pix must have dimension nx,ny ! No gaps or borders
      dimension pix(nx,ny)
      dimension naxlen(7),ifpix(7),ilpix(7)
      data ifpix,ilpix /14*0/

c write the data.
c      write(0,*) 'putsec: ',nunit,naxlen
c      write(0,*) 'putsec: ',pix(1,1),pix(nx,ny)
      ndim=2
      igroup=1
      ifpix(1)=ix1
      ilpix(1)=ix1+nx-1
      ifpix(2)=iy1
      ilpix(2)=iy1+ny-1
      istat=0
      call ftpsse(nunit,igroup,ndim,naxlen,ifpix,ilpix,pix,istat)
      if (istat.ne.0) then
         call ftrprt('stderr',istat)
         write(0,*) 'Putsec: Write data failed: status=',istat
         stop
      endif
      end

