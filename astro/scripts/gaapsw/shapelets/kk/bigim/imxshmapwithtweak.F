      real kercoef(mmpoly2,2*msh-4)
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
c VERSION FOR DOUBLE SHAPELET MAP INCLUDING OUTER TWEAK SHAPELET
c     this version reads and writes entire strips of with m, to curtail disk I/O
c     kernel is a double shapelet, obtained with psfcat2gauskerwithtweak and fitkermaptwk
c     which allows better correction for flux in outer regions of the PSF (typically 5% or so).
      
c read double-shapelet map
 1    read(*,*,err=1) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      nin=(n+1)*(n+2)/2
      do k=1,nin
         read(*,*) (kercoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo
 2    read(*,*,err=2) nout,betaout
      read(*,*) nfitout,xmaxout,ymaxout
      read(*,*)
      do k=1,(nout+1)*(nout+2)/2
         read(*,*) (kercoef(kfit,k+nin),kfit=1,(nfitout+1)*(nfitout+2)/2)
         read(*,*)
      enddo

      write(0,*) 'Read kernel shapelet map'
      write(0,*) 'Order for shape, spatial variation; beta:',n,nfit,beta
      write(0,*) 'These parameters for the outer shapelet: ',nout,nfitout,betaout

c     check whether maps have same spatial order. Would need to change
c     the code to cope with this. (outer loop in shmaptwkft cannot be
c     over all spatial components, would need to have 2 loops, one for
c     inner and one for outer shmap; or pad the kercoefs with zeroes).
      
      if (nfit .ne. nfitout) then
         write(0,*) 'Inner and outer spatial orders differ.'
         stop
      endif

c set width of taper edge
      ntap=nint(4*betaout)

c calculate width of anti-alias border
      nmask=nint(beta*(5.+sqrt(1.*n)))+ntap
      if (nmask.gt.m/4) then
         write(0,*) 'Shapelet map too fat for Fourier convolution.'
         stop
      endif

c     prepare the filters in Fourier domain
      call shmaptwkft(kercoef,beta,n,betaout,nout,nfit,filt,m)

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
c     restore the clipped outliers by adding in the kernel by hand
            if (nclip.gt.0) then
                call addkertwk(out,mx,my,ix1,iy1,nclip,zclip,iclip,jclip,n,beta,nout,betaout,nfit,kercoef,xmax,ymax)
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

      subroutine addkertwk(zz,mx,my,ix1,iy1,nclip,zclip,iclip,jclip,n,beta,nout,betaout,nfit,kercoef,xmax,ymax)
      dimension zz(mx,my),zclip(nclip),iclip(nclip),jclip(nclip)
      real kercoef(mmpoly2,2*msh-4),shker(msh),shkerout(msh)
      nin=(n+1)*(n+2)/2
      do ic=1,nclip
          xc=iclip(ic)
          yc=jclip(ic)
          zc=zclip(ic)
          x=ix1+xc-1
          y=iy1+yc-1
          call interpolpsf(x,y,xmax,ymax,kercoef,n,beta,nfit,shker)
          call interpolpsf(x,y,xmax,ymax,kercoef(1,1+nin),nout,betaout,nfit,shkerout)
          r=(5+sqrt(1.*n))*beta
          i1=max(1,nint(xc-r))
          i2=min(mx,nint(xc+r))
          j1=max(1,nint(yc-r))
          j2=min(my,nint(yc+r))
          do j=j1,j2
              do i=i1,i2
                  if ((i-xc)**2+(j-yc)**2 .le. r*r) then
                      zz(i,j)=zz(i,j)+zc*(shapeletpix(i-xc,j-yc,shker)+shapeletpix(i-xc,j-yc,shker2))
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


      subroutine shmaptwkft(shmap,beta,n,betaout,nout,nfit,filt,npix)
      parameter(twopi=2*3.1415926535)
      complex filt(npix/2,npix,mmpoly2),tfilt
      real xtab(n+1),ytab(n+1),xtabout(nout+1),ytabout(nout+1),shmap(mmpoly2,2*msh-4)

c     construct the FT of a double shapelet map
c     one filter per spatial polynomial term
c     output is a complex array that can be multiplied with output of
c     rlft3 FT of an m x m array   

      nin=(n+1)*(n+2)/2
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
               fynormout=(mod(iy-1+npix/2,npix)-npix/2)*twopi*betaout/npix
               call htab(fynorm,n,ytab)
               call htab(fynormout,nout,ytabout)
               do ix=1,npix/2
                  fxnorm=(mod(ix-1+npix/2,npix)-npix/2)*twopi*beta/npix
                  fxnormout=(mod(ix-1+npix/2,npix)-npix/2)*twopi*betaout/npix
                  call htab(fxnorm,n,xtab)
                  call htab(fxnormout,nout,xtabout)
                  fnorm2=fxnorm*fxnorm+fynorm*fynorm
                  fnorm2out=fxnormout*fxnormout+fynormout*fynormout
                  tfilt=0
                  if (fnorm2.gt.fmax2) goto 6 
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
 6                if (fnorm2out.gt.fmax2) goto 7 
                  fac=exp(-0.5*fnorm2out)/sqrt(3.1415926535)*betaout*twopi/(npix*npix/2)
                  k=0
                  do nn=0,nout
                     do j=0,nn
                        i=nn-j
                        k=k+1
                        if (j.ne.0) then
                           fac=fac/sqrt(j/real(i+1))
                        else 
                           if (i.ne.0) fac=fac/sqrt(real(i))
                        endif
                        tfilt=tfilt+(0.,1.)**nn*shmap(icpt,k+nin)*fac*xtabout(i+1)*ytabout(j+1)
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


      
