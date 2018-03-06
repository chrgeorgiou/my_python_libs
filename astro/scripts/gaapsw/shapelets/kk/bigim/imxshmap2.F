      real kercoef(mmpoly2,msh-2)
      parameter (m=512,maxnx=30000)
      dimension zz(m,m),out(m,m),zbufin(m*maxnx),zbufout(m*maxnx)
      dimension work1(m*m),work2(m*m)
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
            call shmapxpix(filt,m,nfit,zz,out,ix1,iy1,xmax,ymax,work1,work2)
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
     :       ' = pix ',iy1+nmask,' - ',iy1+my-nmask-1
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

      parameter(twopi=2*3.1415926535)
      complex filt(m/2,m,mmpoly2)
      complex speq(m),cpix(m/2,m),cout(m/2,m)
      complex ftim(m/2,m),ctmp(m/2,m)

c Make Fourier transform of pixel array ftim (nondestructively!)
      do iy=1,m
         do ix=1,m/2
            ftim(ix,iy)=cpix(ix,iy)
         enddo
      enddo
      call rlft3(ftim,speq,m,m,1,1)

c filter pixel array with each PSF map polynomial term in turn, add in result
      icpt=0
      do kcpt=0,nfit
         do iypow=0,kcpt
            ixpow=kcpt-iypow
            icpt=icpt+1
            do iy=1,m
               do ix=1,m/2
                  ctmp(ix,iy)=filt(ix,iy,icpt)*ftim(ix,iy)
               enddo
c     set the nyquist freq. component to zero. (if it was not zero there is
c     a problem anyway!)
               speq(iy)=0
            enddo
c     transform back
            call rlft3(ctmp,speq,m,m,1,-1)
c     add result into the output pixel array with proper weight
c     write out full array including padded parts - could trim this
c     write(0,*) 'Adding in part for x,y powers ',ixpow,iypow
            do iy=1,m
               if (iypow.eq.0) then
                  yfac=1
               else
                  yfac=(2*(iy+iy1-1.)/ym-1)**iypow
               endif
               do ix=1,m/2
                  z1=real(ctmp(ix,iy))
                  x1=2*ix-1 + ix1-1
                  z2=aimag(ctmp(ix,iy))
                  x2=2*ix + ix1-1
                  if (ixpow.eq.0) then
                     xfac1=1
                     xfac2=1
                  else
                     xfac1=(2*x1/xm-1)**ixpow
                     xfac2=(2*x2/xm-1)**ixpow
                  endif
                  if (icpt.eq.1) then
                     cout(ix,iy)=yfac*complex(xfac1*z1,xfac2*z2)
                  else
                     cout(ix,iy)=cout(ix,iy)+yfac*complex(xfac1*z1,xfac2*z2)
                  endif
               enddo
            enddo
c     complete the loops over polynomial coefs of the convolution kernel
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


      
