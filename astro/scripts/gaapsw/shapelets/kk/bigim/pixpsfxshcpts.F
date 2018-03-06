c convolve the PSF polynomial terms with sheared circular shapelets,
c save as fits file.

      parameter (mstar=15,mstar2=2*mstar+1,maxbeta=25,maxround=3)
      parameter (mft=64)
      dimension psfcpts(mstar2,mstar2,mmpoly2)
      complex ftpsfcpts(mft/2,mft,mmpoly2),speq(mft*mft/2)
      complex filt(mft/2,mft,9)
      dimension out(mstar2,mstar2,mmpoly2*(3+maxbeta*maxround*3*3))
      dimension sh(msh),sh2(msh),sh3(msh),sh4(msh)

      call get3('pixpsf.fits',nx,ny,ncpt,psfcpts,mstar2,mstar2,mmpoly2)
      write(*,*) 'ok'
      if (nx.ne.mstar2 .or. ny.ne.mstar2) then
         write(0,*) 'pixpsf.fits does not have x,y dimensions ',mstar2,
     :        ' square!'
         stop
      endif

c compute FT's of padded PSF components
      do icpt=1,ncpt
         call copypadsq(psfcpts(1,1,icpt),ftpsfcpts(1,1,icpt),
     :        mstar2,mft)
         do i=1,mft
            speq(i)=0
         enddo
         call rlft3(ftpsfcpts(1,1,icpt),speq,mft,mft,1,1)
      enddo
      write(0,*) 'Fourier tranformed PSF cpts'

c now start to build the model components.

c first components are the raw PSF
      do icpt=1,ncpt
         call copyslice(psfcpts(1,1,icpt),out(1,1,icpt),mstar2)
      enddo
      iout=ncpt
c     INCLUDE SHIFTED PSF COMPONENTS AS WELL.
c     use shiftpatch(pix,nx,ny,dx,dy) with step of 0.1 pix.

      do icpt=1,ncpt
         iout=iout+1
         call copyslice(psfcpts(1,1,icpt),out(1,1,iout),mstar2)
         call shiftpatch(out(1,1,iout),mstar2,mstar2,-0.05,0.)
         do j=1,mstar2
            do i=1,mstar2
               out(i,j,iout)=20.*(out(i,j,iout)-psfcpts(i,j,icpt))
            enddo
         enddo
      enddo
      do icpt=1,ncpt
         iout=iout+1
         call copyslice(psfcpts(1,1,icpt),out(1,1,iout),mstar2)
         call shiftpatch(out(1,1,iout),mstar2,mstar2,0.,-0.05)
         do j=1,mstar2
            do i=1,mstar2
               out(i,j,iout)=20.*(out(i,j,iout)-psfcpts(i,j,icpt))
            enddo
         enddo
      enddo
         
      
      beta=1.
c for harmonic range of Gaussian radii, build 3 circular shapelets,
c shear them, and convolve with PSF.
c set shapelet order such that sheared terms not truncated!
c Note that translations are in units of beta!
      do ibeta=0,maxbeta-1
         sh(1)=beta
         do ic=0,2*(maxround-1),2
            do i=2,msh
               sh(i)=0
            enddo
            if (ic.eq.0) then
               sh(2)=3
               sh(3)=1/sqrt(beta)
            elseif (ic.eq.2) then
               sh(2)=5
               sh(6)=-sqrt(0.5/beta)
               sh(8)=sh(6)
            elseif (ic.eq.4) then
               sh(2)=7
               sh(13)=sqrt(0.375/beta)
               sh(15)=0.5/sqrt(beta)
               sh(17)=sh(13)
            endif
            call shfilt(sh,filt(1,1,1),mft)
            call shear1(sh,sh2)
            call shfilt(sh2,filt(1,1,2),mft)
            call shear2(sh,sh3)
            call shfilt(sh3,filt(1,1,3),mft)
            call trans1(sh,sh4)
            call shcmult(sh4,1./beta)
            call shfilt(sh4,filt(1,1,4),mft)
            call trans1(sh2,sh4)
            call shcmult(sh4,1./beta)
            call shfilt(sh4,filt(1,1,5),mft)
            call trans1(sh3,sh4)
            call shcmult(sh4,1./beta)
            call shfilt(sh4,filt(1,1,6),mft)
            call trans2(sh,sh4)
            call shcmult(sh4,1./beta)
            call shfilt(sh4,filt(1,1,7),mft)
            call trans2(sh2,sh4)
            call shcmult(sh4,1./beta)
            call shfilt(sh4,filt(1,1,8),mft)
            call trans2(sh3,sh4)
            call shcmult(sh4,1./beta)
            call shfilt(sh4,filt(1,1,9),mft)

            do ifilt=1,9
               do icpt=1,ncpt
                  iout=iout+1
                  call convolve(ftpsfcpts(1,1,icpt),mft,filt(1,1,ifilt),
     :                 out(1,1,iout),mstar2)
               enddo
            enddo
         enddo
         write(0,*) 'Completed calculation for beta=',beta,
     :        '  --  mod count ',iout
         beta=beta*2.**0.125
      enddo

      call put3('pixpsfshearcpts.fits',mstar2,mstar2,iout,out,mstar2,mstar2,iout-1)

      end


      subroutine copyslice(in,out,mstar2)
      real in(mstar2*mstar2),out(mstar2*mstar2)
c copy 2D array 'in' to 2D array 'out'. No padding. More efficient than copypadsq.
      do ipix=1,mstar2*mstar2
         out(ipix)=in(ipix)
      enddo
      return
      end


      subroutine copypadsq(in,out,nin,nout)
      real in(nin,nin),out(nout,nout)
c copy square 2d array 'in' to 2d array 'out'.
c if in is smaller than out, pad out with zeros; if it is larger truncate.
      do j=1,min(nin,nout)
         do i=1,min(nin,nout)
            out(i,j)=in(i,j)
         enddo
         do i=nin+1,nout
            out(i,j)=0
         enddo
      enddo
      do j=nin+1,nout
         do i=1,nout
            out(i,j)=0
         enddo
      enddo
      return
      end


      subroutine shfilt(sh,filt,mft)
c compute 2D Fourier filter that corresponds to convolving with shapelet sh.
      parameter(twopi=2*3.1415926535)
      dimension xtab(mm+1),ytab(mm+1),sh(msh)
      complex filt(mft/2,mft),sum

      beta=sh(1)
      n=sh(2)
c construct shapelet function at frequency fxnorm, fynorm
c Hn(x/beta)exp(-0.5*x**2/beta**2) --> i^n Hn(f*beta) exp(-0.5*f**2*beta**2)
      fmax2=(sqrt(n*1.)+10)**2
      do iy=1,mft
         fynorm=(mod(iy-1+mft/2,mft)-mft/2)*twopi*beta/mft
         call htab(fynorm,n,ytab)
         do ix=1,mft/2
            fxnorm=(mod(ix-1+mft/2,mft)-mft/2)*twopi*beta/mft
            call htab(fxnorm,n,xtab)
            fnorm2=fxnorm*fxnorm+fynorm*fynorm
            sum=0
            if (fnorm2.gt.fmax2) goto 7 
            fac=exp(-0.5*fnorm2)/sqrt(3.1415926535)*beta*twopi
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
                  if (sh(k+2).ne.0.) then
                     sum=sum+(0,1)**nn*sh(k+2)*fac*xtab(i+1)*ytab(j+1)
                  endif
               enddo
               fac=fac/sqrt(2.)
            enddo
c     put result in filt array. Later when arrays are FT'd back, result of
c     inverse FT is multiplied by mft*mft/2 - catch this factor here
 7          filt(ix,iy)=sum/(mft*mft/2)
         enddo
      enddo
      return
      end


      subroutine convolve(ftin,mft,filt,out,m)
c convolve slice iin of 3D array in with shapelet sh, result in array
c out. Input the FT of in, after zero-padding to dimension mft*mft,
c computed with rlft3.  output has dimension m*m.
      
      parameter(twopi=2*3.1415926535)
      complex ftin(mft/2,mft),speq(mft),filt(mft/2,mft)
      dimension out(m,m),outpadded(128*128)

c transfer FT to working array, treated as 2D mft*mft array.
      call copyslice(ftin,outpadded,mft)

c     Filter this FT with the filt array. 

      call multc(ftin,filt,outpadded,mft)

c set Nyquist freq. to zero and inverse transform.
      do j=1,mft
         speq(j)=0
      enddo
      call rlft3(outpadded,speq,mft,mft,1,-1)
c copy result into the output array, clipping padded borders.
      call copypadsq(outpadded,out,mft,m)

      return
      end


      subroutine multc(a,b,out,m)
c complex multiply arrays a*b=out
      complex a(m/2,m),b(m/2,m),out(m/2,m)
      do j=1,m
         do i=1,m/2
            out(i,j)=a(i,j)*b(i,j)
         enddo
      enddo
      return
      end
