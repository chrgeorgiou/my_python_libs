
      parameter (mstar=15,mstar2=2*mstar+1,maxbeta=25,maxround=3)
      parameter (maxmod=mmpoly2*(3+maxbeta*maxround*3*3))
      real psfcpts(mstar2,mstar2,maxmod),patch(mstar2,mstar2)
      real chi2(maxbeta),dx(2,maxbeta),cn(4,maxbeta),e(2,maxbeta),esig(3,maxbeta)


      common /mrqfunc/ cpts(mstar2*mstar2,30)

c read source list for image inimage.fits, and pixellated PSF model
c convolved with shapelets from pixpsfshearcpts.fits. 
c
c Fit the sources to sheared, circular shapelet models at different sizes.


c read in the PSF-convolved shapelets
      call get3('pixpsfshearcpts.fits',nx,ny,ncpt,psfcpts,mstar2,mstar2,maxmod)
c reverse engineer number of cpts. Should really put this in header information
      mmod=ncpt/(3+maxbeta*maxround*9)
      nfit=0
      do while ((nfit+1)*(nfit+2)/2.lt.mmod)
         nfit=nfit+1
      enddo
      write(0,*) 'Polynomial order for PSF variation, no.cpts:',nfit,mmod
c open the image containing the sources
      call openfits('inimage.fits',51,nx,ny)
      bgimage=get_bg_sec(51,nx*ny)
      write(0,*) 'BG=',bgimage
      pixerr=guess_errors_sec(51,nx*ny)
      xm=nx
      ym=ny


c now go through the source list.
      do isource=1,1000000
 1       read(*,*,err=1,end=2) x,y,f,fe,fw,a,b,pa,id,iflg
c get image patch, subtract image background
         isok=igetsec(51,nx,ny,nint(x)-mstar,nint(y)-mstar,
     :        mstar2,mstar2,patch,mstar2*mstar2)
         if (isok.ne.0) goto 1
         do iy=1,mstar2
            do ix=1,mstar2
               patch(ix,iy)=patch(ix,iy)-bgimage
            enddo
         enddo
c assign fitting radius based on catalogue
         fitrad=min(2*fw,mstar-1.)
c align the catalogued position with the center of the grid
         xshift=x-nint(x)
         yshift=y-nint(y)
         call shiftpatch(patch,mstar2,mstar2,xshift,yshift)
c put the data to be fitted, within radius fitrad of source, into a 1d list
         call loaddata(patch,fitrad,pixerr)
c go through the list of beta's and find the best-fit model for
c each. catalogue chi^2 values at peak and solution.
         ibestbeta=0
         bestch2=1e31
         do ibeta=1,maxbeta
c tabulate PSF-convolved components using PSF at source position.
            call loadcpts(psfcpts,x,y,xm,ym,nfit,ibeta-1)
c solve for best-fit coefficients.
            call bestfit(cn(1,ibeta),dx(1,ibeta),e(1,ibeta),chi2(ibeta),esig(1,ibeta))
            if (chi2(ibeta).lt.bestch2) then
               ibestbeta=ibeta
               bestch2=chi2(ibeta)
            endif
c            write(0,'(i4,g11.4,5f8.4,x,2f7.3)') ibeta,chi2(ibeta),
c     :           (e(k,ibeta),k=1,2),(esig(k,ibeta),k=1,3),(dx(k,ibeta),k=1,2)
         enddo
c apply shift of best fit so far, and repeat with that beta. Iterate 3x
         ibeta=ibestbeta
         beta=2.**((ibeta-1)*0.125)
         do itercenter=1,3
            call shiftpatch(patch,mstar2,mstar2,dx(1,ibeta),dx(2,ibeta))
            call loaddata(patch,fitrad,pixerr)
            call loadcpts(psfcpts,x,y,xm,ym,nfit,ibeta-1)
            call bestfit(cn(1,ibeta),dx(1,ibeta),e(1,ibeta),chi2(ibeta),esig(1,ibeta))
c            write(0,'(i4,g11.4,5f8.4,x,2f7.3)') ibeta,chi2(ibeta),
c     :           (e(k,ibeta),k=1,2),(esig(k,ibeta),k=1,3),(dx(k,ibeta),k=1,2)
         enddo
c         write(0,*)
         write(*,1002) x,y,f,fe,fw,a,b,pa,(e(k,ibeta),k=1,2),
     :        (esig(k,ibeta),k=1,3),id,iflg,beta,(cn(k,ibeta),k=1,4),
     :        (dx(k,ibeta),k=1,2)
 1002    format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,4g12.4,f8.4,i10,i5,
     :        f7.2,6g12.4)
         if (mod(isource,100).eq.0) write(0,*) 'Done source ',isource
      enddo
 2    continue
      nsource=isource-1

      end

      subroutine bestfit(cn,dx,e,chi2,esig)
      parameter (mstar=15,mstar2=2*mstar+1)
      real cn(4),dx(2),e(2),esig(3)
      real coefs(8),cov(8,8),alpha(8,8)
      external psfxsh

      common /pixtofit/ source(mstar2*mstar2),x(mstar2*mstar2),
     :     sig(mstar2*mstar2),ntofit

c assume the common block has been loaded with the right coefs, do a chi^2 fit 

      coefs(1)=1
      coefs(2)=1
      coefs(3)=1
      coefs(4)=1
      coefs(5)=0
      coefs(6)=0
      coefs(7)=0
      coefs(8)=0
      alam=-1
      iter=0
      do while (alam.lt.10000000. .and. iter.lt.50)
         call mrqmin(x,source,sig,ntofit,coefs,8,cov,alpha,8,chisq,
     :        psfxsh,alam)
c         write(0,'(10g11.3)') chisq,alam,(coefs(i),i=1,8)
         iter=iter+1
      enddo
      alam=0
      call mrqmin(x,source,sig,ntofit,coefs,8,cov,alpha,8,chisq,
     :     psfxsh,alam)
c      write(0,'(5g12.4,4f8.4)') chisq,(coefs(i),i=1,8)
c      write(0,*) 'Correlation matrix:'
c      do i=1,8
c         write(0,'(12f7.3)') (cov(i,j)/sqrt(cov(i,i)*cov(j,j)),j=1,8)
c      enddo
      do i=1,4
         cn(i)=coefs(i)
      enddo
      e(1)=coefs(5)
      e(2)=coefs(6)
      dx(1)=coefs(7)
      dx(2)=coefs(8)
      esig(1)=sqrt(cov(5,5))
      esig(2)=sqrt(cov(6,6))
      esig(3)=cov(5,6)/(esig(1)*esig(2))
      chi2=chisq/ntofit

c .........
      return
      end


      subroutine loaddata(pix,r,zerr)
c make list of pixel index and values to fit, within radius r of center of patch
      parameter (mstar=15,mstar2=2*mstar+1)
      real pix(-mstar:mstar,-mstar:mstar)
      common /pixtofit/ source(mstar2*mstar2),x(mstar2*mstar2),
     :     sig(mstar2*mstar2),ntofit

      npix=0
      do iy=-mstar,mstar
         do ix=-mstar,mstar
            if (ix*ix+iy*iy*1..le.r*r) then
               npix=npix+1
               x(npix)=(iy+mstar)*mstar2+ix+mstar+1
               source(npix)=pix(ix,iy)
               sig(npix)=zerr
            endif
         enddo
      enddo
      ntofit=npix
c      write(0,*) 'Loaded data for ',ntofit,' pixels.'
c      write(0,*) (x(i),source(i),sig(i),i=1,20)
      return
      end


      subroutine loadcpts(psfcpts,x,y,xm,ym,nfit,ibeta)
      parameter (mstar=15,mstar2=2*mstar+1)
      real psfcpts(mstar2*mstar2,*),wts((nfit+1)*(nfit+2)/2)
      real legendremodel

      common /mrqfunc/ cpts(mstar2*mstar2,30)

      npoly=(nfit+1)*(nfit+2)/2
      ipoly=0
      do k=0,nfit
         do m=0,k
            l=k-m
            ipoly=ipoly+1
            wts(ipoly)=legendremodel(x,y,xm,ym,l,m)
         enddo
      enddo
c      write(0,*) x,y,xm,ym
c      write(0,*) (wts(ipoly),ipoly=1,npoly)
c      write(0,*)
      call zeroarray(cpts,mstar2*mstar2*30)
c set the PSF, and its x and y translations in first three slices
      do i=1,npoly
         do iterm=0,2
            do ipix=1,mstar2*mstar2
               cpts(ipix,iterm+1)=cpts(ipix,iterm+1) + wts(i)*
     :              psfcpts(ipix,i+iterm*npoly)
            enddo
         enddo
      enddo
c now load the shapelet-convolved PSFs for the beta bin ibeta. First
c beta bin has ibeta=0
      do i=1,npoly
         do iterm=0,26
            do ipix=1,mstar2*mstar2
               cpts(ipix,iterm+4)=cpts(ipix,iterm+4) + wts(i)*
     :              psfcpts(ipix,i+3*npoly+iterm*npoly+ibeta*27*npoly)
            enddo
         enddo
      enddo
      return
      end


      subroutine zeroarray(a,m)
      real a(m)
      do i=1,m
         a(i)=0
      enddo
      return
      end


      subroutine psfxsh(pix,coefs,ymod,ygrad,ma)
      dimension coefs(ma),ygrad(ma)

      parameter (mstar=15,mstar2=2*mstar+1)
      common /mrqfunc/ cpt(mstar2*mstar2,30)

c routine that evaluates model (shapelet x PSF) at pixel 'pix'
c
c common block must be preloaded with the PSF-convolved components, as follows:
c
c PSF, T1.PSF, T2.PSF
c PSF-convolved C0, S1.C0, S2.C0, T1.C0, T1.S1.C0, T1.S2.C0, T2.C0, T2.S1.C0, T2.S2.C0
c ditto for C2
c ditto for C4
c
c the order of the coefficients is cp, c0, c2, c4, e1, e2, dx, dy
c
c and the model is 
c
c cp [PSF+dx T1.PSF+dy T2.PSF] + PSF * [1+dx T1+dy T2].[1+e1 S1+e2 S2].[c0C0+c2C2+c4C4]

      i=pix

      cp=coefs(1)
      c0=coefs(2)
      c2=coefs(3)
      c4=coefs(4)
      e1=coefs(5)
      e2=coefs(6)
      d1=coefs(7)
      d2=coefs(8)
      
      psf=      cp*cpt(i,1)
      t1psf=    cp*cpt(i,2)
      t2psf=    cp*cpt(i,3)

      circ=     c0*cpt(i,4)  + c2*cpt(i,13) + c4*cpt(i,22)
      s1circ=   c0*cpt(i,5)  + c2*cpt(i,14) + c4*cpt(i,23)
      s2circ=   c0*cpt(i,6)  + c2*cpt(i,15) + c4*cpt(i,24)
      t1circ=   c0*cpt(i,7)  + c2*cpt(i,16) + c4*cpt(i,25)
      t2circ=   c0*cpt(i,8)  + c2*cpt(i,17) + c4*cpt(i,26)
      t1s1circ= c0*cpt(i,9)  + c2*cpt(i,18) + c4*cpt(i,27)
      t1s2circ= c0*cpt(i,10) + c2*cpt(i,19) + c4*cpt(i,28)
      t2s1circ= c0*cpt(i,11) + c2*cpt(i,20) + c4*cpt(i,29)
      t2s2circ= c0*cpt(i,12) + c2*cpt(i,21) + c4*cpt(i,30)

      ymod= psf + d1*t1psf + d2*t2psf
     :     + circ + e1*s1circ + e2*s2circ
     :     + d1*t1circ + d2*t2circ + d1*e1*t1s1circ 
     :     + d1*e2*t1s2circ  + d2*e1*t2s1circ  + d2*e2*t2s2circ 
c dydcp
      ygrad(1)= cpt(i,1) + d1*cpt(i,2) + d2*cpt(i,3)
c dydc0
      ygrad(2)= cpt(i,4) + e1*cpt(i,5) + e2*cpt(i,6)
     :     + d1*cpt(i,7) + d2*cpt(i,8) + d1*e1*cpt(i,9)
     :     + d1*e2*cpt(i,10) + d2*e1*cpt(i,11) + d2*e2*cpt(i,12)
c dydc2
      ygrad(3)= cpt(i,13) + e1*cpt(i,14) + e2*cpt(i,15)
     :     + d1*cpt(i,16) + d2*cpt(i,17) + d1*e1*cpt(i,18)
     :     + d1*e2*cpt(i,19) + d2*e1*cpt(i,20) + d2*e2*cpt(i,21)
c dydc4
      ygrad(4)= cpt(i,22) + e1*cpt(i,23) + e2*cpt(i,24)
     :     + d1*cpt(i,25) + d2*cpt(i,26) + d1*e1*cpt(i,27)
     :     + d1*e2*cpt(i,28) + d2*e1*cpt(i,29) + d2*e2*cpt(i,30)
c dyde1
      ygrad(5)= s1circ + d1*t1s1circ + d2*t2s1circ
c dyde2
      ygrad(6)= s2circ + d1*t1s2circ + d2*t2s2circ
c dydd1
      ygrad(7)= t1psf + t1circ + e1*t1s1circ + e2*t1s2circ
c dydd2
      ygrad(8)= t2psf + t2circ + e1*t2s1circ + e2*t2s2circ

      return
      end



c========================================================
c DEBUGGING - NO LONGER NEEDED

c            call testfunc(481.,0.2,0.3,0.4,0.5,0.1,0.05,0.4,0.8)


      subroutine testfunc(pix,cp,c0,c2,c4,e1,e2,d1,d2)
      real coef(8)
      coef(1)=cp
      coef(2)=c0
      coef(3)=c2
      coef(4)=c4
      coef(5)=e1
      coef(6)=e2
      coef(7)=d1
      coef(8)=d2

      do i=1,8
         call compgrad(pix,coef,i,0.001)
      enddo
      return
      end


      subroutine compgrad(x,c,i,del)
      real dy(8),dy0(8),c(8)

      call psfxsh(x,c,y0,dy,8)
      c(i)=c(i)+del
      call psfxsh(x,c,y1,dy0,8)
      c(i)=c(i)-del
      write(0,*) x,i,y0,dy(i),(y1-y0)/del
      return
      end




      function ch2one(zz,fitrad,zmod,ix,iy,m)
c compute chi^2 for data zz, using zmod shifted by ix,iy, within radius fitrad

      real zz(-m:m,-m:m),zmod(-m:m,-m:m)
      double precision summm,summd,ch2

      n=min(m-1,int(fitrad+1))
      summm=0
      summd=0
      do j=-n,n
         do i=-n,n
            if (real(i*i+j*j).le.fitrad*fitrad) then
               summm=summm+zmod(i-ix,j-iy)**2
               summd=summd+zmod(i-ix,j-iy)*zz(i,j)
            endif
         enddo
      enddo
      amp=summd/summm
      ch2=0
      do j=-n,n
         do i=-n,n
            if (i*i+j*j.le.fitrad*fitrad) then
               ch2=ch2+(amp*zmod(i-ix,j-iy)-zz(i,j))**2
            endif
         enddo
      enddo
      ch2one=ch2
      return
      end



