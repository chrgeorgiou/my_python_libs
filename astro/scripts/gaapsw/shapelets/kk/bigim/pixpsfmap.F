      parameter (mstar=15,mstar2=2*mstar+1,m2=32,nstarmax=10000)
      dimension starpatch(-mstar:mstar,-mstar:mstar+2*mstar2,nstarmax)
c      dimension zz(m2,m2)
c      dimension starmod(-mstar:mstar,-mstar:mstar,nstarmax)
      dimension starcoefs(-mstar:mstar,-mstar:mstar,mmpoly2)
      dimension x(nstarmax),y(nstarmax),beta(nstarmax),flx(nstarmax),
     :     dx(nstarmax),dy(nstarmax),ch2(nstarmax),id(nstarmax)
      dimension coefs(mmpoly2),tofit(nstarmax),bestfit(nstarmax),
     :     w(nstarmax),nlost(nstarmax)
      dimension shape(msh)

      nim=51
      call openfits('inimage.fits',nim,nx,ny)
      bgimage=get_bg_sec(nim,nx*ny)
      write(0,*) 'BG=',bgimage
c      err=guess_errors_sec(51,nx*ny)
      call pgbeg(0,'pixpsfstars.ps/ps',8,6)
      call pgask(.false.)

      do istar=1,nstarmax
    1     read(*,*,err=1,end=2) x(istar),y(istar),(ttt,i=3,8),id(istar)
          ix=nint(x(istar)-mstar)
          iy=nint(y(istar)-mstar)
          xc=x(istar)-ix+1
          yc=y(istar)-iy+1
          isok=igetsec(nim,nx,ny,ix,iy,mstar2,mstar2,
     :        starpatch(-mstar,-mstar,istar),mstar2*mstar2)
          if (isok.ne.0) then
              write(0,*) isok, x(istar),y(istar)
              call ftrprt('stderr',isok)
              goto 1
          endif
          do iy=-mstar,mstar
             do ix=-mstar,mstar
                starpatch(ix,iy,istar)=starpatch(ix,iy,istar)-bgimage
             enddo
          enddo
          bg=1.e-31
c          call fitbg(starpatch(-mstar,-mstar,istar),
c     :        mstar2,mstar2,mstar2,xc,yc,mstar*1.,bg)
          call fitbeta(starpatch(-mstar,-mstar,istar),
     :        mstar2,mstar2,mstar2,xc,yc,mstar*1.,bg,beta(istar))
      enddo
    2 nstar=istar-1
      write(0,*) 'Read ',nstar,' stars.'

      betapsf=fmedian(nstar,beta)*1.3
      write(0,*) 'Adopted beta ',betapsf
      nsh=8
      do istar=1,nstar
c assume background has been subtracted properly.
c otherwise fit it separately
c         call fitbg(starpatch(-mstar,-mstar,istar),
c     :        mstar2,mstar2,mstar2,mstar+1.,mstar+1.,mstar*1.,bg)
c iterate a low-order shapelets fit, to determine centroid shift to apply.
         delx0=x(istar)-nint(x(istar))
         dely0=y(istar)-nint(y(istar))
         delx=delx0
         dely=dely0
         do it=1,4
            xc=delx+mstar+1
            yc=dely+mstar+1
            call fitshapelets(starpatch(-mstar,-mstar,istar),
     :           mstar2,mstar2,mstar2,xc,yc,mstar*1.,bg,shape,2,betapsf)
            dxsh=betapsf*1.4142*(shape(4)/shape(3))
            dysh=betapsf*1.4142*(shape(5)/shape(3))
c            write(0,'(i7,6f12.4)') istar,x(istar),y(istar),delx,dely,dxsh,dysh
            delx=delx+dxsh
            dely=dely+dysh
         enddo
         dx(istar)=delx-delx0
         dy(istar)=dely-dely0

c shift star to position determined this way - enforces zero gauss-wt 1st moment
         call shiftpatch(starpatch(-mstar,-mstar,istar),mstar2,mstar2,delx,dely)
         xc=mstar+1
         yc=mstar+1
         
         call fitshapelets(starpatch(-mstar,-mstar,istar),
     :        mstar2,mstar2,mstar2,xc,yc,mstar*1.,bg,
     :        shape,nsh,betapsf)
         flx(istar)=seriessum(shape)
         do iy=-mstar,mstar
            do ix=-mstar,mstar
               starpatch(ix,iy,istar)=(starpatch(ix,iy,istar)-bg)/flx(istar)
            enddo
         enddo
         pk=starpatch(0,0,istar)
c display about 50 randomly selected stars on one page.
         if (mod(istar,(nstar/48+1)).eq.0) then
            call showpatch(starpatch(-mstar,-mstar,istar),
     :           mstar2,mstar2,mstar2,mstar2,xc,yc,14.,-0.1*pk,pk)
         endif
      enddo
      call pgend

c fit normalized stars as function of field position.

      norder=min(4,int(nstar/20))
      niter=3
      clip=4.
      open(12,file='pixmap.par',status='old',err=3)
      read(12,*,err=3) norder,niter,clip
      close(12)
 3    write(0,*) 'fit order, no. iters, sigclip: ',norder,niter,clip
      nmod=(norder+1)*(norder+2)/2
      do is=1,nstar
         nlost(is)=0
      enddo
      do iy=-mstar,mstar
         do ix=-mstar,mstar
            do is=1,nstar
               tofit(is)=starpatch(ix,iy,is)
               w(is)=sqrt(max(0.,flx(is)))
c               if (err.eq.0.) then
c                  w(is)=flx(is)
c               else
c                  w(is)=flx(is)**2/(err**2+max(0.,tofit(is))*flx(is))
c               endif
            enddo
c            write(0,*) 'Fitting pixel ',ix,iy
c            write(0,'(3(i8,4g10.3))') (is,x(is),y(is),tofit(is),w(is),is=1,nstar)
            call fitlegendre(nstar,nx*1.,ny*1.,x,y,tofit,w,norder,niter,clip,coefs,bestfit)
            do imod=1,nmod
               starcoefs(ix,iy,imod)=coefs(imod)
            enddo
c write out best fit to each star as middle panel of starpatch images
            do is=1,nstar
               starpatch(ix,iy+mstar2,is)=bestfit(is)
               if (w(is).eq.0.) nlost(is)=nlost(is)+1
            enddo
         enddo
         if (mod(iy,5).eq.0) write(0,*) 'Fitted row ',iy
      enddo

c show the final PSF components
      
      if (norder.eq.0) then
         nxpanel=1
         nypanel=1
      else
         nxpanel=norder+1
         nypanel=max(2,norder)
         if (norder.eq.2) nypanel=3
      endif
c      npanel=int(sqrt(nmod-0.001))+1
      call pgbeg(0,'pixpsfmodel.ps/ps',nxpanel,nypanel)
      call pgask(.false.)
      pk=starcoefs(0,0,1)
      zhi=pk
      zlo=-0.1*pk
      imod=0
      ix=0
      do n=0,norder
         if (ix+n+1.gt.nxpanel) then
            do k=1,nxpanel-ix
               call pgadvance
            enddo
            ix=0
         endif
         do n1=0,n
            imod=imod+1
            ix=ix+1
            call showpatch(starcoefs(-mstar,-mstar,imod),
     :           mstar2,mstar2,mstar2,mstar2,xc,yc,14.,zlo,zhi)
            zhi=0.1*pk
         enddo
         if (ix.lt.nxpanel) then
            call pgadvance
            ix=ix+1
         endif
         ix=mod(ix,nxpanel)
      enddo
      call pgend


c write out the PSF components as a 3D fits image.
      call put3('pixpsf.fits',mstar2,mstar2,nmod,
     :     starcoefs,mstar2,mstar2,nmod)

c make images of residuals between model and stars, as top panel of starpatch

      do istar=1,nstar
         pk=starpatch(0,0,istar)
         ch2tmp=0
         do iy=-mstar,mstar
            do ix=-mstar,mstar
               starpatch(ix,iy+2*mstar2,istar)=starpatch(ix,iy,istar)
     :              -starpatch(ix,iy+mstar2,istar)
               ch2tmp=ch2tmp+starpatch(ix,iy+2*mstar2,istar)**2
            enddo
         enddo
         ch2(istar)=sqrt(ch2tmp)*flx(istar)
      enddo

c write out PSF stars, models, and residuals, to a fits stack
      call put3('psfstars.fits',mstar2,3*mstar2,nstar,
     :     starpatch,mstar2,3*mstar2,nstar)

      write(*,*) '#',betapsf,norder,niter,clip
      do is=1,nstar
         write(*,1001) is,x(is),y(is),dx(is),dy(is),flx(is),starpatch(0,0,is),
     :        starpatch(0,2*mstar2,is),ch2(is),nlost(is),id(is)
      enddo
 1001 format(i7,2f11.3,2f8.4,g14.6,2f10.5,g14.5,i6,i10)
      end

