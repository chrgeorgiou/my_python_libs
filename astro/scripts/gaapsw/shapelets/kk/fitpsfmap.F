      real sh(msh),sharr(mmcat,msh-2),x(mmcat),y(mmcat),
     :    dx(mmcat),dy(mmcat),coefs(mmpoly2),shm(msh)
      real w(mmcat),w0(mmcat),bg(mmcat),totsum(mmcat)
      real allcoef(mmpoly2,msh-2)
      real bestfit(mmcat),resid(mmcat)
      parameter (maxstamp=800,maxstamp2=maxstamp*maxstamp)
      real stamp(maxstamp2)
      integer id(mmcat)
      character*11 toplbl

      call openfits('inimage.fits',51,nx,ny)

      xmax=nx
      ymax=ny
      do i=1,mmcat
    1     read(*,*,err=1,end=2) x(i),y(i),f,fe,fw,a,b,th,id(i),iflg
          if (x(i).gt.xmax .or. y(i).gt.ymax) then
              write(0,*) 'PSF STAR OUTSIDE IMAGE!!!  *** STOP ***'
              stop
          endif
          read(*,*) totsum(i),dx(i),dy(i),bg(i),tt,tt
          read(*,*) beta,n
          read(*,*) (sharr(i,k),k=1,(n+1)*(n+2)/2)
          read(*,*)
      enddo
    2 continue
      npsf=i-1

      call readorders(n2,nfit)
      n=min(n,n2)

      write(*,*) n,beta
      write(*,*) nfit,xmax,ymax
      write(*,*)
      call pgbeg(0,'fitpsfmap.ps/cps',4,3)
      call pgsch(1.5)
      k=0
c     set initial weights to 1. Later will only use those stars that
c     were used for the 00 fit.
      do i=1,npsf
          w(i)=1
      enddo
      do kk=0,n
        do m=0,kk
          l=kk-m
          k=k+1
          nrejiter=3
          clip=4.
c if this is not the 00 component, set the weights to those of the 00 cpt 
          if (k.ne.1) then 
              do i=1,npsf
                  w(i)=w0(i)
              enddo
          endif
          call fitpoly(npsf,xmax,ymax,x,y,sharr(1,k),w,nfit,nrejiter,clip,
     :        coefs,bestfit)
c preserve the weights (indicating what has been clipped) for 00 cpt fit
          if (k.eq.1) then
              do i=1,npsf
                  w0(i)=w(i)
              enddo
          endif
          write(*,*) (coefs(kfit),kfit=1,(nfit+1)*(nfit+2)/2)
          write(*,*)
          do kfit=1,(nfit+1)*(nfit+2)/2
              allcoef(kfit,k)=coefs(kfit)
          enddo
          do i=1,npsf
              resid(i)=sharr(i,k)-bestfit(i)
          enddo
          tmin=-0.0002
          tmax= 0.0002
          do i=1,npsf
              tmin=min(tmin,sharr(i,k))
              tmax=max(tmax,sharr(i,k))
              tmin=min(tmin,resid(i))
              tmax=max(tmax,resid(i))
          enddo
          dt=0.1*(tmax-tmin)
          tmin=tmin-dt
          tmax=tmax+dt
c bound the plots to something sensible
          if (kk.eq.0) then
              tmax=min(1.,tmax)
              tmin=max(-0.5,tmin)
          else
              tmax=min(tmax,0.5)
              tmin=max(tmin,-0.5)
          endif
          
          call pgsch(1.5)
          call pgenv(0.,xmax,tmin,tmax,0,0)
          call pgmove(0.,0.)
          call pgdraw(xmax,0.)
          call pgsch(2.)
          write(toplbl,'(''l='',i2,''   m='',i2)') l,m
          call pglab('X',' ',toplbl)
          do i=1,npsf
              call pgsch(1.+2*w(i))
              call pgsci(2)
              call pgpt1(x(i),sharr(i,k),5)
              call pgsci(1)
              call pgpt1(x(i),resid(i),21)
          enddo

          call pgsch(1.5)
          call pgenv(0.,ymax,tmin,tmax,0,0)
          call pgmove(0.,0.)
          call pgdraw(ymax,0.)
          call pgsch(2.)
          call pglab('Y',' ',toplbl)
          do i=1,npsf
              call pgsch(1.+2*w(i))
              call pgsci(2)
              call pgpt1(y(i),sharr(i,k),5)
              call pgsci(1)
              call pgpt1(y(i),resid(i),21)
          enddo
          call pgsch(1.5)
        enddo
      enddo
      call pgend

c     now make plots for all stars of the shapelet expansion, and the
c     interpolated model.

      nplot=max(10,nint(6*beta))
      call pgbeg(0,'psfstars.ps/cps',11,8)
      call pgsch(1.5)
      call pgscir(16,49)
      fitradius=max(10.,min(40.,8*beta))
      do i=1,npsf
          sh(1)=beta
          sh(2)=n
          do k=1,(n+1)*(n+2)/2
              sh(k+2)=sharr(i,k)
          enddo
          if (mod(i,2).eq.0) call pgadvance
c display the patch around the star as a greyscale
          pk=totsum(i)*sh(3)/beta/sqrt(3.1416)
          npix=3*nplot+1
          ix1=nint(x(i))-(npix-1)/2
          ix1=max(1,min(ix1,nx-npix+1))
          iy1=nint(y(i))-(npix-1)/2
          iy1=max(1,min(iy1,ny-npix+1))
          xoff=x(i)-ix1+1
          yoff=y(i)-iy1+1
          isok=igetsec(51,nx,ny,ix1,iy1,npix,npix,stamp,maxstamp2)
          if (isok.eq.0) then
             call showpatch(stamp,npix,npix,npix,npix,xoff,yoff,real(3*nplot),bg(i),bg(i)+pk)
c     subtract the shapelet model from image, display as greyscale
             shm(1)=sh(1)
             shm(2)=sh(2)
             do k=1,(n+1)*(n+2)/2
                shm(k+2)=-totsum(i)*sh(k+2)
             enddo
             call addtoimage(shm,xoff-dx(i),yoff-dy(i),fitradius,stamp,npix,npix)
             call showpatch(stamp,npix,npix,npix,npix,xoff,yoff,real(3*nplot),bg(i)-0.1*pk,bg(i)+0.1*pk)
          else
             call pgadvance
             call pgadvance
          endif
c show contour plot of the star's shapelet model
          call plotshape(sh,nplot,'REAL')
c show contour plot of the interpolated PSF model
          call interpolpsf(x(i),y(i),xmax,ymax,allcoef,n,beta,nfit,sh)
          call plotshape(sh,nplot,'MODEL')

          call pgsch(0.5)
          call pgenv(0.,1.,0.,1.,1,-2)
          call pgsch(3.)
          if (w0(i).eq.0.) then
              call pgsci(2)
              call pgptxt(0.4,0.5,45.,0.5,'ZAPPED!')
              call pgsci(1)
          endif
          write(toplbl,'(''ID='',i8)') id(i)
          call pgptxt(0.4,0.7,0.,0.4,toplbl)
          write(toplbl,'(''X ='',f8.1)') x(i)
          call pgptxt(0.4,0.4,0.,0.4,toplbl)
          write(toplbl,'(''Y ='',f8.1)') y(i)
          call pgptxt(0.4,0.2,0.,0.4,toplbl)
          call pgarro(0.4,0.6,0.,0.6)
          call pgsch(1.5)
      enddo
      call pgend

      call pgbeg(0,'psfstarpos.ps/ps',1,1)
      call pgsch(1.5)
      call pgslw(2)
      call pgenv(0.,xmax*1.05,0.,ymax*1.05,1,0)
      call pglab('X (pix)','Y (pix)','Location of PSF stars')
      call pgsch(3.)
      do i=1,npsf
          isymbol=nint(5+16*w0(i))
          call pgpt1(x(i),y(i),isymbol)
      enddo
      call pgend

      write(0,*) 'Fitpsfmap found ',npsf,' stars, fitted to x+y order ',nfit


      end


