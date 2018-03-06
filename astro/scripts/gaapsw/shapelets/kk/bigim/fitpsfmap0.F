       real sh(msh),sharr(mmcat,msh-2),x(mmcat),y(mmcat),
     :    dx(mmcat),dy(mmcat),coefs(mmpoly2)
      real w(mmcat),w0(mmcat),bg(mmcat),totsum(mmcat)
      real allcoef(mmpoly2,msh-2)
      real bestfit(mmcat,msh-2),resid(mmcat)
      parameter (pi=3.1415926535)
      integer id(mmcat)
      character*11 toplbl

      do i=1,mmcat
    1     read(*,*,err=1,end=2) x(i),y(i),f,fe,fw,a,b,th,id(i),iflg
          read(*,*) totsum(i),dx(i),dy(i),bg(i),tt,tt
          read(*,*) beta,n
          read(*,*) (sharr(i,k),k=1,(n+1)*(n+2)/2)
          read(*,*)
      enddo
    2 continue
      npsf=i-1

      call readorders(n2,nfit)
      n=min(n,n2)

c estimate the extent of the images from the maximum x and y values found. Details not important
      xmax=-1e30
      ymax=-1e30
      do i=1,npsf
         xmax=max(xmax,x(i))
         ymax=max(ymax,y(i))
      enddo
      xmax=xmax/(1-0.5/npsf)
      ymax=ymax/(1-0.5/npsf)
         
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
     :        coefs,bestfit(1,k))
c preserve the weights (indicating what has been clipped) for 00 cpt fit
          if (k.eq.1) then
              do i=1,npsf
                  w0(i)=w(i)
              enddo
          endif
          write(*,1002) (coefs(kfit),kfit=1,(nfit+1)*(nfit+2)/2)
          write(*,*)
 1002     format (mmpoly2(g15.6))
          do kfit=1,(nfit+1)*(nfit+2)/2
              allcoef(kfit,k)=coefs(kfit)
          enddo
          do i=1,npsf
              resid(i)=sharr(i,k)-bestfit(i,k)
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

c write out a file listing the stars that were used in the fit.
      open(17,file='starsused.txt',status='unknown')
      do i=1,npsf
          if (w0(i).gt.0.) then
              write(17,*) id(i)
          endif
      enddo
      close(17)

      call pgbeg(0,'psfsticks.ps/ps',2,1)
      call pgsch(1.5)
      call pgslw(2)
      call pgenv(-0.05*xmax,xmax*1.05,-0.05*ymax,ymax*1.05,1,0)
      call pglab('X (pix)','Y (pix)','PSF ellipticities--raw')
      scale=0.3*xmax
      call pgsci(3)
      call pgmove(0.12*xmax,1.02*ymax)
      call pgdraw(0.12*xmax+0.05*scale,1.02*ymax)
      call pgtext(0.13*xmax+0.05*scale,1.02*ymax,'5%')
      do i=1,npsf
          call pgsci(nint(1+w0(i)))
          a00=sharr(i,1)
          a10=sharr(i,2)
          a01=sharr(i,3)
          a20=sharr(i,4)
          a11=sharr(i,5)
          a02=sharr(i,6)
c old values - unweighted 2nd moments?
          txx=beta**3*sqrt(pi)*(2*a00+sqrt(2.)*(5*a20+a02))
          txy=beta**3*sqrt(pi)*4*a11
          tyy=beta**3*sqrt(pi)*(2*a00+sqrt(2.)*(a20+5*a02))
c these are the second moments, weighted by gaussian wt fn of width beta
c     first three not needed here, but they are wt'ed 0th and 1st moments
c          t=beta*sqrt(pi)*a00
c          tx=beta**2*sqrt(pi)*sqrt(0.5)*a10
c          ty=beta**2*sqrt(pi)*sqrt(0.5)*a01
          txx=beta**3*sqrt(pi)*(0.5*a00+sqrt(0.5)*a20)
          txy=beta**3*sqrt(pi)*0.5*a11
          tyy=beta**3*sqrt(pi)*(0.5*a00+sqrt(0.5)*a02)
          e1=(txx-tyy)/(txx+tyy)
          e2=2*txy/(txx+tyy)
          ell=sqrt(e1*e1+e2*e2)
          ang=atan2(e2,e1)/2
          dxe=ell*cos(ang)*scale/2
          dye=ell*sin(ang)*scale/2
c          write(0,*) txx,txy,tyy,e1,e2,ell,ang
          call pgmove(x(i)-dxe,y(i)-dye)
          call pgdraw(x(i)+dxe,y(i)+dye)
      enddo
c now plot the residuals from the fit as sticks as well.
      call pgenv(-0.05*xmax,xmax*1.05,-0.05*ymax,ymax*1.05,1,0)
      call pglab('X (pix)','Y (pix)','PSF ellipticities--resid')
      scale=xmax
      call pgsci(3)
      call pgmove(0.12*xmax,1.02*ymax)
      call pgdraw(0.12*xmax+0.01*scale,1.02*ymax)
      call pgtext(0.13*xmax+0.01*scale,1.02*ymax,'1%')
      do i=1,npsf
          call pgsci(nint(1+w0(i)))
          a00=sharr(i,1)
          a10=sharr(i,2)
          a01=sharr(i,3)
          a20=sharr(i,4)
          a11=sharr(i,5)
          a02=sharr(i,6)
          call interpolpsf(x(i),y(i),xmax,ymax,allcoef,n,beta,nfit,sh)
          f00=sh(3)
          f10=sh(4)
          f01=sh(5)
          f20=sh(6)
          f11=sh(7)
          f02=sh(8)
c old values - unweighted 2nd moments?
c          txx=beta**3*sqrt(pi)*(2*a00+sqrt(2.)*(5*a20+a02))
c          txy=beta**3*sqrt(pi)*4*a11
c          tyy=beta**3*sqrt(pi)*(2*a00+sqrt(2.)*(a20+5*a02))
c these are the second moments, weighted by gaussian wt fn of width beta
c     first three not needed here, but they are wt'ed 0th and 1st moments
c          t=beta*sqrt(pi)*a00
c          tx=beta**2*sqrt(pi)*sqrt(0.5)*a10
c          ty=beta**2*sqrt(pi)*sqrt(0.5)*a01
c     2nd moments from fit
          txx=beta**3*sqrt(pi)*(0.5*f00+sqrt(0.5)*f20)
          txy=beta**3*sqrt(pi)*0.5*f11
          tyy=beta**3*sqrt(pi)*(0.5*f00+sqrt(0.5)*f02)
c     original 2nd moments
          pxx=beta**3*sqrt(pi)*(0.5*a00+sqrt(0.5)*a20)
          pxy=beta**3*sqrt(pi)*0.5*a11
          pyy=beta**3*sqrt(pi)*(0.5*a00+sqrt(0.5)*a02)
          e1=(pxx-txx-pyy+tyy)/(txx+tyy)
          e2=2*(pxy-txy)/(txx+tyy)
          ell=sqrt(e1*e1+e2*e2)
          ang=atan2(e2,e1)/2
          dxe=ell*cos(ang)*scale/2
          dye=ell*sin(ang)*scale/2
c          write(0,*) txx,txy,tyy,e1,e2,ell,ang
          call pgmove(x(i)-dxe,y(i)-dye)
          call pgdraw(x(i)+dxe,y(i)+dye)
      enddo
      call pgend
      end


