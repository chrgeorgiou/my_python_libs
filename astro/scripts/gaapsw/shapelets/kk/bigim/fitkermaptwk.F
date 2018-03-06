c     fit to double-shapelet kernel fits from psfcat2gauskerwithtweak
c     for each star there is an inner + outer shapelet fit
c     this program fits a polynomial to the variation of the coeffs and writes it out as a map
      
      parameter (mmpsf=20000)
      real sharr(mmpsf,2*msh-4),x(mmpsf),y(mmpsf),
     :    dx(mmpsf),dy(mmpsf),coefs(mmpoly2),sh(msh)
      real w(mmpsf),w0(mmpsf),bg(mmpsf),totsum(mmpsf)
      real allcoef(mmpoly2,2*msh-4)
      real bestfit(mmpsf,2*msh-4),resid(mmpsf)
      parameter (pi=3.1415926535)
      integer id(mmpsf)
      character*11 toplbl

      call openfits('inimage.fits',51,nx,ny)

      xmax=nx
      ymax=ny
      do i=1,mmpsf
    1     read(*,*,err=1,end=2) x(i),y(i),f,fe,fw,a,b,th,id(i),iflg
          if (x(i).gt.xmax .or. y(i).gt.ymax) then
              write(0,*) 'PSF STAR OUTSIDE IMAGE!!!  *** STOP ***'
              stop
          endif
          read(*,*) totsum(i),dx(i),dy(i),bg(i),tt,tt
          read(*,*) beta,n
          nin=(n+1)*(n+2)/2
          read(*,*) (sharr(i,k),k=1,nin)
          read(*,*) betaout,nout
          read(*,*) (sharr(i,k),k=1+nin,nin+(nout+1)*(nout+2)/2)
          read(*,*)
c     report shapelet sum for each star found
          sh(1)=beta
          sh(2)=n
          do k=1,nin
             sh(k+2)=sharr(i,k)
          enddo
          sum=seriessum(sh)
          sh(1)=betaout
          sh(2)=nout
          do k=1,(nout+1)*(nout+2)/2
             sh(k+2)=sharr(i,k+nin)
          enddo
          sum=sum+seriessum(sh)
          write(0,*) x(i),y(i),sum
      enddo
    2 continue
      npsf=i-1

      call readorders(n2,nfit)

      write(*,*) n,beta
      write(*,*) nfit,xmax,ymax
      write(*,*)
      call pgbeg(0,'fitkermaptwk.ps/cps',4,3)
      call pgsch(1.5)
      k=0
c     set initial weights to 1. Later will only use those stars that
c     were used for the 00 fit.
      do i=1,npsf
          w(i)=1
      enddo
      kk=0
      l=0
      m=0
      nrejiter=8
      clip=4.
      do k=1,nin+(nout+1)*(nout+2)/2 
c     if this is not the 00 component, set the weights to those of the 00 cpt 
         if (k.ne.1) then 
            do i=1,npsf
               w(i)=w0(i)
            enddo
         endif
         call fitpoly(npsf,xmax,ymax,x,y,sharr(1,k),w,nfit,nrejiter,clip,
     :        coefs,bestfit(1,k))
c     preserve the weights (indicating what has been clipped) for 00 cpt fit
         if (k.eq.1) then
            do i=1,npsf
               w0(i)=w(i)
            enddo
c no separate rejections once 00 component has been accepted. change?
            nrejiter=1
         endif
         write(*,1002) (coefs(kfit),kfit=1,(nfit+1)*(nfit+2)/2)
         write(*,*)
 1002    format (mmpoly2(g15.6))
         do kfit=1,(nfit+1)*(nfit+2)/2
            allcoef(kfit,k)=coefs(kfit)
         enddo
         do i=1,npsf
            resid(i)=sharr(i,k)-bestfit(i,k)
         enddo
         tmin=-0.0002
         tmax= 0.0002
         do i=1,npsf
            if (w(i).gt.0.) then
               tmin=min(tmin,sharr(i,k))
               tmax=max(tmax,sharr(i,k))
               tmin=min(tmin,resid(i))
               tmax=max(tmax,resid(i))
            endif
         enddo
         dt=0.1*(tmax-tmin)
         tmin=tmin-dt
         tmax=tmax+dt
c     bound the plots to something sensible
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
c     update counters for x and y orders. kk is sum of both. x-order counts down, y-order counts up           
         if (l.gt.0) then
            l=l-1
            m=m+1
         else
            kk=kk+1
            l=kk
            m=0
         endif
c     once hit nin models, reset to zero and start on the outer shapelet.
c     use the same nfit, xmax, ymax as for the inner shapelet (but write them out again - in principle could be different)         
         if (k.eq.nin) then
            kk=0
            l=0
            m=0
            write(*,*) '# OUTER SHAPELET COEFFS:'
            write(*,*) nout,betaout
            write(*,*) nfit,xmax,ymax
            write(*,*)
         endif
      enddo
      call pgend

      write(0,*) 'Fitpsfmap found ',npsf,' stars, fitted to x+y order ',nfit

      end


