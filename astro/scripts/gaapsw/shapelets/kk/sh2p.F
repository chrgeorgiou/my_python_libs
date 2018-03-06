c     from a shape catalogue, form the estimators that are like Pgamma
c     and e in KSB

c     they are:
c     P11 = S1 . Pinv . G | 20-02
c     P12 = S1 . Pinv . G | 11
c etc
c     and e1=Pinv . G | 20-02
c         e2=Pinv . G | 11
c where | indicates the coefficients to take.

c this corresponds to solving, for g,

c  (1 - g1 S1 - g2 S2 ) Pinv G = 0 in these coefs
c ie asking which shear needs to be applied to the deconvolved 
c galaxy Pinv G to make the 2nd-order cpts round.

c     errors on G can be propagated in solving this equation to give a
c     covariance for e1 and e2.

      real sh(msh),shpsf(msh),psfcoef(mmpoly2,msh-2),psfmat(msh-2,msh-2)
      real sh1(msh),sh2(msh)
      integer indx(msh-2)
c      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
c      real coefs(mm/2+2),covar(mm/2+2,mm/2+2)
      real p(0:mm),frac(2:6)
c      common /shearfns/ round,sheared1,sheared2,nshear

      open(10,file='psf.map',status='old')
      read(10,*) npsf,betapsf
      read(10,*) nfit,xmax,ymax
      read(10,*)
      do k=1,(npsf+1)*(npsf+2)/2
         read(10,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(10,*)
      enddo
      close(10)
      write(0,*) 'Read psf model from psf.map'
      write(0,*) 'Order for shape, spatial variation; beta:',npsf,nfit,betapsf

c set cuts on maximum fraction of power that can be found at order 2-6
      shiftmax=9
      do i=2,6
          frac(i)=1
      enddo
      open(10,file='shapecuts.dat',status='old',err=3)
      read(10,*) shiftmax,(frac(i),i=2,6)
      close(10)
    3 write(0,'('' Shift, shape (2..6) cuts: '',6f6.3)') 
     :    shiftmax,(frac(i),i=2,6)

      nread=0
      do i=1,mmcat
    1     read(*,*,err=1,end=99) x,y,f,fe,fw,a,b,th,id,iflg
          read(*,*) sum,dx,dy,bg,fitradius,sherr
          read(*,*) beta,n
          read(*,*) (sh(k),k=3,(n+1)*(n+2)/2+2)
          read(*,*)
          sh(1)=beta
          sh(2)=n
          nsh=(n+1)*(n+2)/2
c apply filters to the catalogue - set flags if shape looks dubious
c penalty flag is 100 for each shape constraint violated; 
c penalty is 1000 if fwhm is smaller than 1.1*PSF fwhm
          k=2
          sump=0
          do nn=0,n
              p(nn)=0
              do j=0,nn
                  k=k+1
                  p(nn)=p(nn)+sh(k)**2
              enddo
              sump=sump+p(nn)
          enddo
          if (dx*dx+dy*dy .gt. shiftmax*shiftmax) iflg=iflg+100
          do nn=2,min(n,6)
              if (p(nn).gt.frac(nn)*sump) iflg=iflg+100
          enddo
          if (fw.lt.1.1*betapsf*2.3) iflg=iflg+1000

          call interpolpsf(x,y,xmax,ymax,psfcoef,npsf,betapsf,nfit,shpsf)
c          call writeshapen('PSF',shpsf)
c          call writeshapen('Target',sh)

c deconvolve the galaxy shape - use same beta as for observed shape
          call setpsfmat(shpsf,beta,beta,psfmat)
          call ludcmp(psfmat,nsh,msh-2,indx,det)
          call lubksb(psfmat,nsh,msh-2,indx,sh(3))
          call shear1(sh,sh1)
          call shear2(sh,sh2)
          pg11=sh1(6)-sh1(8)
          pg12=sh2(6)-sh2(8)
          pg21=sh1(7)*sqrt(2.)
          pg22=sh2(7)*sqrt(2.)
          e1=sh(6)-sh(8)
          e2=sh(7)*sqrt(2.)
          do j=1,nsh+2
              sh(j)=0
              sh1(j)=0
              sh2(j)=0
          enddo
          sh(6)=1
          sh1(7)=1
          sh2(8)=1
          call lubksb(psfmat,nsh,msh-2,indx,sh(3))
          call lubksb(psfmat,nsh,msh-2,indx,sh1(3))
          call lubksb(psfmat,nsh,msh-2,indx,sh2(3))
          var66=0
          var77=0
          var88=0
          cov67=0
          cov68=0
          cov78=0
          do j=3,nsh+2
              var66=var66+sh(j)**2
              var77=var77+sh1(j)**2
              var88=var88+sh2(j)**2
              cov67=cov67+sh(j)*sh1(j)
              cov68=cov68+sh(j)*sh2(j)
              cov78=cov78+sh1(j)*sh2(j)
          enddo
          vare1=(var66-2*cov68+var88)*sherr**2
          vare2=2*var77*sherr**2
          cov=sqrt(2.)*(cov67-cov78)*sherr**2
          e1sig=sqrt(vare1)
          e2sig=sqrt(vare2)
          nread=nread+1
          write(*,1002) x,y,f,fe,fw,a,b,th,e1,e2,e1sig,e2sig,cov/e1sig/e2sig,pg11,pg12,pg21,pg22,id,iflg
 1002     format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,4g12.4,f8.4,4g12.4,i6,i5)
      enddo
      write(0,*) 'Sh2p stopped after ',mmcat,' sources'
   99 continue
      write(0,*) 'Ps and es calculated for ',nread,' sources'

      end
      
