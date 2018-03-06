c read a shape catalogue (produced with cat2sh), and a psf map (in psf.map)
c call fitshear for each shape to deduce best-fit g1, g2 and cov matrix.

      parameter (mm012=5*(mm/2)-1)
      real sh(msh),shpsf(msh),psfcoef(mmpoly2,msh-2),psfmat(msh-2,msh-2)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :    ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2)
      real coefs(mm/2+4),covar(mm/2+4,mm/2+4)
      real p(0:mm),frac(2:6)
      common /shearshiftfns012/ round,sheared1,sheared2,shifted1,shifted2,nshear

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
          call interpolpsf(x,y,xmax,ymax,psfcoef,npsf,betapsf,nfit,shpsf)
          beta2=sqrt(max(beta,1.1*betapsf)**2-betapsf**2)
          call setpsfmat(shpsf,beta,beta2,psfmat)
          nshear=n
          call setroundshearshift_012(psfmat,round,sheared1,sheared2,shifted1,shifted2)
          call fitshearshift_err_012(sh,sherr,g1,g2,g1sig,g2sig,coefs,covar,.false.)
          nread=nread+1
          write(*,1002) x,y,f,fe,fw,a,b,th,g1,g2,g1sig,g2sig,
     :        max(-2.,min(2.,covar(n/2+1,n/2+2)/g1sig/g2sig)),id,iflg,
     :        (coefs(ic),ic=1,n/2)
 1002     format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,4g12.4,f8.4,i10,i5,9g12.4)
      enddo
      write(0,*) 'Sh2gerrshift012 stopped after ',mmcat,' sources'
   99 continue
      write(0,*) 'Shears calculated for ',nread,' sources'

      end
      
