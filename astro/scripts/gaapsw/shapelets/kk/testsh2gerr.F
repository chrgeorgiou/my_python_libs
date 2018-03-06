c read a shape catalogue (produced with cat2sh), and a psf map (in psf.map)
c call fitshear for each shape to deduce best-fit g1, g2 and cov matrix.

      parameter (mm012=5*(mm/2)-1)
      real sh(msh),shpsf(msh),psfcoef(mmpoly2,msh-2),psfmat(msh-2,msh-2)
      real shearmat4(msh-2,15,0:mm/2),
     :     shifted1(msh-2,0:mm/2),shifted2(msh-2,0:mm/2)
      real coefs(mm/2+4),covar(mm/2+4,mm/2+4)
      real p(0:mm),frac(2:6)
      common /shearshiftfns_o4/ shearmat4,shifted1,shifted2,nshear

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
c write out a perfect PSF convolved sheared circular shapelet as a test
      do i=1,1
    1     read(*,*,err=1,end=99) x,y,f,fe,fw,a,b,th,id,iflg
          read(*,*) sum,dx,dy,bg,fitradius,sherr
          read(*,*) beta,n
          read(*,*) (sh(k),k=3,(n+1)*(n+2)/2+2)
          read(*,*)
          sh(1)=beta
          sh(2)=n
          nsh=(n+1)*(n+2)/2
          call interpolpsf(x,y,xmax,ymax,psfcoef,npsf,betapsf,nfit,shpsf)
          beta2=sqrt(max(beta,1.1*betapsf)**2-betapsf**2)
          call setpsfmat(shpsf,beta,beta2,psfmat)
          nshear=n
          call setroundshearshift_o4(psfmat,shearmat4,shifted1,shifted2)
          do k=1,nsh
             sh(k+2)=shearmat4(k,1,1)+0.07*shearmat4(k,2,1)-0.03*shearmat4(k,3,1)
          enddo
          write(*,*) x,y,f,fe,fw,a,b,th,id,iflg
          write(*,*) sum,dx,dy,bg,fitradius,sherr
          write(*,*) beta,n
          write(*,*) (sh(k),k=3,(n+1)*(n+2)/2+2)
          write(*,*)
      enddo
   99 continue
      end
      
