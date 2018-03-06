c read a shape catalogue (produced with cat2sh), and a psf map (in psf.map)
c call fitshear for each shape to deduce best-fit g1, g2 and cov matrix.

      real sh(msh),shpsf(msh),psfcoef(mmpoly2,msh-2),psfmat(msh-2,msh-2)
      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
      real coefs(mm/2+2),covar(mm/2+2,mm/2+2)
      common /shearfns/ round,sheared1,sheared2,nshear

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

      nread=0
      do i=1,mmcat
    1     read(*,*,err=1,end=99) x,y,f,fe,fw,a,b,th,id,iflg
          read(*,*) sum,dx,dy,bg,fitradius
          read(*,*) beta,n
          read(*,*) (sh(k),k=3,(n+1)*(n+2)/2+2)
          read(*,*)
          sh(1)=beta
          sh(2)=n
          do k=3,(n+1)*(n+2)/2+2
c clooge: reduce signal in shapelets so that the error bars on shear come
c out bigger.
              sh(k)=sh(k)*sum/100
          enddo
          call interpolpsf(x,y,xmax,ymax,psfcoef,npsf,betapsf,nfit,shpsf)
c          call writeshapen('PSF',shpsf)
c          call writeshapen('Target',sh)
          call setpsfmat(shpsf,beta,beta,psfmat)
          nshear=n
          call setroundshear(psfmat,round,sheared1,sheared2)
          call fitshear(sh,g1,g2,g1sig,g2sig,coefs,covar,.false.)
          nread=nread+1
          write(*,1002) x,y,f,fe,fw,a,b,th,g1,g2,g1sig,g2sig,covar(n/2+1,n/2+2)/g1sig/g2sig,id,iflg
 1002     format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,4g12.4,f8.4,i6,i5)
      enddo
      write(0,*) 'Sh2g stopped after ',mmcat,' sources'
   99 continue
      write(0,*) 'Shears calculated for ',nread,' sources'

      end
      
