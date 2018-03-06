      real psf(30,30),shpsf(msh),z(30,30),shape(msh)
      real pmat(msh-2,msh-2)
      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
      real coefs(mm/2+2),covar(mm/2+2,mm/2+2)
      common /shearfns/ round,sheared1,sheared2,nshear
    
      call pgbeg(0,'/xs',4,3)
      do i=1,30
          x=i-15.5
          do j=1,30
              y=j-15.5
c              psf(i,j)=exp(-0.5*(x*x+y*y)/4)
              psf(i,j)=(1+0.1*x*x + 0.1*y**2)**(-2.5)
          enddo
      enddo

      n=10
      betaout=2.
      beta=betaout
c      call fitshapelets(psf,30,30,30,15.5,15.5,10.,0.,shpsf,n,0.)
      call encode(psf,30,30,30,15.5,15.5,10.,0.,shpsf,n,0.)
      betapsf=shpsf(1)
      call plotshape(shpsf,-15,'PSF')
      call writeshapen('PSF',shpsf)
      call setpsfmat(shpsf,betaout,beta,pmat)
      call setroundshear(pmat,round,sheared1,sheared2)
      nshear=n
      
      do kk=1,10
          vxx=3.5+0.1*kk
          vxy=0
          vyy=3.5-0.1*kk
          do i=1,30
              x=i-15.5
              do j=1,30
                  y=j-15.5
                  z(i,j)=exp(-0.5*(x*x/vxx+y*y/vyy))
              enddo
          enddo
          call fitshapelets(z,30,30,30,15.5,15.5,10.,0.,shape,n,betaout)
          call plotshape(shape,15,'SHAPE')
          call fitshear(shape,g1,g2,g1sig,g2sig,coefs,covar,.false.)
          write(0,'(4f12.5)') g1,g2,g1sig,g2sig
      enddo
      call pgend
      end
