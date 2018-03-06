      real zpsf(50,50),shpsf(msh),cmat(msh-2,msh-2,msh-2)
      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
      real sh(msh),sh2(msh),sh3(msh),pmat(msh-2,msh-2)
      double precision pp
      common /convcoef/ cmat,b1,b2,bout,norder

c get PSF

      do i=1,31
          x=i-16
          do j=1,31
              y=j-16
              r2=x*x*0.5+y*y
              r=sqrt(r2)+0.5*x
              zpsf(i,j)=(2/(exp(-r/4)+exp(r/4)))**2
          enddo
      enddo

      npsf=6
      betapsf=2.5
      call fitshapelets(zpsf,50,31,31,16.,16.,10.,0.,shpsf,npsf,betapsf)
      
      
      call pgbeg(0,'/xs',4,4)
      call plotshape(shpsf,41,'PSF')
      call writeshape('PSF',shpsf)
      
      beta=4
      n=shpsf(2)
      betaout=beta*1.2

      call setpsfmat(shpsf,betaout,beta,pmat)

      do j=1,8
          write(0,'(8f10.5)') (pmat(i,j),i=1,8)
      enddo

c do a test where convolve a shape with PSF directly, and with PSF matrix

      do i=1,31
          x=i-16
          do j=1,31
              y=j-16
              r2=x*x*0.5+y*y
              r=sqrt(r2)
              zpsf(i,j)=(2/(exp(-r/4)+exp(r/4)))**2
          enddo
      enddo
      call fitshapelets(zpsf,50,31,31,16.,16.,10.,0.,sh,n,beta)
      call plotshape(sh,41,'SH')
      call writeshape('SH',sh)

      sh2(1)=betaout
      sh2(2)=n
      do i=1,(n+1)*(n+2)/2
          sh2(i+2)=0
          do j=1,(n+1)*(n+2)/2
              sh2(i+2)=sh2(i+2)+pmat(i,j)*sh(j+2)
          enddo
      enddo
      call plotshape(sh2,41,'conv w pmat')
      call writeshapen('CONV W PSFMAT',sh2)

      call convolve(shpsf,sh,betaout,sh3)
      call plotshape(sh3,41,'conv w sum')
      call writeshapen('CONV w SUM',sh3)

      call setroundshear(pmat,round,sheared1,sheared2)

c display the constructed shapes
      sh(1)=betaout
      sh(2)=n
      nsh=(n+1)*(n+2)/2
      do m=0,n,2
          do k=1,nsh
              sh(k+2)=round(k,m/2)
          enddo
          call plotshape(sh,41,'Round')
          call writeshapen('Round',sh)
          do k=1,nsh
              sh(k+2)=sheared1(k,m/2)
          enddo
          call plotshape(sh,41,'Sheared 1')
          call writeshapen('Sheared 1',sh)
          do k=1,nsh
              sh(k+2)=sheared2(k,m/2)
          enddo
          call plotshape(sh,41,'Sheared 2')
          call writeshapen('Sheared 2',sh)
      enddo

      call pgend
      end
      
