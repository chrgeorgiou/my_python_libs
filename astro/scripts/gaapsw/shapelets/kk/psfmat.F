      real zpsf(50,50),shpsf(msh),cmat(msh-2,msh-2,msh-2)
      real sh(msh),round(msh-2,0:mm),sh2(msh),sh3(msh),pmat(msh-2,msh-2)
      common /convcoef/ cmat,b1,b2,bout,norder

c get PSF

      do i=1,31
          x=i-16
          do j=1,31
              y=j-16
              zpsf(i,j)=(1+0.2*x**2+0.3*y**2)**(-2.)
          enddo
      enddo
      npsf=10
      betapsf=1.5
      call encode(zpsf,50,31,31,16.,16.,10.,0.,shpsf,npsf,betapsf)
      
      
      call pgbeg(0,'/xs',4,3)
      call plotshape(shpsf,16,'PSF')
      call writeshape('PSF',shpsf)
      
      beta=3.
      n=shpsf(2)
      betaout=beta

c skip calculation of the coeffs if same as what is stored in common block
      if (b1.eq.betapsf .and. b2.eq.beta .and. bout.eq.betaout 
     :    .and. norder.eq.n) goto 99
       call setcmatrix(cmat,beta1,beta2,betaout,n)
c update info in common block on which parameters are in there     
      b1=betapsf
      b2=beta
      bout=betaout
      norder=n
c make the PSF matrix by contracting with cmat
      do k3=1,msh-2
          do k2=1,msh-2
              pp=0
              do k1=1,msh-2
                  pp=pp+cmat(k3,k1,k2)*shpsf(k1+2)
              enddo
              pmat(k3,k2)=pp
          enddo
      enddo

   99 continue
      sh(1)=beta
      sh(2)=n
      k=2
      do m=0,n,2
          do k=3,msh
              sh(k)=0
          enddo
          k0=2+m*(m+1)/2 + 1
          cc=1
          sh(k0)=cc
          do kk=2,m,2
             cc=cc*(m/2-kk/2+1)/(kk/2)*sqrt((kk-1)*kk /real((m-kk+2)*(m-kk+1)))
             sh(k0+kk)=cc
          enddo
          call plotshape(sh,31,'Unconvolved')
          call writeshapen('unconvolved ',sh)
          call convolve(shpsf,sh,betaout,sh2)
          call plotshape(sh2,31,'Convolved')
          call writeshapen('convolved ',sh2)
          sh3(1)=betaout
          sh3(2)=n
          do k3=1,msh-2
              pp=0
              do k2=1,msh-2
                  pp=pp+sh(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call plotshape(sh3,31,'Times Pmat')
          call writeshapen('mult by pmat ',sh3)
          write(0,*)
          do k=1,(n+1)*(n+2)/2
              round(k,m/2)=sh2(k+2)
          enddo
      enddo
      




      call pgend
      end
