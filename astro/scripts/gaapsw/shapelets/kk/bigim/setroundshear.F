
c pmat is the PSF matrix
c round, sheared1, sheared2 become a set of round objects (one per order)
c and shear operators acting on round objects, then convolved with PSF

      subroutine setroundshear(pmat,round,sheared1,sheared2)

      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
      real sh(msh),sh2(msh),pmat(msh-2,msh-2),cmat(msh-2,msh-2,msh-2)
      double precision pp

      common /convcoef/ cmat,b1,b2,bout,norder

      
c     Set up the basis for the fit of an object as a PSF-convolved,
c     sheared, round object:
c
c     for each even order, load the `round' array with a PSF-convolved
c     round object of order n; load the `sheared1' array with PSF x
c     shear1 operator acting on round object; and load the `sheared2'
c     array with PSF x shear2 operator acting on round object
c
c     the model will then be of the form 
c             PSF x (1+g1 S1 + g2 S2) x arbitrary round shape
c

c take beta and order from what is specified in common for cmat
      beta=b2
      n=norder
      nsh=(n+1)*(n+2)/2

c now make series of round shapes at each even order
      sh(1)=beta
      sh(2)=n
      do m=0,n,2
          do k=3,msh
              sh(k)=0
          enddo
c make a round shape at order m
          k0=2+m*(m+1)/2 + 1
          cc=1
          sh(k0)=cc
          do kk=2,m,2
             cc=cc*(m/2-kk/2+1)/(kk/2)*sqrt((kk-1)*kk /real((m-kk+2)*(m-kk+1)))
             sh(k0+kk)=cc
          enddo
c convolve shape with the PSF via matrix
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh(k2+2)*pmat(k3,k2)
              enddo
              round(k3,m/2)=pp
          enddo
c act on round object with shear1, then convolve with PSF
          call shear1(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sheared1(k3,m/2)=pp
          enddo
c act on round object with shear2, then convolve with PSF
          call shear2(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sheared2(k3,m/2)=pp
          enddo
      enddo
      
      return
      end


      subroutine setroundshearnopsf(beta,n,round,sheared1,sheared2)

      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
      real sh(msh),sh2(msh),pmat(msh-2,msh-2),cmat(msh-2,msh-2,msh-2)

      
c     Set up the basis for the fit of an object as a PSF-convolved,
c     sheared, round object:
c
c     for each even order, load the `round' array with a PSF-convolved
c     round object of order n; load the `sheared1' array with PSF x
c     shear1 operator acting on round object; and load the `sheared2'
c     array with PSF x shear2 operator acting on round object
c
c     the model will then be of the form 
c             PSF x (1+g1 S1 + g2 S2) x arbitrary round shape
c

c take beta and order from what is specified in common for cmat
      nsh=(n+1)*(n+2)/2

c now make series of round shapes at each even order
      sh(1)=beta
      sh(2)=n
      do m=0,n,2
          do k=3,msh
              sh(k)=0
          enddo
c make a round shape at order m
          k0=2+m*(m+1)/2 + 1
          cc=1
          sh(k0)=cc
          do kk=2,m,2
             cc=cc*(m/2-kk/2+1)/(kk/2)*sqrt((kk-1)*kk /real((m-kk+2)*(m-kk+1)))
             sh(k0+kk)=cc
          enddo
c convolve shape with the PSF via matrix --- NO PSF
          do k3=1,nsh
              round(k3,m/2)=sh(k3+2)
          enddo
c act on round object with shear1, then convolve with PSF --- NO PSF
          call shear1(sh,sh2)
          do k3=1,nsh
              sheared1(k3,m/2)=sh2(k3+2)
          enddo
c act on round object with shear2, then convolve with PSF --- NO PSF
          call shear2(sh,sh2)
          do k3=1,nsh
              sheared2(k3,m/2)=sh2(k3+2)
          enddo
      enddo
      
      return
      end
