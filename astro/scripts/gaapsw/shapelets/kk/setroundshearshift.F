
c pmat is the PSF matrix
c round, sheared1, sheared2 become a set of round objects (one per order)
c and shear operators acting on round objects, then convolved with PSF


      subroutine setroundshearshift_012_o4(pmat,shearmat4,shifted1,shifted2)

      parameter (mmo4=8,msho4=(mmo4+1)*(mmo4+2)/2+2,mm012=5*(mm/2)-1)
      real shearmat4(mm012,15,0:mm/2),shearedcirc4(msho4-2,15,0:mmo4/2),
     :    shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2)
      real sh(msh),sh2(msh),sh3(msh),pmat(msh-2,msh-2),cmat(msh-2,msh-2,msh-2)
      double precision pp

      common /convcoef/ cmat,b1,b2,bout,norder

c Read in the 4th-order coefficients, shearedcirc4(ish,g-order,ic)
c generated automatically with writehighordershearcoefs


#include "shearedcirc4.inc"


c take beta and order from what is specified in common for cmat
      beta=b2
      n=norder
      nsh=(n+1)*(n+2)/2
c take polar order to one level lower than shapelets supplied
      npolar=n-2
c      write(0,*) 'setroundshearshift: ',n,nsh
c catch case where shapelet order is higher than coefs given above
      if (nsh.gt.msho4-2) then
          write(0,*) '4th-ORDER SHEAR FITS ONLY IMPLEMENTED TO N=',mmo4
          stop
      endif

c now make series of round shapes at each even order
      sh(1)=beta
      sh(2)=n
      sh3(1)=beta
      sh3(2)=n
      do m=0,n,2
          do k=3,msh
              sh(k)=0
          enddo
c read off the round shape at order m
          do k=1,nsh
              sh(k+2)=shearedcirc4(k,1,m/2)/sqrt(beta)
          enddo
c act on round object with trans1, then convolve with PSF
          call trans1(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,shifted1(1,m/2),npolar)
c act on round object with trans2, then convolve with PSF
          call trans2(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,shifted2(1,m/2),npolar)
c convolve shape itself with the PSF, write to shearedcirc4
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,shearmat4(1,1,m/2),npolar)
c make the PSF-convolved higher-order terms (g-order kg= 2..15)
          do kg=2,15
              do k=1,nsh
                  sh(k+2)=shearedcirc4(k,kg,m/2)/sqrt(beta)
              enddo
              do k3=1,nsh
                  pp=0
                  do k2=1,nsh
                      pp=pp+sh(k2+2)*pmat(k3,k2)
                  enddo
                  sh3(k3+2)=pp
              enddo
              call to012(sh3,shearmat4(1,kg,m/2),npolar)
          enddo
c          write(0,*) 'SHEARMAT FOR CIRC ORDER',m
c          do k=1,nsh
c              sh(k+2)=shearmat4(k,1,m/2)
c          enddo
c          call writeshape('00',sh)
c          do k=1,nsh
c              sh(k+2)=shearmat4(k,2,m/2)
c          enddo
c          call writeshape('10',sh)
c          do k=1,nsh
c              sh(k+2)=shearmat4(k,3,m/2)
c          enddo
c          call writeshape('01',sh)
      enddo


      return
      end




      subroutine setroundshearshift_o4(pmat,shearmat4,shifted1,shifted2)

      parameter (mmo4=8,msho4=(mmo4+1)*(mmo4+2)/2+2)
      real shearmat4(msh-2,15,0:mm/2),shearedcirc4(msho4-2,15,0:mmo4/2),
     :    shifted1(msh-2,0:mm/2),shifted2(msh-2,0:mm/2)
      real sh(msh),sh2(msh),pmat(msh-2,msh-2),cmat(msh-2,msh-2,msh-2)
      double precision pp

      common /convcoef/ cmat,b1,b2,bout,norder

c Read in the 4th-order coefficients, shearedcirc4(ish,g-order,ic)
c generated automatically with writehighordershearcoefs


#include "shearedcirc4.inc"


c take beta and order from what is specified in common for cmat
      beta=b2
      n=norder
      nsh=(n+1)*(n+2)/2
c      write(0,*) 'setroundshearshift: ',n,nsh
c catch case where shapelet order is higher than coefs given above
      if (nsh.gt.msho4-2) then
          write(0,*) '4th-ORDER SHEAR FITS ONLY IMPLEMENTED TO N=',mmo4
          stop
      endif

c now make series of round shapes at each even order
      sh(1)=beta
      sh(2)=n
      do m=0,n,2
          do k=3,msh
              sh(k)=0
          enddo
c read off the round shape at order m
          do k=1,nsh
              sh(k+2)=shearedcirc4(k,1,m/2)/sqrt(beta)
          enddo
c act on round object with trans1, then convolve with PSF
          call trans1(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              shifted1(k3,m/2)=pp
          enddo
c act on round object with trans2, then convolve with PSF
          call trans2(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              shifted2(k3,m/2)=pp
          enddo
c convolve shape itself with the PSF, write to shearedcirc4
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh(k2+2)*pmat(k3,k2)
              enddo
              shearmat4(k3,1,m/2)=pp
          enddo
c make the PSF-convolved higher-order terms (g-order kg= 2..15)
          do kg=2,15
              do k=1,nsh
                  sh(k+2)=shearedcirc4(k,kg,m/2)/sqrt(beta)
              enddo
              do k3=1,nsh
                  pp=0
                  do k2=1,nsh
                      pp=pp+sh(k2+2)*pmat(k3,k2)
                  enddo
                  shearmat4(k3,kg,m/2)=pp
              enddo
          enddo
c          write(0,*) 'SHEARMAT FOR CIRC ORDER',m
          do k=1,nsh
              sh(k+2)=shearmat4(k,1,m/2)
          enddo
c          call writeshape('00',sh)
          do k=1,nsh
              sh(k+2)=shearmat4(k,2,m/2)
          enddo
c          call writeshape('10',sh)
          do k=1,nsh
              sh(k+2)=shearmat4(k,3,m/2)
          enddo
c          call writeshape('01',sh)
      enddo


      return
      end







      subroutine setroundshearshift_012(pmat,round,sheared1,sheared2,shifted1,shifted2)

      parameter (mm012=5*(mm/2)-1)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :    ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2)
      real sh(msh),sh2(msh),sh3(msh),pmat(msh-2,msh-2),cmat(msh-2,msh-2,msh-2)
      double precision pp

      common /convcoef/ cmat,b1,b2,bout,norder

      
c     Set up the basis for the fit of an object as a PSF-convolved,
c     sheared, shifted, round object:
c
c     for each even order, load the `round' array with a PSF-convolved
c     round object of order n; load the `sheared1' array with PSF x
c     shear1 operator acting on round object; and load the `sheared2'
c     array with PSF x shear2 operator acting on round object
c     the shifted1 and shifted2 arrays are translated objects
c
c     the model will then be of the form 
c             PSF x (1+g1 S1 + g2 S2) x arbitrary round shape
c

c take beta and order from what is specified in common for cmat
      beta=b2
      n=norder
      npolar=n-2

      nsh=(n+1)*(n+2)/2
c now make series of round shapes at each even order
      sh(1)=beta
      sh(2)=n
      sh3(1)=beta
      sh3(2)=n
      do m=0,n,2
          do k=3,msh
              sh(k)=0
              sh3(k)=0
          enddo
c make a round shape at order m
          k0=2+m*(m+1)/2 + 1
          cc=1
          sh(k0)=cc
          do kk=2,m,2
             cc=cc*(m/2-kk/2+1)/(kk/2)*sqrt((kk-1)*kk /real((m-kk+2)*(m-kk+1)))
             sh(k0+kk)=cc
          enddo
c normalize the shape
          sum=seriessum(sh)
          do kk=0,m,2
              sh(k0+kk)=sh(k0+kk)/sum
          enddo
c convolve shape with the PSF via matrix
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
c          write(0,*) 'SHEARMAT FOR CIRC ORDER',m
c          call writeshape('00',sh3)
c write object out as a polar m=0,1,2 shapelet up to order npolar
          call to012(sh3,round(1,m/2),npolar)
c act on round object with shear1, then convolve with PSF
          call shear1(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
c          call writeshape('10',sh3)
          call to012(sh3,sheared1(1,m/2),npolar)
c act on round object with shear2, then convolve with PSF
          call shear2(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
c          call writeshape('01',sh3)
          call to012(sh3,sheared2(1,m/2),npolar)
c act on round object with trans1, then convolve with PSF
          call trans1(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,shifted1(1,m/2),npolar)
c act on round object with trans2, then convolve with PSF
          call trans2(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,shifted2(1,m/2),npolar)
      enddo
      
      return
      end



      subroutine setroundshearshift_012_o2(pmat,round,sheared1,sheared2,shifted1,shifted2,sheared11,sheared12,sheared22)

      parameter (mm012=5*(mm/2)-1)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :    ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2),
     :     sheared11(mm012,0:mm/2),sheared12(mm012,0:mm/2),sheared22(mm012,0:mm/2)
      real sh(msh),sh2(msh),sh3(msh),sh4(msh),pmat(msh-2,msh-2),cmat(msh-2,msh-2,msh-2)
      double precision pp

      common /convcoef/ cmat,b1,b2,bout,norder

      
c     Set up the basis for the fit of an object as a PSF-convolved,
c     sheared, shifted, round object:
c
c     for each even order, load the `round' array with a PSF-convolved
c     round object of order n; load the `sheared1' array with PSF x
c     shear1 operator acting on round object; and load the `sheared2'
c     array with PSF x shear2 operator acting on round object
c     the shifted1 and shifted2 arrays are translated objects
c     sheared11 etc are the 2nd-order operators
c
c     the model will then be of the form 
c             PSF x (1+g1 S1 + g2 S2 + g1^2/2 S11 + g1g2 S1S2 + g2^2/2 S22 + dx T1 + dy T2) 
c                         x arbitrary round shape
c

c take beta and order from what is specified in common for cmat
      beta=b2
      n=norder
      npolar=n-2

      nsh=(n+1)*(n+2)/2
c now make series of round shapes at each even order
      sh(1)=beta
      sh(2)=n
      sh3(1)=beta
      sh3(2)=n
      do m=0,n,2
          do k=3,msh
              sh(k)=0
              sh3(k)=0
          enddo
c make a round shape at order m
          k0=2+m*(m+1)/2 + 1
          cc=1
          sh(k0)=cc
          do kk=2,m,2
             cc=cc*(m/2-kk/2+1)/(kk/2)*sqrt((kk-1)*kk /real((m-kk+2)*(m-kk+1)))
             sh(k0+kk)=cc
          enddo
c normalize the shape *** COULD BE DONE ANALYTICALLY FASTER!! ***
          sum=seriessum(sh)
          do kk=0,m,2
              sh(k0+kk)=sh(k0+kk)/sum
          enddo
c convolve shape with the PSF via matrix
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
c write object out as a polar m=0,1,2 shapelet up to order npolar
          call to012(sh3,round(1,m/2),npolar)
c act on round object with shear1, then convolve with PSF
          call shear1(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared1(1,m/2),npolar)
c similarly construct sheared11 
          call shear1(sh2,sh4)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh4(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared11(1,m/2),npolar)
c and sheared12
          call shear2(sh2,sh4)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh4(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared12(1,m/2),npolar)
c act on round object with shear2, then convolve with PSF
          call shear2(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared2(1,m/2),npolar)
c similarly construct sheared22
          call shear2(sh2,sh4)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh4(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared22(1,m/2),npolar)
c act on round object with trans1, then convolve with PSF
          call trans1(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,shifted1(1,m/2),npolar)
c act on round object with trans2, then convolve with PSF
          call trans2(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,shifted2(1,m/2),npolar)
      enddo
      
      if (.false.) then
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Round ',2*i
              write(0,*) 'm=0:',(round(k,i),k=1,norder/2)
              write(0,*) 'cos:',(round(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(round(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(round(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(round(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared1 ',2*i
              write(0,*) 'm=0:',(sheared1(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared1(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared1(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared1(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared1(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared2 ',2*i
              write(0,*) 'm=0:',(sheared2(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared2(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared2(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared2(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared2(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared11 ',2*i
              write(0,*) 'm=0:',(sheared11(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared11(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared11(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared11(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared11(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared12 ',2*i
              write(0,*) 'm=0:',(sheared12(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared12(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared12(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared12(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared12(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared12 ',2*i
              write(0,*) 'm=0:',(sheared22(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared22(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared22(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared22(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared22(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Shifted1 ',2*i
              write(0,*) 'm=0:',(shifted1(k,i),k=1,norder/2)
              write(0,*) 'cos:',(shifted1(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(shifted1(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(shifted1(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(shifted1(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Shifted2 ',2*i
              write(0,*) 'm=0:',(shifted2(k,i),k=1,norder/2)
              write(0,*) 'cos:',(shifted2(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(shifted2(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(shifted2(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(shifted2(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
       endif
      return
      end


      subroutine setroundshearshift_012_o2_nopsf(beta,n,round,sheared1,sheared2,shifted1,shifted2,sheared11,sheared12,sheared22)

      parameter (mm012=5*(mm/2)-1)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :    ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2),
     :     sheared11(mm012,0:mm/2),sheared12(mm012,0:mm/2),sheared22(mm012,0:mm/2)
      real sh(msh),sh2(msh),sh3(msh),sh4(msh)

      
c     Set up the basis for the fit of an object as a PSF-convolved,
c     sheared, shifted, round object:
c
c     for each even order, load the `round' array with a PSF-convolved
c     round object of order n; load the `sheared1' array with PSF x
c     shear1 operator acting on round object; and load the `sheared2'
c     array with PSF x shear2 operator acting on round object
c     the shifted1 and shifted2 arrays are translated objects
c     sheared11 etc are the 2nd-order operators
c
c     the model will then be of the form 
c             PSF x (1+g1 S1 + g2 S2 + g1^2/2 S11 + g1g2 S1S2 + g2^2/2 S22 + dx T1 + dy T2) 
c                         x arbitrary round shape
c

c take beta and order from what is specified in common for cmat
      npolar=n-2

      nsh=(n+1)*(n+2)/2
c now make series of round shapes at each even order
      sh(1)=beta
      sh(2)=n
      sh3(1)=beta
      sh3(2)=n
      do m=0,n,2
          do k=3,msh
              sh(k)=0
              sh3(k)=0
          enddo
c make a round shape at order m
          k0=2+m*(m+1)/2 + 1
          cc=1
          sh(k0)=cc
          do kk=2,m,2
             cc=cc*(m/2-kk/2+1)/(kk/2)*sqrt((kk-1)*kk /real((m-kk+2)*(m-kk+1)))
             sh(k0+kk)=cc
          enddo
c normalize the shape *** COULD BE DONE ANALYTICALLY FASTER!! ***
          sum=seriessum(sh)
          do kk=0,m,2
              sh(k0+kk)=sh(k0+kk)/sum
          enddo
c write object out as a polar m=0,1,2 shapelet up to order npolar
          call to012(sh,round(1,m/2),npolar)
c act on round object with shear1, then convolve with PSF
          call shear1(sh,sh2)
          call to012(sh2,sheared1(1,m/2),npolar)
c similarly construct sheared11 
c          call shear1(sh2,sh4)
          call shear11(sh,sh4)
          call to012(sh4,sheared11(1,m/2),npolar)
c and sheared12
          call shear2(sh2,sh4)
          call to012(sh4,sheared12(1,m/2),npolar)
c act on round object with shear2, then convolve with PSF
          call shear2(sh,sh2)
          call to012(sh2,sheared2(1,m/2),npolar)
c similarly construct sheared22
          call shear2(sh2,sh4)
          call to012(sh4,sheared22(1,m/2),npolar)
c act on round object with trans1, then convolve with PSF
          call trans1(sh,sh2)
          call to012(sh2,shifted1(1,m/2),npolar)
c act on round object with trans2, then convolve with PSF
          call trans2(sh,sh2)
          call to012(sh2,shifted2(1,m/2),npolar)
      enddo
      
      if (.false.) then
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Round ',2*i
              write(0,*) 'm=0:',(round(k,i),k=1,norder/2)
              write(0,*) 'cos:',(round(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(round(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(round(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(round(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared1 ',2*i
              write(0,*) 'm=0:',(sheared1(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared1(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared1(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared1(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared1(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared2 ',2*i
              write(0,*) 'm=0:',(sheared2(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared2(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared2(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared2(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared2(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared11 ',2*i
              write(0,*) 'm=0:',(sheared11(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared11(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared11(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared11(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared11(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared12 ',2*i
              write(0,*) 'm=0:',(sheared12(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared12(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared12(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared12(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared12(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared12 ',2*i
              write(0,*) 'm=0:',(sheared22(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared22(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared22(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared22(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared22(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Shifted1 ',2*i
              write(0,*) 'm=0:',(shifted1(k,i),k=1,norder/2)
              write(0,*) 'cos:',(shifted1(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(shifted1(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(shifted1(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(shifted1(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Shifted2 ',2*i
              write(0,*) 'm=0:',(shifted2(k,i),k=1,norder/2)
              write(0,*) 'cos:',(shifted2(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(shifted2(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(shifted2(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(shifted2(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
       endif
      return
      end




      subroutine setroundshearshift_012_o3(pmat,round,sheared1,sheared2,shifted1,shifted2,sheared11,sheared12,sheared22,sheared111,sheared112,sheared122,sheared222)

      parameter (mm012=5*(mm/2)-1)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :    ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2),
     :     sheared11(mm012,0:mm/2),sheared12(mm012,0:mm/2),
     :     sheared22(mm012,0:mm/2),
     :     sheared111(mm012,0:mm/2),sheared112(mm012,0:mm/2),
     :     sheared122(mm012,0:mm/2),sheared222(mm012,0:mm/2)

      real sh(msh),sh2(msh),sh3(msh),sh4(msh),sh5(msh),
     :     pmat(msh-2,msh-2),cmat(msh-2,msh-2,msh-2)
      double precision pp

      common /convcoef/ cmat,b1,b2,bout,norder

      
c     Set up the basis for the fit of an object as a PSF-convolved,
c     sheared, shifted, round object:
c
c     for each even order, load the `round' array with a PSF-convolved
c     round object of order n; load the `sheared1' array with PSF x
c     shear1 operator acting on round object; and load the `sheared2'
c     array with PSF x shear2 operator acting on round object
c     the shifted1 and shifted2 arrays are translated objects
c     sheared11 etc are the 2nd-order operators
c
c     the model will then be of the form 
c             PSF x (1+g1 S1 + g2 S2 + g1^2/2 S11 + g1g2 S1S2 + g2^2/2 S22 + dx T1 + dy T2) 
c                         x arbitrary round shape
c

c take beta and order from what is specified in common for cmat
      beta=b2
      n=norder
      npolar=n-2

      nsh=(n+1)*(n+2)/2
c now make series of round shapes at each even order
      sh(1)=beta
      sh(2)=n
      sh3(1)=beta
      sh3(2)=n
      do m=0,n,2
          do k=3,msh
              sh(k)=0
              sh3(k)=0
          enddo
c make a round shape at order m
          k0=2+m*(m+1)/2 + 1
          cc=1
          sh(k0)=cc
          do kk=2,m,2
             cc=cc*(m/2-kk/2+1)/(kk/2)*sqrt((kk-1)*kk /real((m-kk+2)*(m-kk+1)))
             sh(k0+kk)=cc
          enddo
c normalize the shape *** COULD BE DONE ANALYTICALLY FASTER!! ***
          sum=seriessum(sh)
          do kk=0,m,2
              sh(k0+kk)=sh(k0+kk)/sum
          enddo
c convolve shape with the PSF via matrix
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
c write object out as a polar m=0,1,2 shapelet up to order npolar
          call to012(sh3,round(1,m/2),npolar)
c act on round object with shear1, then convolve with PSF
          call shear1(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared1(1,m/2),npolar)
c similarly construct sheared11 
          call shear1(sh2,sh4)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh4(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared11(1,m/2),npolar)
c construct sheared111
          call shear1(sh4,sh5)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh5(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared111(1,m/2),npolar)
c construct sheared112
          call shear2(sh4,sh5)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh5(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared112(1,m/2),npolar)
c and sheared12
          call shear2(sh2,sh4)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh4(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared12(1,m/2),npolar)
c act on round object with shear2, then convolve with PSF
          call shear2(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared2(1,m/2),npolar)
c similarly construct sheared22
          call shear2(sh2,sh4)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh4(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared22(1,m/2),npolar)
c construct sheared122 
          call shear1(sh4,sh5)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh5(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared122(1,m/2),npolar)
c construct sheared222
          call shear2(sh4,sh5)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh4(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,sheared222(1,m/2),npolar)
c act on round object with trans1, then convolve with PSF
          call trans1(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,shifted1(1,m/2),npolar)
c act on round object with trans2, then convolve with PSF
          call trans2(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              sh3(k3+2)=pp
          enddo
          call to012(sh3,shifted2(1,m/2),npolar)
      enddo
      
      if (.false.) then
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Round ',2*i
              write(0,*) 'm=0:',(round(k,i),k=1,norder/2)
              write(0,*) 'cos:',(round(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(round(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(round(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(round(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared1 ',2*i
              write(0,*) 'm=0:',(sheared1(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared1(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared1(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared1(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared1(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared2 ',2*i
              write(0,*) 'm=0:',(sheared2(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared2(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared2(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared2(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared2(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared11 ',2*i
              write(0,*) 'm=0:',(sheared11(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared11(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared11(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared11(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared11(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared12 ',2*i
              write(0,*) 'm=0:',(sheared12(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared12(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared12(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared12(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared12(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared12 ',2*i
              write(0,*) 'm=0:',(sheared22(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared22(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared22(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared22(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared22(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Shifted1 ',2*i
              write(0,*) 'm=0:',(shifted1(k,i),k=1,norder/2)
              write(0,*) 'cos:',(shifted1(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(shifted1(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(shifted1(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(shifted1(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Shifted2 ',2*i
              write(0,*) 'm=0:',(shifted2(k,i),k=1,norder/2)
              write(0,*) 'cos:',(shifted2(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(shifted2(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(shifted2(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(shifted2(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
       endif
      return
      end



      subroutine setroundshearshift_012_nopsf(beta,n,round,sheared1,sheared2,shifted1,shifted2)

      parameter (mm012=5*(mm/2)-1)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :    ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2)
      real sh(msh),sh2(msh)

c AS SETROUNDSHEARSHIFT_012, BUT WITHOUT PSF CONVOLUTION

c     Set up the basis for the fit of an object as a PSF-convolved,
c     sheared, shifted, round object:
c
c     for each even order, load the `round' array with a PSF-convolved
c     round object of order n; load the `sheared1' array with PSF x
c     shear1 operator acting on round object; and load the `sheared2'
c     array with PSF x shear2 operator acting on round object
c     the shifted1 and shifted2 arrays are translated objects
c
c     the model will then be of the form 
c             PSF x (1+g1 S1 + g2 S2) x arbitrary round shape
c



      npolar=n-2

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
c normalize the shape
          sum=seriessum(sh)
          do kk=0,m,2
              sh(k0+kk)=sh(k0+kk)/sum
          enddo
c write object out as a polar m=0,1,2 shapelet up to order npolar
          call to012(sh,round(1,m/2),npolar)
c act on round object with shear1, don't convolve with PSF
          call shear1(sh,sh2)
          call to012(sh2,sheared1(1,m/2),npolar)
c act on round object with shear2, don't convolve with PSF
          call shear2(sh,sh2)
          call to012(sh2,sheared2(1,m/2),npolar)
c act on round object with trans1, don't convolve with PSF
          call trans1(sh,sh2)
          call to012(sh2,shifted1(1,m/2),npolar)
c act on round object with trans2, don't convolve with PSF
          call trans2(sh,sh2)
          call to012(sh2,shifted2(1,m/2),npolar)
      enddo
      
      return
      end


      subroutine setroundshearshift(pmat,round,sheared1,sheared2,shifted1,shifted2)

      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
     :    ,shifted1(msh-2,0:mm/2),shifted2(msh-2,0:mm/2)
      real sh(msh),sh2(msh),pmat(msh-2,msh-2),cmat(msh-2,msh-2,msh-2)
      double precision pp

      common /convcoef/ cmat,b1,b2,bout,norder

      
c     Set up the basis for the fit of an object as a PSF-convolved,
c     sheared, shifted, round object:
c
c     for each even order, load the `round' array with a PSF-convolved
c     round object of order n; load the `sheared1' array with PSF x
c     shear1 operator acting on round object; and load the `sheared2'
c     array with PSF x shear2 operator acting on round object
c     the shifted1 and shifted2 arrays are translated objects
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
c normalize the shape
          sum=seriessum(sh)
          do kk=0,m,2
              sh(k0+kk)=sh(k0+kk)/sum
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
c act on round object with trans1, then convolve with PSF
          call trans1(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              shifted1(k3,m/2)=pp
          enddo
c act on round object with trans2, then convolve with PSF
          call trans2(sh,sh2)
          do k3=1,nsh
              pp=0
              do k2=1,nsh
                  pp=pp+sh2(k2+2)*pmat(k3,k2)
              enddo
              shifted2(k3,m/2)=pp
          enddo
      enddo
      
      return
      end


      subroutine setroundshearshiftnopsf(beta,n,round,sheared1,sheared2,shifted1,shifted2)

      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
     :    ,shifted1(msh-2,0:mm/2),shifted2(msh-2,0:mm/2)
      real sh(msh),sh2(msh)

      
c     Set up the basis for the fit of an object as a 
c     sheared, shifted, round object:
c
c     for each even order, load the `round' array with a
c     round object of order n; load the `sheared1' array with
c     shear1 operator acting on round object; and load the `sheared2'
c     array with shear2 operator acting on round object
c     the shifted1 and shifted2 arrays are translated objects
c
c     the model will then be of the form 
c             (1+g1 S1 + g2 S2) x arbitrary round shape
c


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
c convolve shape with the PSF via matrix  --- NO PSF
          do k3=1,nsh
              round(k3,m/2)=sh(k3+2)
          enddo
c act on round object with shear1, then convolve with PSF  --- NO PSF
          call shear1(sh,sh2)
          do k3=1,nsh
              sheared1(k3,m/2)=sh2(k3+2)
          enddo
c act on round object with shear2, then convolve with PSF  --- NO PSF
          call shear2(sh,sh2)
          do k3=1,nsh
              sheared2(k3,m/2)=sh2(k3+2)
          enddo
c act on round object with trans1, then convolve with PSF  --- NO PSF
          call trans1(sh,sh2)
          do k3=1,nsh
              shifted1(k3,m/2)=sh2(k3+2)
          enddo
c act on round object with trans2, then convolve with PSF  --- NO PSF
          call trans2(sh,sh2)
          do k3=1,nsh
              shifted2(k3,m/2)=sh2(k3+2)
          enddo
      enddo
      
      return
      end
