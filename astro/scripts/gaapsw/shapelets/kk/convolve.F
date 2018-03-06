c
c convolution ignoring pixellation effects.
c
      subroutine convolve(shape1,shape2,betaout,result)
      real shape1(msh),shape2(msh),result(msh)
      real cmat(msh-2,msh-2,msh-2)

      common /convcoef/ cmat,b1,b2,bout,norder

      beta1=shape1(1)
      n1=shape1(2)
      beta2=shape2(1)
      n2=shape2(2)

      n=max(n1,n2)
      result(1)=betaout
      result(2)=n

c skip calculation of the coeffs if same as what is stored in common block
      if (b1.eq.beta1 .and. b2.eq.beta2 .and. bout.eq.betaout 
     :    .and. norder.eq.n) goto 99

c
c construct the 1-D C matrix (Refregier et at paper II) 
c 

      call setcmatrix(cmat,beta1,beta2,betaout,n)
c update info in common block on which parameters are in there     
      b1=beta1
      b2=beta2
      bout=betaout
      norder=n

c land here if there was no need to calculate the coefficients
c
c operate the convolution matrix on shape 1 and shape2 to give result

   99 continue

      do k3=1,(n+1)*(n+2)/2
          sum=0
          do k1=1,(n1+1)*(n1+2)/2
              do k2=1,(n2+1)*(n2+2)/2
                  sum=sum+cmat(k3,k1,k2)*shape1(k1+2)*shape2(k2+2)
              enddo
          enddo
          result(k3+2)=sum
      enddo

      return
      end


c
c correlation function 
c similar to convolution except change the sign of the odd orders in shape2
c
      subroutine correlate(shape1,shape2,betaout,result)
      real shape1(msh),shape2(msh),shape3(msh),result(msh)

c change the sign of all odd coefs in shape2
      beta2=shape2(1)
      n2=shape2(2)
      shape3(1)=beta2
      shape3(2)=n2
      kk=3
      s=1
      do k=0,n2
          do l=0,k
              shape3(kk)=shape2(kk)*s
              kk=kk+1
          enddo
          s=-s
      enddo
      call convolve(shape1,shape3,betaout,result)
      return
      end

          





      subroutine setcmatrix(cmat,beta1,beta2,betaout,n)

      real c1d(0:mm,0:mm,0:mm),cmat(msh-2,msh-2,msh-2),el(-1:mm,-1:mm,-1:mm)

      write(0,'('' Calculating Clmn for n='',i2,'' betas'',3f7.3,''... ''$)') 
     :    n,betaout,beta1,beta2
      a1=1/betaout
      a2=1/beta1
      a3=1/beta2
      anu=sqrt(1./(1./a1**2+1./a2**2+1./a3**2))
      a=sqrt(2.)*anu/a1
      b=sqrt(2.)*anu/a2
      c=sqrt(2.)*anu/a3
      a2min1=a*a-1
      b2min1=b*b-1
      c2min1=c*c-1
      ab=a*b
      bc=b*c
      ac=a*c
c first construct the edges (only 1 nonzero cpt) of the el matrix = L(a,b,c)
      el(0,0,0)=1
      do i=0,n
          do j=0,n
              el(-1,i,j)=0
              el(i,-1,j)=0
              el(i,j,-1)=0
          enddo
      enddo
      do isum=2,3*n,2
          do lll=0,min(n,isum)
              do mmm=max(0,isum-lll-n),min(n,isum-lll)
                  nnn=isum-lll-mmm
                  if (nnn.gt.0) then
c apply recipe 3
                      el(lll,mmm,nnn)=2*lll*ac*el(lll-1,mmm,nnn-1)
     :                    +2*mmm*bc*el(lll,mmm-1,nnn-1)
     :                    +2*(nnn-1)*c2min1*el(lll,mmm,nnn-2)
                  elseif (mmm.gt.0) then
c apply recipe 2 with nnn=0
                      el(lll,mmm,0)=2*lll*ab*el(lll-1,mmm-1,0)
     :                    +2*(mmm-1)*b2min1*el(lll,mmm-2,0)
                  else
c apply recipe 1 with mmm=nnn=0
                      el(lll,0,0)=2*(lll-1)*a2min1*el(lll-2,0,0)
                  endif
              enddo
          enddo
      enddo

c
c now turn the el matrix into the 1-D c1d matrices
c

      coefi=1
      do i=0,n
          if (i.gt.0) coefi=coefi*i*2
          coefj=1
          do j=0,n
              if (j.gt.0) coefj=coefj*j*2
              coefk=1
              do k=0,n
                  if (k.gt.0) coefk=coefk*k*2
                  if (mod(i+j+k,2).eq.1) then
                      c1d(i,j,k)=0
                  else
                      c1d(i,j,k)=sqrt(2*3.1415926535)*el(i,j,k)*anu
     :                 /sqrt(coefi*coefj*coefk*0.5*sqrt(3.1415926535)*a1*a2*a3)
                      if (mod((i+j+k)/2,2).eq.1) c1d(i,j,k)=-c1d(i,j,k)
                      if (mod(i,2).eq.1) c1d(i,j,k)=-c1d(i,j,k)
                  endif
c                  write(0,*) i,j,k,c1d(i,j,k)
              enddo
          enddo
      enddo


c
c construct the convolution matrix cmat, which gives the relation
c result(k3)= Sum cmat(k3,k1,k2) shape1(k1) shape2(k2)
c
      k1=0
      do nn1=0,n
          do j1=0,nn1
              i1=nn1-j1
              k1=k1+1
              k2=0
              do nn2=0,n
                  do j2=0,nn2
                      i2=nn2-j2
                      k2=k2+1
                      k3=0
                      do nn3=0,n
                          do j3=0,nn3
                              i3=nn3-j3
                              k3=k3+1
                              cmat(k3,k1,k2)=c1d(i3,i1,i2)*c1d(j3,j1,j2)
                          enddo
                      enddo
                  enddo
              enddo
          enddo
      enddo
      write(0,'('' Done.'')')
      return
      end
