c take PSF shape, turn it into a PSF matrix ready to operate on other shapes
c with scale beta. Result's scale is betaout

      subroutine setpsfmat(shpsf,betaout,beta,psfmat)

      real shpsf(msh),psfmat(msh-2,msh-2),cmat(msh-2,msh-2,msh-2)
      double precision pp
      common /convcoef/ cmat,b1,b2,bout,norder

      betapsf=shpsf(1)
      n=shpsf(2)
      if (b1.ne.betapsf .or. b2.ne.beta .or. bout.ne.betaout 
     :    .or. norder.ne.n) then
          call setcmatrix(cmat,betapsf,beta,betaout,n)
c update info in common block on which parameters are in there     
          b1=betapsf
          b2=beta
          bout=betaout
          norder=n
      endif

c make the PSF matrix by contracting with cmat
      nsh=(n+1)*(n+2)/2
      do k3=1,nsh
          do k2=1,nsh
              pp=0
              do k1=1,nsh
                  pp=pp+cmat(k3,k1,k2)*shpsf(k1+2)
              enddo
              psfmat(k3,k2)=pp
          enddo
      enddo
c      write(0,*) 'Made PSF matrix for n, betas ',n,betaout,betapsf,beta
c      write(0,*)

c      do i=1,nsh
c         write(0,'(15g13.5)') (psfmat(i,j),j=1,nsh)
c      enddo
c      write(0,*)
      return
      end
