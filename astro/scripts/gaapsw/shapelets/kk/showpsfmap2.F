      real sh(msh),psfcoefrotated(msh-2,mmpoly2)

      read(*,*) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(n+1)*(n+2)/2
         read(*,*) (psfcoefrotated(k,kfit),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo
      nplot=max(10,nint(6*beta))
      np=max(4,2*nfit)
      call pgbeg(0,'psfmap2.ps/ps',nfit+1,nfit+1)
      call pgslw(2)
      call pgsch(1.5)
      kfit=1
      sh(1)=beta
      sh(2)=n
      do ix=0,nfit
         do iy=ix,0,-1
            do k=1,(n+1)*(n+2)/2
               sh(k+2)=psfcoefrotated(k,kfit)
            enddo
            call plotshape(sh,nplot,'PSF')
            kfit=kfit+1
         enddo
         do iy=ix+1,nfit
            call pgadvance
         enddo
      enddo
      call pgend
      end

