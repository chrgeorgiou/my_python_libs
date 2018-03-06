      real sh(msh),sh2(msh),psfcoef(mmpoly2,2*msh-4)

 1    read(*,*,err=1) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      nin=(n+1)*(n+2)/2
      do k=1,nin
         read(*,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo
 2    read(*,*,err=2) nout,betaout
      read(*,*) nfitout,xmaxout,ymaxout
      read(*,*)
      do k=1,(nout+1)*(nout+2)/2
         read(*,*) (psfcoef(kfit,k+nin),kfit=1,(nfitout+1)*(nfitout+2)/2)
         read(*,*)
      enddo
      nplot=max(10,nint(10*beta))
      np=max(4,2*nfit)
      call pgbeg(0,'psfmap.ps/cps',np,np)
      call pgslw(2)
      call pgsch(1.5)
      do iy=np,1,-1
         yfield=ymax/np*(iy-0.5)
         do ix=1,np
            xfield=xmax/np*(ix-0.5)
            call interpolpsf(xfield,yfield,xmax,ymax,psfcoef,n,beta,nfit,sh)
            call interpolpsf(xfield,yfield,xmaxout,ymaxout,psfcoef(1,nin+1),nout,betaout,nfitout,sh2)
            call displayshapesum(sh,sh2,nplot,'PSF')
         enddo
      enddo
      call pgend
      end

