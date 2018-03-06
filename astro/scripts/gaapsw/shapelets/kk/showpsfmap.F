      real sh(msh),psfcoef(mmpoly2,msh-2)

      read(*,*) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(n+1)*(n+2)/2
         read(*,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo
      nplot=max(10,nint(6*beta))
      np=max(4,2*nfit)
      call pgbeg(0,'psfmap.ps/cps',np,np)
      call pgslw(2)
      call pgsch(1.5)
      do iy=np,1,-1
         yfield=ymax/np*(iy-0.5)
         do ix=1,np
            xfield=xmax/np*(ix-0.5)
            call interpolpsf(xfield,yfield,xmax,ymax,psfcoef,n,beta,nfit,sh)
            call plotshape(sh,nplot,'PSF')
         enddo
      enddo
      call pgend
      end

