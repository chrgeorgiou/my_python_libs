c test whether kernel integrates to 1 everywhere.
c [can be extended with calculation of centroid as well.]

c read in kernel file from stdin

      real sh(msh),psfcoef(mmpoly2,msh-2)

      read(*,*) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(n+1)*(n+2)/2
         read(*,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo
      np=max(4,2*nfit)

c over a grid (density determined by spatial variation polynomial order
c nfit) reconstruct kernel and tabulate integral in format that can be used by showim:
c label for plot
c nx,ny
c x,y,value
c ...

      write(*,*)
      write(*,*) np,np
      do iy=1,np
         yfield=ymax/np*(iy-0.5)
         do ix=1,np
            xfield=xmax/np*(ix-0.5)
            call interpolpsf(xfield,yfield,xmax,ymax,psfcoef,n,beta,nfit,sh)
            write(*,*) ix,iy,seriessum(sh)
         enddo
      enddo
      end

