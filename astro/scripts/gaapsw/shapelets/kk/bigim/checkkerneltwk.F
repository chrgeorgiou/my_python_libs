c     test whether kernel integrates to 1 everywhere.
c     this version for double-shapelet map (from psfcat2gauskertwk and fitpsfmaptwk)
c     [can be extended with calculation of centroid as well.]

c read in kernel file from stdin

      real sh(msh),psfcoef(mmpoly2,2*msh-4)

 1    read(*,*,err=1) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      nin=(n+1)*(n+2)/2
      do k=1,nin
         read(*,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo
      write(0,*) 'Inner shapelet has n,beta=',n,beta
 2    read(*,*,err=2) nout,betaout
      read(*,*) nfitout,xmaxout,ymaxout
      read(*,*)
      do k=1,(nout+1)*(nout+2)/2
         read(*,*) (psfcoef(kfit,k+nin),kfit=1,(nfitout+1)*(nfitout+2)/2)
         read(*,*)
      enddo
      write(0,*) 'Outer shapelet has n,beta=',nout,betaout

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
            sum=seriessum(sh)
            call interpolpsf(xfield,yfield,xmaxout,ymaxout,psfcoef(1,1+nin),nout,betaout,nfitout,sh)
            sum=seriessum(sh)+sum
            write(*,*) ix,iy,sum
         enddo
      enddo
      end

