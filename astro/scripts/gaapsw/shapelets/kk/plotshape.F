      subroutine plotshape(shape,nbinin,title)
      parameter (maxbin=200)
      real shape(msh),z(maxbin,maxbin),tr(6)
      character title*(*)

      nbin=abs(nbinin)
      if (nbin.gt.maxbin) then
          write(0,*) 'Cannot plot more than ',maxbin,' bins'
          nbin=maxbin
      endif
      call decode(shape,z,maxbin,sign(nbinin,nbin))
      zmax=-1e30
      zmin=-zmax
      do i=1,nbin
          do j=1,nbin
              zmax=max(zmax,z(i,j))
              zmin=min(zmin,z(i,j))
          enddo
      enddo
      tr(1)=-0.5*nbin-0.5
      tr(2)=1
      tr(3)=0
      tr(4)=tr(1)
      tr(5)=0
      tr(6)=1
      w=0.5*nbin
      call pgenv(-w,w,-w,w,1,1)
      call pgqch(ch)
      call pgsch(1.5*ch)
      call pglab(' ',' ',title)
      call pgsch(ch)
      do i=1,10
          c=0.1*i*zmax
          call pgcons(z,maxbin,maxbin,1,nbin,1,nbin, c,1,tr)
      enddo
      call pgsci(3)
      do i=1,10
          c=-0.1*i*zmax
          call pgcons(z,maxbin,maxbin,1,nbin,1,nbin, c,-1,tr)
      enddo
      call pgsci(1)
      end


      subroutine displayshape(shape,nbinin,title)
      parameter (maxbin=200)
      real shape(msh),z(maxbin,maxbin),tr(6),r(3),g(3),b(3),clev(3)
      character title*(*)
      data r,g,b,clev /1.,1.,0., 0.,1.,0., 0.,1.,1., 0.,0.090909,1./

      nbin=abs(nbinin)
      if (nbin.gt.maxbin) then
          write(0,*) 'Cannot plot more than ',maxbin,' bins'
          nbin=maxbin
      endif
      call decode(shape,z,maxbin,sign(nbinin,nbin))
      zmax=-1e30
      zmin=-zmax
      do i=1,nbin
          do j=1,nbin
              zmax=max(zmax,z(i,j))
              zmin=min(zmin,z(i,j))
          enddo
      enddo
      zmax=max(abs(zmax),abs(zmin))
      zmin=-0.1*zmax
      tr(1)=-0.5*nbin-0.5
      tr(2)=1
      tr(3)=0
      tr(4)=tr(1)
      tr(5)=0
      tr(6)=1
      w=0.5*nbin
      call pgenv(-w,w,-w,w,1,1)
      call pgqch(ch)
      call pgsch(1.5*ch)
      call pglab(' ',' ',title)
      call pgsch(ch)

      call pgscir(16,255)
      call pgctab(clev,r,g,b,3,1.,0.5)
      call pgimag(z,maxbin,maxbin,1,nbin,1,nbin,zmin,zmax,tr)

      end



      subroutine compare(z,mz,nx,ny,xc,yc,bg,r,sh)
      parameter (maxbin=200)
      real sh(msh),z(mz,mz),tr(6),zz(maxbin,maxbin)

      ix1=max(1,int(xc-r))
      ix2=min(nx,int(xc+r)+1)
      iy1=max(1,int(yc-r))
      iy2=min(ny,int(yc+r)+1)

      zmax=-1e30
      zmin=-zmax
      do ix=ix1,ix2
          do iy=iy1,iy2
              zmin=min(z(ix,iy),zmin)
              zmax=max(z(ix,iy),zmax)
          enddo
      enddo
      tr(1)=-1
      tr(2)=1
      tr(3)=0
      tr(4)=-1
      tr(5)=0
      tr(6)=1
      call pgenv(xc-r,xc+r,yc-r,yc+r,1,0)
      call pggray(z,mz,mz,ix1,ix2,iy1,iy2,zmin,zmax,tr)
      call pglab(' ',' ','IMAGE')

      nplot=min(199,1+2*nint(r))
      call decode(sh,zz,maxbin,nplot)
      call pgenv(1.,real(nplot),1.,real(nplot),1,0)
      call pggray(zz,maxbin,maxbin,1,nplot,1,nplot,zmin-bg,zmax-bg,tr)
      call pglab(' ',' ','FIT')

      idx=nint(xc-nplot/2)-1
      idy=nint(yc-nplot/2)-1
      do ix=ix1,ix2
          do iy=iy1,iy2
              zz(ix-idx,iy-idy)=z(ix,iy)-zz(ix-idx,iy-idy)
          enddo
      enddo
      call pgenv(1.,real(nplot),1.,real(nplot),1,0)
      range=0.1*(zmax-zmin)
      call pggray(zz,maxbin,maxbin,1,nplot,1,nplot,
     :    bg-range,bg+range,tr)
      call pglab(' ',' ','DIFF')
      
      return
      end


      subroutine showpatch(z,mx,my,nx,ny,x,y,rplot,bg,fg)
      real z(mx,my),tr(6)

      tr(1)=0
      tr(2)=1
      tr(3)=0
      tr(4)=0
      tr(5)=0
      tr(6)=1

      call pgenv(x-0.5*rplot-2,x+0.5*rplot+2,y-0.5*rplot-2,y+0.5*rplot+2,1,0)
      nx1=max(1,int(x-0.5*rplot))
      nx2=min(nx,1+int(x+0.5*rplot))
      ny1=max(1,int(y-0.5*rplot))
      ny2=min(ny,1+int(y+0.5*rplot))
      call pggray(z,mx,my,nx1,nx2,ny1,ny2,fg,bg,tr)
      return
      end

