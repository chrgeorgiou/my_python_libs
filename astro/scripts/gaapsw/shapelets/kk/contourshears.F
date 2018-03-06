c read in a list of g1,g2,w; make contour plot of the distribution
      real z(200,200),tr(6)
      dg=0.0025
      tr(1)=-0.05-dg*0.5
      tr(2)=dg
      tr(3)=0
      tr(4)=-0.05-dg*0.5
      tr(5)=0
      tr(6)=dg
      nx=80
      ny=40
      zmax=0
      do i=1,999999999
          read(*,*,end=99) g1,g2,w
          ix=nint((g1-tr(1))/tr(2))
          iy=nint((g2-tr(4))/tr(6))
          if (ix.ge.1 .and. ix.le.nx .and.iy.ge.1 .and.iy.le.ny) then
              z(ix,iy)=z(ix,iy)+w
              zmax=max(zmax,z(ix,iy))
          endif
      enddo
   99 write(0,*) 'Read ',i-1,' sources. Zmax=',zmax
      call pgbeg(0,'/xw',1,1)
      call pgenv(tr(1)-0.5*dg,tr(1)+(nx+1.5)*dg,tr(4)-0.5*tr(6),tr(4)+(ny+1.5)*dg,1,1)
      do i=1,10
          call pgcont(z,200,200,1,nx,1,ny,zmax*i*0.1,1,tr)
      enddo
      call pgend
      end
