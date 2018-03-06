      parameter (mmax=40)
      real ztot(mmax,mmax),zg1(mmax,mmax),zg2(mmax,mmax),
     :    zg1s(mmax,mmax),zg2s(mmax,mmax),zgg(mmax,mmax)
      real tr(6),tr2(6)

      tr(1)=0.
      tr(2)=0.25
      tr(3)=0
      tr(4)=1
      tr(5)=0.
      tr(6)=0.1
      tr2(2)=0.3/mmax
      tr2(1)=-0.15-tr2(2)/2
      tr2(3)=0
      tr2(4)=tr2(1)
      tr2(5)=0
      tr2(6)=tr2(2)

      fwmin=1e30
      fmin=1e30
      call zero(ztot,mmax*mmax)
      call zero(zg1,mmax*mmax)
      call zero(zg2,mmax*mmax)
      call zero(zg1s,mmax*mmax)
      call zero(zg2s,mmax*mmax)
      call zero(zgg,mmax*mmax)

    1 read(*,*,err=1,end=2) x,y,f,fe,fw,a,b,pa,g1,g2,sg1,sg2,covgg,id,iflg
      if (iflg.ne.0) goto 1
      ix=nint((fw-tr(1))/tr(2))
      iy=nint((log10(f)-tr(4))/tr(6))
      if (ix.ge.1 .and. ix.le.mmax .and. iy.ge.1 .and. iy.le.mmax) then
          fwmin=min(fwmin,fw)
          fmin=min(fmin,log10(f))
          ztot(ix,iy) =ztot(ix,iy) +1
          zg1(ix,iy)  =zg1(ix,iy)  +g1
          zg1s(ix,iy) =zg1s(ix,iy) +g1*g1
          zg2(ix,iy)  =zg2(ix,iy)  +g2
          zg2s(ix,iy) =zg2s(ix,iy) +g2*g2
          wt=1./(2*0.05**2+sg1**2+sg2**2)
          ig1=nint((g1-tr2(1))/tr2(2))
          ig2=nint((g2-tr2(4))/tr2(6))
          if (ig1.le.mmax.and.ig2.le.mmax.and.ig1.ge.1.and.ig2.ge.1)  zgg(ig1,ig2)=zgg(ig1,ig2)+wt
      endif
      
      goto 1

    2 continue

      do j=1,mmax
          do i=1,mmax
              zz=ztot(i,j)
              if (zz.ge.10) then 
                  zg1(i,j)=zg1(i,j)/zz
                  zg2(i,j)=zg2(i,j)/zz
                  zg1s(i,j)=zz/(0.05**2+zg1s(i,j)/zz-zg1(i,j)**2)
                  zg2s(i,j)=zz/(0.05**2+zg2s(i,j)/zz-zg2(i,j)**2)
              else
                  zg1(i,j)=-999
                  zg2(i,j)=-999
                  zg1s(i,j)=0
                  zg2s(i,j)=0
              endif
          enddo
      enddo

      call pgbeg(0,'/xs',3,2)

      call pgscir(16,99)
      call pgscr(16,0.,0.,0.)
      do i=17,99
         t=(i-16)/83.*3.1415
c         cr=sin(t/2)
c         cg=sin(t)
c         cb=1-cr
         cr=sin(0.5*t)**2
         cg=sin(1.*t)**2
         cb=sin(1.5*t)**2
         call pgscr(i,cr,cg,cb)
      enddo

      call pgenv(fwmin,tr(1)+mmax*(tr(2)),fmin,tr(4)+mmax*tr(6),0,0)
      call getrange(ztot,mmax,mmax,zmin,zmax)
      call pgimag(ztot,mmax,mmax,1,mmax,1,mmax,0.,zmax,tr)
      call pgwedg('RI',1.5,3.,0.,zmax,' ')
      call pglab ('FWHM','Lg(Flx)','NUMBER OF GALAXIES')

      call pgenv(fwmin,mmax*(tr(2)),fmin,tr(4)+mmax*tr(6),0,0)
      call getrange(zg1,mmax,mmax,zmin,zmax)
      zmin=-0.05
      zmax=0.05
      call pgimag(zg1,mmax,mmax,1,mmax,1,mmax,zmin,zmax,tr)
      call pgwedg('RI',1.5,3.,zmin,zmax,' ')
      call pglab ('FWHM','Lg(Flx)','MEAN G1')

      call pgenv(fwmin,mmax*(tr(2)),fmin,tr(4)+mmax*tr(6),0,0)
      call getrange(zg2,mmax,mmax,zmin,zmax)
      zmin=-0.05
      zmax=0.05
      call pgimag(zg2,mmax,mmax,1,mmax,1,mmax,zmin,zmax,tr)
      call pgwedg('RI',1.5,3.,zmin,zmax,' ')
      call pglab ('FWHM','Lg(Flx)','MEAN G2')

      call pgenv(-0.17,.17,-.17,.17,1,0)
      call getrange(zgg,mmax,mmax,zmin,zmax)
      write(*,*) zmin,zmax
      zmin=0
      call pgimag(zgg,mmax,mmax,1,mmax,1,mmax,zmin,zmax,tr2)
      call pgwedg('RI',1.5,3.,zmin,zmax,' ')
      call pglab ('G1','G2','WEIGHT')

      call pgenv(fwmin,mmax*(tr(2)),fmin,tr(4)+mmax*tr(6),0,0)
      call getrange(zg1s,mmax,mmax,zmin,zmax)
      zmin=0
      call pgimag(zg1s,mmax,mmax,1,mmax,1,mmax,zmin,zmax,tr)
      call pgwedg('RI',1.5,3.,zmin,zmax,' ')
      call pglab ('FWHM','Lg(Flx)','WEIGHT G1')

      call pgenv(fwmin,mmax*(tr(2)),fmin,tr(4)+mmax*tr(6),0,0)
      call getrange(zg2s,mmax,mmax,zmin,zmax)
      zmin=0
      call pgimag(zg2s,mmax,mmax,1,mmax,1,mmax,zmin,zmax,tr)
      call pgwedg('RI',1.5,3.,zmin,zmax,' ')
      call pglab ('FWHM','Lg(Flx)','WEIGHT G2')

      call pgend

      end


      subroutine getrange(z,n,m,zmin,zmax)
      real z(n,m)
      zmax=z(1,1)
      zmin=zmax
      do j=1,m
          do i=1,n
              zmax=max(zmax,z(i,j))
              zmin=min(zmin,z(i,j))
          enddo
      enddo
      return
      end


      subroutine zero(a,n)
      real a(n)
      do i=1,n
          a(i)=0
      enddo
      return
      end
