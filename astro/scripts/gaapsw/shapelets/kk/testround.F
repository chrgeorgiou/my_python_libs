      real s1(msh),s2(msh),z(30,30),tr(6)

      beta=3
      n=10
      s1(1)=beta
      s1(2)=n
      do k=3,msh
         s1(k)=0
      enddo
      do nn=0,n,2
         write(0,*) 'Coef c',nn
         read(*,*) cn
         k0=2+(nn)*(nn+1)/2 + 1
         cc=1
c         do i=1,nn
c             cc=cc*i
c         enddo
c         cc=sqrt(cc)
         s1(k0)=cc*cn
         do kk=2,nn,2
            cc=cc*(nn/2-kk/2+1)/(kk/2)*sqrt((kk-1)*kk /real((nn-kk+2)*(nn-kk+1)))
            s1(k0+kk)=cc*cn 
         enddo
      enddo
      ktot=2+(n+1)*(n+2)/2
      call writeshapen('Shape',s1)

      call decode(s1,z,30,25)

      zmax=-1e30
      zmin=1e30
      do i=1,25
          zmin=min(zmin,z(i,13))
          zmax=max(zmax,z(i,13))
      enddo

      call pgbeg(0,'/xs',2,2)
      call pgenv(-15.,15.,-15.,15.,1,1)
      tr(1)=-13.
      tr(2)=1.
      tr(3)=0
      tr(4)=-13
      tr(5)=0
      tr(6)=1
      do i=1,20
          c=0.05*i*zmax
          call pgcons(z,30,30,1,25,1,25,c,1,tr)
      enddo

      call pgenv(0.,13.,zmin,zmax,0,1)
      call pgmove(tr(1)+tr(2),z(1,13))
      do i=2,25
          call pgdraw(tr(1)+i*tr(2),z(i,13))
      enddo

      call shear1(s1,s2)
      do k=3,2+(n+1)*(n+2)/2
          s2(k)=s1(k)+0.05*s2(k)
      enddo
      call decode(s2,z,30,25)
      call pgenv(-15.,15.,-15.,15.,1,1)
      do i=1,20
          c=0.05*i*zmax
          call pgcons(z,30,30,1,25,1,25,c,1,tr)
      enddo
      

      call shear2(s1,s2)
      do k=3,2+(n+1)*(n+2)/2
          s2(k)=s1(k)+0.05*s2(k)
      enddo
      call decode(s2,z,30,25)
      call pgenv(-15.,15.,-15.,15.,1,1)
      do i=1,20
          c=0.05*i*zmax
          call pgcons(z,30,30,1,25,1,25,c,1,tr)
      enddo

      call pgend
      end
