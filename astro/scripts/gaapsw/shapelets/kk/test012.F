c test of decomposition into polars

      real z1(41,41),z2(41,41),sh1(msh),sh2(msh),sh3(msh),sh4(msh),sh5(msh)
      parameter (mm012=5*(mm/2)-1)
      real sh1_012(mm012)

      call pgbeg(0,'/xs',3,2)
      sig=2.5
      sig2=sig*sig
      n=10
      do i=1,41
          x=i-21
          do j=1,41
              y=j-21
              z1(i,j)=exp(-(x**2/0.8**2+y**2/1.2**2)/2/sig2)/(2*3.1415926535*sig2)
          enddo
      enddo
      call encode(z1,41,41,41,21.,21.,20.,0.,sh1,n,sig)
      call writeshapen('shape1',sh1)
      call plotshape(sh1,40,'1')
      call to012(sh1,sh1_012,10)

c now build some cpts, see if they look OK.
      norder=16
      do nn=0,norder,2
         do k=3,msh
            sh1(k)=0
         enddo
         sh1(1)=15.
         sh1(2)=norder
         m=nn/2
         k0=3+nn*(nn+1)/2
         cc=1
         sh1(k0)=cc
         write(0,*) 'rnd',nn,0,cc
         sum2=cc**2
         do i=1,m
            i2=2*i
            mi2=2*(m-i)
            cc=cc*(m-i+1)/i*sqrt(i2*(i2-1)/(mi2+2.)/(mi2+1.))
            sum2=sum2+cc**2
            sh1(k0+i2)=cc
            write(0,*) 'rnd',nn,i2,cc
         enddo
         call plotshape(sh1,150,' ')
         call to012(sh1,sh1_012,norder)
         write(0,*)
      enddo
c now the cos cpts
      do nn=1,norder-1,2
         do k=3,msh
            sh1(k)=0
         enddo
         sh1(1)=15.
         sh1(2)=norder
         m=(nn-1)/2
         k0=3+nn*(nn+1)/2
         cc=1
         sh1(k0)=cc
         write(0,*) 'cos',nn,0,cc
         sum2=cc**2
         do i=1,m
            i2=2*i
            mi2=2*(m-i)
            cc=cc*(m-i+1)/i*sqrt(i2*(i2-1)/(mi2+3.)/(mi2+2.))
            sum2=sum2+cc**2
            sh1(k0+i2)=cc
            write(0,*) 'cos',nn,i2,cc
         enddo
         call plotshape(sh1,150,' ')
         call to012(sh1,sh1_012,norder)
         write(0,*)
      enddo
c now the sin cpts
      do nn=1,norder-1,2
         do k=3,msh
            sh1(k)=0
         enddo
         sh1(1)=15.
         sh1(2)=norder
         m=(nn-1)/2
         k0=3+nn*(nn+1)/2
         cc=1
         sh1(k0+1)=cc
         write(0,*) 'sin',nn,1,cc
         sum2=cc**2
         do i=1,m
            i2=2*i
            mi2=2*(m-i)
            cc=cc*(m-i+1)/i*sqrt((i2+1)*i2/(mi2+2.)/(mi2+1.))
            sum2=sum2+cc**2
            sh1(k0+i2+1)=cc
            write(0,*) 'sin',nn,i2+1,cc
         enddo
         call plotshape(sh1,150,' ')
         call to012(sh1,sh1_012,norder)
         write(0,*)
      enddo
c now the cos2 cpts
      do nn=2,norder-2,2
         do k=3,msh
            sh1(k)=0
         enddo
         sh1(1)=15.
         sh1(2)=norder
         m=(nn-2)/2
         k0=3+nn*(nn+1)/2
         cc=1
         sh1(k0)=cc
         write(0,*) 'cos-2th',nn,0,cc
         sum2=cc**2
         cplus=cc
         fac=1
         do i=1,m+1
            i2=2*i
            mi2=2*(m-i)
            cmin=cplus
            cplus=cplus*(m-i+1)/i
            fac=fac*(i2-1)*i2/(mi2+4.)/(mi2+3.)
            cc=(cplus-cmin)*sqrt(fac)
            sum2=sum2+cc**2
            sh1(k0+i2)=cc
            write(0,*) 'cos-2th',nn,i2,cc
         enddo
         call plotshape(sh1,150,' ')
         call to012(sh1,sh1_012,norder)
         write(0,*)
      enddo
c now the sin2 cpts
      do nn=2,norder-2,2
         do k=3,msh
            sh1(k)=0
         enddo
         sh1(1)=15.
         sh1(2)=norder
         m=(nn-2)/2
         k0=3+nn*(nn+1)/2
         cc=1
         sh1(k0+1)=cc
         write(0,*) 'sin-2th',nn,1,cc
         sum2=cc**2
         do i=1,m
            i2=2*i
            mi2=2*(m-i)
            cc=cc*(m-i+1)/i*sqrt((i2+1)*i2/(mi2+3.)/(mi2+2.))
            sh1(k0+i2+1)=cc
            write(0,*) 'sin-2th',nn,i2+1,cc
         enddo
         call plotshape(sh1,150,' ')
         call to012(sh1,sh1_012,norder)
         write(0,*)
      enddo


      call pgend

      end
