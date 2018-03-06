      real pix(mmpix,mmpix)
      real sh(msh)

      call pgbeg(0,'/xs',4,3)
      call pgsch(2.)
      call get2('inimage.fits',nx,ny,pix,mmpix,mmpix)
      do i=1,nx
          do j=1,ny
              pix(i,j)=pix(i,j)-7
          enddo
      enddo
      open(10,file='test.cat',status='old')
      do i=1,999999
    1     read(10,*,err=1,end=99) xo,yo,fo,foe,fwo,ao,bo,tho,id,iflg
          beta=fwo/2.3
          fitradius=max(10.,min(40.,8*beta))
          xo=nint(xo)
          yo=nint(yo)
          write(0,*)
          write(0,*) xo,yo,beta
c          call encode(pix,mmpix,nx,ny,xo,yo,ro,0.,sh,10,0.)
          bg=0.
          call fitshapelets(pix,mmpix,nx,ny,xo,yo,fitradius,bg,sh,-10,beta)
          call writeshapen('Sh:',sh)
c          write(0,'(''Sh: '',2f8.4/4x,3f8.4)') (sh(k)/sh(3),k=4,8)
          rplot=max(7.,sh(1)*4)
          call plotshape(sh,nint(2*rplot),'Shape')
          call compare(pix,mmpix,nx,ny,xo,yo,bg,1.5*rplot,sh)
      enddo
   99 continue
      call pgend
      end
