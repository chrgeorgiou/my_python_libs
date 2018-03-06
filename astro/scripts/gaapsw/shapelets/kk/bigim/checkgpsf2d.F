c     read list of star positions, measure deviation of PSF from the
c     target Gaussian whose dispersion is given in GPSFSIG keyword or in
c     separate file gpsfsig.dat
c produce plot which shows these profile as a 2d grid.

      parameter (maxpix=20000,pi=3.1415926535)
      logical fitshavekey
      character*5 siglab
      real pix(maxpix),x(mmcat),y(mmcat)

      call openfits('inimage.fits',51,nx,ny)
      if (fitshavekey(51,'GPSFSIG')) then
          gpsfsig=fitskeyr(51,'GPSFSIG')
          goto 2
      else
          open(10,file='gpsfsig.dat',status='old',err=1)
          read(10,*) gpsfsig
          write(0,*) '***Did not find GPSFSIG value in fits header.'
          write(0,*) 'Took value from file gpsfsig.dat'
          goto 2
      endif
      
    1 write(0,*) '***Did not find GPSFSIG value in fits header.'
      write(0,*) '***You can specify one in text file gpsfsig.dat'
      stop
    2 continue      
      write(0,*) 'Input image wants Gaussian PSF of width ',gpsfsig

      np=nint(8*gpsfsig)
      pk=1./(2.*pi*gpsfsig**2)

      npanel=5
      call pgbeg(0,'checkgpsf2d.ps/vcps',npanel,npanel)

      do i=1,mmcat
 3       read(*,*,err=3,end=4) x(i),y(i)
      enddo
 4    continue
      ncat=i-1
      
      do iypanel=npanel,1,-1
          do ixpanel=1,npanel
            call pgenv(0.,5*gpsfsig,-0.1*pk,1.3*pk,0,0)
            call pglab('Radius [pix]','Normalized Star',' ')
            call pgmove(0.,0.)
            do ir=0,50
               call pgdraw(ir*0.1*gpsfsig,pk*exp(-0.5*ir*ir*0.01))
            enddo
            if (ixpanel.eq.npanel .and. iypanel.eq.npanel) then
                write(siglab,'(f5.2)') gpsfsig
                call pgsch(npanel*1.)
                call pgptxt(4.5*gpsfsig,1.1*pk,0.,1.,'SIG='//siglab)
                call pgsch(1.)
            endif
            do i=1,ncat
               if (int(npanel*x(i)/nx).eq.ixpanel-1 .and. int(npanel*y(i)/ny).eq.iypanel-1) then
                  ix=nint(x(i))
                  iy=nint(y(i))
                  call igetsecpad(51,nx,ny,ix-np,iy-np,2*np+1,2*np+1,pix,maxpix,0.)
                  call fitbg(pix,2*np+1,2*np+1,2*np+1,x(i)-ix+1,y(i)-iy+1,8*gpsfsig,bg)
                  do ip=1,(2*np+1)**2
                     pix(ip)=pix(ip)-bg
                  enddo
                  call plotpsf(pix,2*np+1,2*np+1,x(i)-ix+np+1,y(i)-iy+np+1,gpsfsig)
               endif
            enddo
            call pgmove(0.,pk)
            call pgsci(3)
            call pgslw(1)
            do ir=0,50
c               call pgdraw(ir*0.1*gpsfsig,pk*exp(-0.5*ir*ir*0.01))
               call pgpt1(ir*0.1*gpsfsig,pk*exp(-0.5*ir*ir*0.01),24)
            enddo
            call pgsci(1)
            call pgslw(1)
         enddo
      enddo
      call pgend

c repeat plot but now one plot per page
      call pgbeg(0,'checkgpsf2d_pages.ps/cps',1,1)
      call pgask(.false.)
      do iypanel=npanel,1,-1
         do ixpanel=1,npanel
            call pgslw(4)
            call pgsch(2.)
            call pgenv(0.,5*gpsfsig,-0.1*pk,1.3*pk,0,0)
            call pglab('Radius [pix]','Normalized Star',' ')
            call pgsch(1.)
            call pgmove(0.,0.)
            do ir=0,50
               call pgdraw(ir*0.1*gpsfsig,pk*exp(-0.5*ir*ir*0.01))
            enddo
c            if (ixpanel.eq.npanel .and. iypanel.eq.npanel) then
c                write(siglab,'(f5.2)') gpsfsig
c                call pgsch(npanel*1.)
c                call pgptxt(4.5*gpsfsig,1.1*pk,0.,1.,'SIG='//siglab)
c                call pgsch(1.)
c            endif
            call pgslw(1)
            do i=1,ncat
               if (int(npanel*x(i)/nx).eq.ixpanel-1 .and. int(npanel*y(i)/ny).eq.iypanel-1) then
                  ix=nint(x(i))
                  iy=nint(y(i))
                  call igetsecpad(51,nx,ny,ix-np,iy-np,2*np+1,2*np+1,pix,maxpix,0.)
                  call fitbg(pix,2*np+1,2*np+1,2*np+1,x(i)-ix+1,y(i)-iy+1,8*gpsfsig,bg)
                  do ip=1,(2*np+1)**2
                     pix(ip)=pix(ip)-bg
                  enddo
                  call plotpsf(pix,2*np+1,2*np+1,x(i)-ix+np+1,y(i)-iy+np+1,gpsfsig)
               endif
            enddo
            call pgmove(0.,pk)
            call pgsci(3)
            call pgslw(1)
            do ir=0,50
c               call pgdraw(ir*0.1*gpsfsig,pk*exp(-0.5*ir*ir*0.01))
               call pgpt1(ir*0.1*gpsfsig,pk*exp(-0.5*ir*ir*0.01),24)
            enddo
            call pgsci(1)
            call pgslw(1)
         enddo
      enddo
      call pgend

      write(0,*) 'Plotted ',ncat,' objects.'
      end

      subroutine plotpsf(pix,nx,ny,xc,yc,sig)
      real pix(nx,ny)
      double precision sum,sum2
      sum=0
      sum2=0
      do ix=1,nx
          dx=ix-xc
          do iy=1,ny
              dy=iy-yc
              dr2=dx*dx+dy*dy
              if (dr2.lt.(5*sig)**2) then
                  ee=exp(-dr2/(2*sig**2))
                  sum=sum+ee*pix(ix,iy)
                  sum2=sum2+ee*ee
              endif
          enddo
      enddo
      fac=sum/sum2 * (2*3.1415926535) * sig**2
      do ix=1,nx
          dx=ix-xc
          do iy=1,ny
              dy=iy-yc
              dr2=dx*dx+dy*dy
              if (dr2.lt.(5*sig)**2) then
                  ee=exp(-dr2/(2*sig**2))/(2*3.1415926535*sig**2)
                  call pgpt1(sqrt(dr2),pix(ix,iy)/fac,7)
c                  call pgpt1(sqrt(dr2),100*(pix(ix,iy)/fac-ee),-1)
              endif
          enddo
      enddo
      return
      end
