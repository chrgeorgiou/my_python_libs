c     read list of star positions, measure deviation of PSF from the
c     target Gaussian whose dispersion is given in GPSFSIG keyword or in
c     separate file gpsfsig.dat

      parameter (maxpix=20000,pi=3.1415926535)
      logical fitshavekey
      character*5 siglab
      real pix(maxpix)

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
      call pgbeg(0,'checkgpsf.ps/cps',1,1)
      call pgenv(0.,8*gpsfsig,-0.1*pk,1.3*pk,0,0)
      call pglab('Radius [pix]','Normalized Star','PSF vs. Gaussian')
      call pgmove(0.,0.)
      do ir=0,80
          call pgdraw(ir*0.1*gpsfsig,pk*exp(-0.5*ir*ir*0.01))
      enddo
      write(siglab,'(f5.2)') gpsfsig
      call pgptxt(7*gpsfsig,1.1*pk,0.,1.,'PSF SIGMA '//siglab)
      do i=1,mmcat
    3     read(*,*,err=3,end=4) x,y
          ix=nint(x)
          iy=nint(y)
          call igetsecpad(51,nx,ny,ix-np,iy-np,2*np+1,2*np+1,pix,maxpix,0.)
          call fitbg(pix,2*np+1,2*np+1,2*np+1,x-ix+1,y-iy+1,8*gpsfsig,bg)
          do ip=1,(2*np+1)**2
             pix(ip)=pix(ip)-bg
          enddo
          call plotpsf(pix,2*np+1,2*np+1,x-ix+np+1,y-iy+np+1,gpsfsig)
      enddo
    4 continue
      call pgmove(0.,pk)
      call pgsci(3)
      call pgslw(5)
      do ir=0,80
          call pgdraw(ir*0.1*gpsfsig,pk*exp(-0.5*ir*ir*0.01))
      enddo
      call pgend
      write(0,*) 'Plotted ',i-1,' objects.'
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
              if (dr2.lt.(8*sig)**2) then
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
              if (dr2.lt.(8*sig)**2) then
                  ee=exp(-dr2/(2*sig**2))/(2*3.1415926535*sig**2)
                  call pgpt1(sqrt(dr2),pix(ix,iy)/fac,-1)
c                  call pgpt1(sqrt(dr2),100*(pix(ix,iy)/fac-ee),-1)
              endif
          enddo
      enddo
      return
      end
