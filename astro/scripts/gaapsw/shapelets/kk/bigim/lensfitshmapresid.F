      parameter (mpix=32,mstars=10000)
      dimension stars(mpix,mpix,mstars),xs(mstars),ys(mstars)
      dimension pixerr(mpix,mpix),sh(msh)
      dimension patch(mpix,mpix)
      dimension psfcoef(mmpoly2,msh-2)
      dimension xccd(36),yccd(36)
      character cubename*132

c read table with CCD positions
      open(10,file='ccdoffsets.dat')
      read(10,*) (xccd(i),yccd(i),i=1,36)
      close(10)

c read in file names from standard input. Each file contains lensfit stars and position table.
c each filename is followed by a line identifying the CCD (0 if no offset)
      ntot=0
 1    read(*,'(a)',err=1,end=2) cubename
      call get3(cubename,nx,ny,nstars,stars(1,1,ntot+1),mpix,mpix,mstars)
      read(*,*) iccd
      call tabget(cubename,'star info',ns2,'X',xs(ntot+1),mstars)
      call tabget(cubename,'star info',ns2,'Y',ys(ntot+1),mstars)
      if (ns2.ne.nstars) then
          write(0,*) 'Table and image have difft. nos. of stars'
          stop
      endif
      if (iccd.gt.0) then
         do i=1,ns2
            xs(ntot+i)=xs(ntot+i)+xccd(iccd)
            ys(ntot+i)=ys(ntot+i)+yccd(iccd)
         enddo
      endif
      ntot=ntot+nstars
      goto 1
 2    continue
      nstars=ntot

c Get psf map from psf.map file (or link that points to it)

      open(10,file='psf.map',status='old')
      read(10,*) npsf,betapsf
      read(10,*) nfit,xmax,ymax
      read(10,*)
      do k=1,(npsf+1)*(npsf+2)/2
         read(10,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(10,*)
      enddo
      close(10)
      write(0,*) 'Read psf model from psf.map'
      write(0,*) 'Order for shape, spatial variation; beta:',npsf,nfit,betapsf
      
      fitradius=min(nx,ny)/2-1
      do istar=1,nstars
c     repeat shapelet fit to determine x,y offset
          do j=1,ny
              do i=1,nx
                  if (stars(i,j,istar).gt.-0.099) then
                      pixerr(i,j)=1
                  else
                      pixerr(i,j)=0
                  endif
              enddo
          enddo
          call fitshapeletserr(stars(1,1,istar),pixerr,mpix,nx,ny,17.,17.,fitradius,1.e-30,sh,sherr,npsf,betapsf)
          call centershape(sh,sh,dx,dy)
          call interpolpsf(xs(istar),ys(istar),xmax,ymax,psfcoef,npsf,betapsf,nfit,sh)
          do j=1,mpix
              do i=1,mpix
                  patch(i,j)=shapeletpix(i-17.+dx,j-17.+dy,sh)
              enddo
          enddo
          top=0
          bot=0
          do j=1,mpix
              do i=1,mpix
                  if (stars(i,j,istar).gt.-0.099) then
                      top=top+patch(i,j)*stars(i,j,istar)
                      bot=bot+patch(i,j)**2
                  endif
              enddo
          enddo
          fluxfact=top/bot
          do j=1,mpix
              do i=1,mpix
                  if (stars(i,j,istar).gt.-0.099) then
                      stars(i,j,istar)=stars(i,j,istar)-patch(i,j)*fluxfact
                  else
                      stars(i,j,istar)=0
                  endif
              enddo
          enddo
      enddo

      call put3('lensfitmapresid.fits',nx,ny,nstars,stars,mpix,mpix,mstars)
      write(0,*) 'Wrote out residual image cube lensfitmapresid.fits'
      end

