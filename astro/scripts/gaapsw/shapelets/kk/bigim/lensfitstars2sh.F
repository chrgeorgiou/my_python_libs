      parameter (mpix=32,mstars=10000)
      dimension stars(mpix,mpix,mstars),xs(mstars),ys(mstars)
      dimension pixerr(mpix,mpix),sh(msh),sh2(msh)
      dimension betastar(mstars)
      dimension xccd(36),yccd(36)
      character cubename*132

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

c first fit a common Gaussian radius
c then do shapelet fits
c write out result but with the original x,y coordinates

      fitradius=min(nx,ny)/2-1
      call readorders(nsh,nfit)

      do istar=1,nstars
         call fitbeta(stars(1,1,istar),32,32,32,17.,17.,fitradius,0.,betastar(istar))
      enddo
      beta=fmedian(nstars,betastar)*1.3

      write(0,*) 'Adopted beta for stars [pix]: ',beta

      do istar=1,nstars
         do j=1,ny
            do i=1,nx
               if (stars(i,j,istar).gt.-0.099) then
                  pixerr(i,j)=1
               else
                  pixerr(i,j)=0
               endif
            enddo
         enddo
         call fitshapeletserr(stars(1,1,istar),pixerr,mpix,nx,ny,17.,17.,fitradius,1.e-30,sh,sherr,nsh,beta)
         call centershape(sh,sh2,dx,dy)
         sum=seriessumerr(sh2,sherr,sumerr)
         write(*,1001) xs(istar),ys(istar),1.,1e-3,2.3*betastar(istar),1.,1.,0.,istar,0
 1001    format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i10,i5)
         write(*,*) sum,dx,dy,0.,fitradius,sherr/sum
         write(*,*) beta,nsh
         write(*,1002) (sh2(k)/sum,k=3,(nsh+1)*(nsh+2)/2+2)
 1002    format(msh(g15.6))
         write(*,*)
c calculate residuals of non-masked pixels
         do j=1,mpix
            do i=1,mpix
               if (pixerr(i,j).le.0.) then
                  stars(i,j,istar)=0
               else
                  stars(i,j,istar)=stars(i,j,istar)-shapeletpix(i-17.+dx,j-17.+dy,sh2)
               endif
            enddo
         enddo
      enddo

      call put3('lensfitstar_resid.fits',nx,ny,nstars,stars,mpix,mpix,mstars)

      write(0,*) 'lensfitstars2sh read ',nstars,' PSFs, fitted shapelets to order ',nsh
      write(0,*) 'Residuals image cube in lensfitstar_resid.fits .'
      end
