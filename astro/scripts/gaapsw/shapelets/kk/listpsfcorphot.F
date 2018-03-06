c     list-driven GaaP photometry 
c
c     open fits file photimage.fits, and error file photerr.fits
c     then read sets of x,y,beta,n,(q(i),i=1..n) and perform shapelet
C     expansion on position x,y with shapelet scale beta
c     produce list of psfcorphot values for this shapelet at scales q(i)

c     PSF information is gathered from photpsf.map file

      program listpsfcorphot

      parameter (maxpsf=1000)
      dimension shobs(msh),shpsf(msh)
      dimension psfcoef(mmpoly2,msh-2),psfrescoef(mmpoly2)
      dimension pix(mmpix,mmpix),pixerr(mmpix,mmpix)
      
c first get the PSF map coefficients

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

      n=npsf

c open the fits file

      call get2err('inimage.fits','err.fits',nx,ny,pix,pixerr,mmpix,mmpix)
c     read through the list of sources, and process each one:
c     fit beta
c     make a shapelet expansion
c     calculate the PSF
c     compute the GaaP fluxes

      betaprev=-1
      do i=1,mmcat
    1     read(*,*,err=1,end=99) x,y,f,fe,fw,a,b,th,id,iflg,beta,qq
          nread=nread+1
          if (beta.ne.betaprev) then
              call getpsfresidmap(beta,betapsf,qq,psfrescoef,nfit)
              betaprev=beta
              write(0,*) 'Got PSF resid corr factor for beta=',beta
          endif
          fitradius=max(10.,4*beta)
          if (x.gt.fitradius+1 .and. x.lt.nx-fitradius-1
     :        .and. y.gt.fitradius+1 .and. y.lt.ny-fitradius-1) then
              bg=0.
              call fitshapeletserr(pix,pixerr,mmpix,nx,ny,x,y,fitradius
     :            ,bg,shobs,sherr,n,beta)
              call interpolpsf(x,y,xmax,ymax,psfcoef,npsf,betapsf,nfit,shpsf)
              fq=psfcorphot3(shobs,sherr,shpsf,qq,fqerr)
c squared dispersion of gaussian wt fn for residuals
              filtdisp2=max(2*qq*qq-(betapsf/1.3)**2,0.5)
              fqresid=(qq*qq/filtdisp2)*
     :            pixresid(shobs,x,y,pix,mmpix,nx,ny,filtdisp2)
              write(*,1001) x,y,f,fe,fw,a,b,th,id,iflg,
     :            beta,qq,fq+fqresid,fqerr
 1001         format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i6,i5,
     :            2f10.3,3g13.4)
          else
              goto 1
          endif
          
   99 enddo

      end



      subroutine getpsfresidmap(beta,betapsf,qq,coefs,nfit)
      dimension coefs(mmpoly2)

