c     list-driven GaaP photometry 
c
c     open fits file inimage.fits, and error file inerr.fits
c     then read sets of sextractor sources,beta,q and perform shapelet
C     expansion on position x,y with shapelet scale beta
c     produce list of psfcorphot values for this shapelet at scales q

c     PSF information is gathered from psf.map file
c     PSF star list is taken from psf.sh file (for calculating PSF residuals)

      program listpsfcorphot

      parameter (maxpsf=1000)
      dimension shobs(msh),shpsf(msh)
      dimension psfcoef(mmpoly2,msh-2),psfrescoef(mmpoly2)
      dimension pix(mmpix,mmpix),pixerr(mmpix,mmpix)
      dimension xpsf(maxpsf),ypsf(maxpsf),sumpsf(maxpsf),bgpsf(maxpsf),
     :    res(maxpsf),wt(maxpsf),bestfit(maxpsf)
      
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

c get the list of PSF star positions from psf.sh
      open(10,file='psf.sh',status='old')
      do i=1,maxpsf
          read(10,*,end=11) xpsf(i),ypsf(i),(ttt,kkk=1,8)
          read(10,*) sumpsf(i),dxpsf,dypsf,bgpsf(i),(ttt,kkk=1,2)
          xpsf(i)=xpsf(i)-dxpsf
          ypsf(i)=ypsf(i)-dypsf
          read(10,*) ttt,ttt
          read(10,*) (ttt,kkk=3,(n+1)*(n+2)/2+2)
          read(10,*)
c          write(0,'(4g15.5)') xpsf(i),ypsf(i),sumpsf(i),bgpsf(i)
      enddo
   11 npsfstars=i-1
      close(10)
      write(0,*) 'Read ',npsfstars,' PSF star positions'

c open the image and error files

      call get2err('inimage.fits','err.fits',nx,ny,pix,pixerr,mmpix,mmpix)

c     read through the list of sources, and process each one:
c     make a shapelet expansion
c     calculate the PSF
c     compute the GaaP fluxes

      betaprev=-1
      do i=1,mmcat
    1     read(*,*,err=1,end=99) x,y,f,fe,fw,a,b,th,id,iflg,beta,qq
          nread=nread+1
          if (beta.ne.betaprev) then
c fit PSF residuals with this scale radius
              filt2=max(2*qq**2+(beta/1.3)**2-2*(betapsf/1.3)**2,0.5)
              do is=1,npsfstars
                  xx=xpsf(is)
                  yy=ypsf(is)
                  call interpolpsf(xx,yy,xmax,ymax,psfcoef,npsf,betapsf,nfit,shpsf)
                  do k=3,(npsf+1)*(npsf+2)/2+2
                      shpsf(k)=shpsf(k)*sumpsf(is)
                  enddo
                  res(is)=(filt2+(betapsf/1.3)**2)/filt2/sumpsf(is) *
     :             pixresid(shpsf,xx,yy,pix,mmpix,nx,ny,bgpsf(is),filt2)
                  wt(is)=1.
              enddo
              nrejiter=3
              clip=4.
              call fitpoly(npsfstars,xmax,ymax,xpsf,ypsf,res,wt,nfit,
     :            nrejiter,clip,psfrescoef,bestfit)
c              write(0,*) (psfrescoef(k),k=1,(nfit+1)*(nfit+2)/2)
              betaprev=beta
c              write(0,*) 'Got PSF resid corr factor map for beta=',beta
          endif
          fitradius=max(10.,4*beta)
          if (x.gt.fitradius+1 .and. x.lt.nx-fitradius-1
     :        .and. y.gt.fitradius+1 .and. y.lt.ny-fitradius-1) then
              call fitbg(pix,mmpix,nx,ny,x,y,fitradius,bg)
              call fitshapeletserr(pix,pixerr,mmpix,nx,ny,x,y,fitradius
     :            ,bg,shobs,sherr,n,beta)
              call interpolpsf(x,y,xmax,ymax,psfcoef,npsf,betapsf,nfit,shpsf)
              fq=psfcorphot3(shobs,sherr,shpsf,qq,fqerr)
c squared dispersion of gaussian wt fn for residuals
              filtdisp2=max(2*qq*qq-(betapsf/1.3)**2,0.5)
              fqresid=(qq*qq/filtdisp2)*
     :            pixresid(shobs,x,y,pix,mmpix,nx,ny,bg,filtdisp2)
c              write(0,*) 'fqresid factor:',qq*qq/filtdisp2,qq,sqrt(filtdisp2),fq,fqresid
c abuse interpolpsf routine to interpolate the correction factor (in 00 coef)
              call interpolpsf(x,y,xmax,ymax,psfrescoef,0,0.,nfit,shpsf)
              psfresidfac=shpsf(3)
              if (abs(psfresidfac).gt.0.15 
     :            .or. abs(fqresid).gt.0.15*abs(fq)) iflg=iflg+4000
              write(*,1001) x,y,f,fe,fw,a,b,th,id,iflg,
     :            beta,qq,fq,fq+fqresid,(fq+fqresid)/(1+psfresidfac),fqerr,bg
 1001         format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i6,i5,
     :            2f10.3,6g13.4)
          else
              goto 1
          endif
          
   99 enddo

      end



          
