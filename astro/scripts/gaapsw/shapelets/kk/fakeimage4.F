      parameter (mmpsf=500)
      real z(mmpix,mmpix),zerr(mmpix,mmpix),psf(-mmpsf:mmpsf,-mmpsf:mmpsf)
      double precision sum
      parameter (satur=50000., nover=3)

c makes realization of PSF stars and elliptical, randomly oriented galaxies.
c stars are moffat profiles, galaxies exponentials.

c read parameters for PSF stars, and shape of the galaxy
      read(*,*) idum,nx,ny,fs,pa,pb,ppa,pmoffat,ns,fg,ga,gq,g1,g2,ng,bg,idx,idy,f2
      write(0,'('' Image size to make, bg: '',2i6,g9.2)') nx,ny,bg
      write(0,'('' Stars: flux, rms maj, rms min, PA, Moffat idx, no.:'',
     :    g9.2,2f6.2,f6.0,f6.2,i6)') fs,pa,pb,ppa,pmoffat,ns
      write(0,'('' Gals:  flux, exp scl, q, g1, g2, no.:'',
     :    g9.2,4f6.2,i6)') fg,ga,gq,g1,g2,ng
      write(0,'('' PSF asymmetry: '',2f5.1,f5.2)') idx/real(nover),idy/real(nover),f2
      do i=1,mmpix
          do j=1,mmpix
              z(i,j)=bg
          enddo
      enddo
      cpa=cos(ppa*3.14159265/180)
      spa=sin(ppa*3.14159265/180)
      npsf=int(max(pa,pb)*nover) * 10.**(4./pmoffat)+1
      if (npsf.gt.mmpsf) then
         write(*,*) 'PSF array too small - need ',npsf
         stop
      endif
      pa=pa/sqrt(2.**(2./pmoffat)-1)/2
      pb=pb/sqrt(2.**(2./pmoffat)-1)/2
c set psf radius as level where it hits 10^-4
      do i=-mmpsf,mmpsf
         do j=-mmpsf,mmpsf
            psf(i,j)=0
         enddo
      enddo
      do i=-npsf,npsf
         x=i/real(nover)
         do j=-npsf,npsf
            y=j/real(nover)
            xp= x*cpa+y*spa
            yp=-x*spa+y*cpa
            ppp=(1+(xp/pa)**2+(yp/pb)**2)**(-pmoffat/2)
            psf(i,j)=ppp
            psf(i+idx,j+idy)=psf(i+idx,j+idy) + f2*ppp
            sum=sum+psf(i,j)
         enddo
      enddo
      sum=0
      do i=-npsf,npsf
         do j=-npsf,npsf
            sum=sum+psf(i,j)
         enddo
      enddo
      do i=-npsf,npsf
         do j=-npsf,npsf
            psf(i,j)=psf(i,j)/sum
         enddo
      enddo


      do istar=1,ns
         xc=ran3(idum)*nx
         yc=ran3(idum)*ny
         ixc=int(xc)
         iyc=int(yc)
         ixc2=int(nover*(xc-ixc))
         iyc2=int(nover*(yc-iyc))
         flx=fs*(0.01+0.99*ran3(idum))
         do i=-npsf,npsf
            ix=(1000+i-ixc2)/nover-1000/nover+ixc
            if (ix.lt.1 .or. ix.gt.nx) goto 2
            do j=-npsf,npsf
               iy=(1000+j-iyc2)/nover-1000/nover+iyc
               if (iy.lt.1 .or. iy.gt.ny) goto 3
               z(ix,iy)=z(ix,iy)+flx*psf(i,j)
 3          enddo
 2       enddo
         write(*,*) xc,yc,flx,pa,pb,ppa,istar
      enddo
      write(0,*) 'Added ',ns,' stars.'


      do igal=1,ng
         xc=ran3(idum)*nx
         yc=ran3(idum)*ny
         ixc=int(xc)
         iyc=int(yc)
         ixc2=int(nover*(xc-ixc))
         iyc2=int(nover*(yc-iyc))
         flx=fg*(1+ran3(idum))/2.
         gpa=ran3(idum)*3.1415926535*2
         cgpa=cos(gpa)
         sgpa=sin(gpa)
         npix=int(5*ga*nover)+1
         do k=-npix,npix
            x=k/real(nover)
            do l=-npix,npix
               y=l/real(nover)
               xp= x*cgpa+y*sgpa
               yp=-x*sgpa+y*cgpa
               ggg=flx*exp(-abs(xp/ga)-abs(yp/ga/gq))/(nover**2*ga**2*3.1416*gq)
               do i=-npsf,npsf
                  ix=(1000+i+k-ixc2)/nover-1000/nover+ixc
                  if (ix.lt.1 .or. ix.gt.nx) goto 4
                  do j=-npsf,npsf
                     iy=(1000+j+l-iyc2)/nover-1000/nover+iyc
                     if (iy.lt.1 .or. iy.gt.ny) goto 5
                     z(ix,iy)=z(ix,iy)+ggg*psf(i,j)
 5                enddo
 4             enddo
            enddo
         enddo
         write(*,*) xc,yc,flx,gpa*180/3.1415926535,ga,ga*gq,igal+ns
      enddo
      write(0,*) 'Added ',ng,' galaxies.'

      do i=1,nx
          do j=1,ny
              zerr(i,j)=sqrt(z(i,j))
              z(i,j)=min(satur,z(i,j)+gasdev(idum)*zerr(i,j))
          enddo
      enddo

      call put2('fake.fits',nx,ny,z,mmpix,mmpix)
      call put2('err.fits',nx,ny,zerr,mmpix,mmpix)
      end
