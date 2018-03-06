c read in a PSF map and subtract a perfect Gaussian from it - expanded in the scale radius of the supplied psf map.

c input is 
c one line containing the dispersion of the Gaussian to subtract
c then the lines of a PSF map file.

      parameter (mpix=800)
      real zz((2*mpix+1)*(2*mpix+1))
      real sh(msh),psfcoef(mmpoly2,msh-2)

c get sigma of target gaussian
      read(*,*) sigtarget
      npix=min(mpix,int(10*sigtarget+5))
c read PSF map
      read(*,*) npsf,betapsf
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(npsf+1)*(npsf+2)/2
         read(*,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo
      write(0,*) 'Read psf model:'
      write(0,*) 'Order for shape, spatial variation; beta:',npsf,nfit,betapsf
      do ish=1,(npsf+1)*(npsf+2)/2
         sh(ish+2)=psfcoef(1,ish)
      enddo
      sh(1)=betapsf
      sh(2)=npsf
      write(0,*) 'Average PSF component integral is ',seriessum(sh)
c calculate the shapelet coefs for the target gaussian      
      call gaus2shpl(sigtarget,sh,betapsf,npsf,npix,zz)
c subtract this gaussian from the PSF map and write out.
      do ish=1,(npsf+1)*(npsf+2)/2
         psfcoef(1,ish)=psfcoef(1,ish)-sh(ish+2)
      enddo
      
      write(*,*) npsf,betapsf
      write(*,*) nfit,xmax,ymax
      write(*,*)
      do ish=1,(npsf+1)*(npsf+2)/2
         write(*,*) (psfcoef(kfit,ish),kfit=1,(nfit+1)*(nfit+2)/2)
         write(*,*)
      enddo
      end


 




      subroutine gaus2shpl(sig,sh,beta,n,npix,zz)
c express a single gaussian of dispersion sig as a shapelet series (of arbitrary scale radios beta and to order n)      
c do this by constructing a pixel grid and fitting shapelet series
c gaussian is centered on a grid with coords -npix..npix
      real zz(-npix:npix,-npix:npix),sh(msh),sh2(msh)

      nsh=(n+1)*(n+2)/2+2
      npix2=2*npix+1

c for four different pixel phases, construct a gaussian pixel grid and
c calculate the shapelet expansion for it. Average the results

c center cpt
      do j=-npix,npix
         do i=-npix,npix
            r2=(i*i+j*j)/sig**2
            if (r2.lt.30.) then
               zz(i,j)=exp(-0.5*r2)
            else
               zz(i,j)=0
            endif
         enddo
      enddo
      call fitshapelets(zz,npix2,npix2,npix2,npix+1.,npix+1.,npix*1.,1.e-31,sh,n,beta)
      sum=seriessum(sh)*4
      do ish=3,nsh
         sh(ish)=sh(ish)/sum
      enddo
c x offset 0.5 pix
      do j=-npix,npix
         do i=-npix,npix
            r2=((i-0.5)*(i-0.5)+j*j)/sig**2
            if (r2.lt.30.) then
               zz(i,j)=exp(-0.5*r2)
            else
               zz(i,j)=0
            endif
         enddo
      enddo
      call fitshapelets(zz,npix2,npix2,npix2,npix+1.5,npix+1.,npix*1.,1.e-31,sh2,n,beta)
      sum=seriessum(sh2)*4
      do ish=3,nsh
         sh(ish)=sh(ish)+sh2(ish)/sum
      enddo
c y offset 0.5 pix
      do j=-npix,npix
         do i=-npix,npix
            r2=(i*i+(j-0.5)*(j-0.5))/sig**2
            if (r2.lt.30.) then
               zz(i,j)=exp(-0.5*r2)
            else
               zz(i,j)=0
            endif
         enddo
      enddo
      call fitshapelets(zz,npix2,npix2,npix2,npix+1.,npix+1.5,npix*1.,1.e-31,sh2,n,beta)
      sum=seriessum(sh2)*4
      do ish=3,nsh
         sh(ish)=sh(ish)+sh2(ish)/sum
      enddo
c x,y offset 0.5 pix each
      do j=-npix,npix
         do i=-npix,npix
            r2=((i-0.5)*(i-0.5)+(j-0.5)*(j-0.5))/sig**2
            if (r2.lt.30.) then
               zz(i,j)=exp(-0.5*r2)
            else
               zz(i,j)=0
            endif
         enddo
      enddo
      call fitshapelets(zz,npix2,npix2,npix2,npix+1.5,npix+1.5,npix*1.,1.e-31,sh2,n,beta)
      sum=seriessum(sh2)*4
      do ish=3,nsh
         sh(ish)=sh(ish)+sh2(ish)/sum
      enddo
c zero out all odd components
      kk=3
      do k=0,n
         do i=0,k
            if (mod(k,2).eq.1 .or. mod(i,2).eq.1) sh(kk)=0
            kk=kk+1
         enddo
      enddo
c      write(0,*) (sh(k),k=1,nsh)
c      write(0,*) 'sum of target gaussian shapelet is ',seriessum(sh)
      return
      end
