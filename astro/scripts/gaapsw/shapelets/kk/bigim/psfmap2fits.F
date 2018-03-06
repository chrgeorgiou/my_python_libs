      character psfmapname*100,psfname*20
      parameter (mz=128)
      real sh(msh),psfcoef(mmpoly2,msh-2),z(mz,mz)

      write(0,*) 'File with psf map coefs?'
      read(*,'(a)') psfmapname
      open(17,file=psfmapname,status='old')
      read(17,*) n,beta
      read(17,*) nfit,xmax,ymax
      read(17,*)
      do k=1,(n+1)*(n+2)/2
         read(17,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(17,*)
      enddo
      close(17)
      write(0,*) 'Read '//psfmapname

      do i=1,9999
         write(0,*) 'PSF x & y coordinates, size of image to generate?'
         read(*,*,end=10) xfield,yfield,nz
         write(psfname,'(''psf_'',i5.5,''_'',i5.5,''.fits'')') 
     :        int(xfield),int(yfield)
         call interpolpsf(xfield,yfield,xmax,ymax,psfcoef,n,beta,nfit,sh)
         xc=nz/2+xfield-int(xfield)
         yc=nz/2+yfield-int(yfield)
         call decodectr(sh,z,mz,nz,xc,yc)
         call put2(psfname,nz,nz,z,mz,mz)
         write(0,*) 'Wrote fits file '//psfname
      enddo
 10   end

