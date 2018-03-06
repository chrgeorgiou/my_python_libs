      call openfits('inimage.fits',51,nx,ny)
      call copyfits(51,61,'new.fits',5)
      call fitscomment(61,'test comment')
      call fitshistory(61,'PSF Gaussianized')
      call fitswritekey(61,'GPSFSIG',7.43,'SIGMA OF GAUSSIANIZED PSF')
      call closefits(61)
      call closefits(51)
      end
