      parameter (mc=10)
      real cov(-mc:mc,-mc:mc)
      call openfits('inimage.fits',51,nx,ny)
      sig= guess_errors_sec(51,nx*ny)
      call guesscov_sec(51,nx,ny,mc,cov,0.001*sig)
      end
