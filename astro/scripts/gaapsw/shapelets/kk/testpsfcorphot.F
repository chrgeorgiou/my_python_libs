c     set up a gaussian source, stddev sqrt(a^2+p^2), and a PSF with
c     stddev p expand in shapelets with scale bs and bp
c
c     construct PSF matrix P that convolves a source with scale q with a
c     PSF of scale bp to yield a result with scale bs
c
c     invert P, apply to source to give deconvolved source = gaussian of
c     scale a, expressed with scale bq
c
c     construct new PSF = gaussian of scale and dispersion q
c
c     construct new PSF matrix that convolves a source of scale q and a
c     PSF of scale q to a result of scale q
c
c     read off the 00 component of the result, which is the
c     gaussian-weighted aperture photometry of the source convolved with
c     gaussian PSF of scale q.

      dimension z(101,101),psf(101,101)
      dimension shobs(msh),shpsf(msh)
      parameter (pi=3.1415926535)

c a,p,q are gaussian dispersion of intrinsic gal, PSF, target PSF
      a=4
      p=3
      betaobs=6
      betapsf=4
      q=5
      n=8
      t=2
      aa=a*a+p*p/t
      bb=a*a+p*p*t
      nsh=(n+1)*(n+2)/2
      do i=1,101
         x=i-51
         do j=1,101
            y=j-51
            r2=x*x+y*y
            z(i,j)=exp(-(x*x/aa+y*y/bb)/2)/(2*pi*sqrt(aa*bb))
            psf(i,j)=exp(-(t*x*x+y*y/t)/2/(p*p))/(2*pi*p*p)
         enddo
      enddo
      sg=0
      sp=0
      do i=1,101
         do j=1,101
            sg=sg+z(i,j)
            sp=sp+psf(i,j)
         enddo
      enddo
      write(0,*) 'Integral, sig of galaxy and PSF:',sg,a,sp,p

      call fitshapelets(z,101,101,101,51.,51.,4*betaobs,0.,shobs,n,betaobs)
      call writeshape('OBS GAL',shobs)
      call fitshapelets(psf,101,101,101,51.,51.,4*betapsf,0.,shpsf,n,betapsf)
      call writeshape('PSF',shpsf)

      do kkk=1,1000,100
         qq=q+0.001*kkk
         phot=psfcorphot(shobs,shpsf,qq)
         write(0,*) 'Exper, theor phot00 values:',phot,qq*qq/(a*a+2*qq*qq)
      enddo

      end
