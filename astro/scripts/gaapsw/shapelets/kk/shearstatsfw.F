
c from a shear estimator catalogue (from sh2g) estimate a weighted
c mean of the shear.
c weight is gaussian, depends on total g of each object
c iterate this to avoid an aperture bias

      real fw(mmcat),g1(mmcat),g2(mmcat),w(mmcat),eg1(mmcat),eg2(mmcat)
      double precision s,s1,s2
      logical more,rpt
      character lbl*33

      snmin=10
      open(10,file='sn.dat',status='old',err=80)
      read(10,*) snmin
      close(10)

   80 write(0,*) 'Minimum S/N is ',snmin

      psf10=0
      psf01=0
      psf20=0
      psf11=0
      psf02=0
c read the PSF size
      open(10,file='psf.map',status='old')
      read(10,*) npsf,betapsf
      close(10)
c only include unflagged objects with fwhm above 1.1 x that of the PSF 
c and S/N above the limit
      nread=0
      do i=1,mmcat
    1     read(*,*,err=1,end=99) x,y,f,fe,fw(i),a,b,pa,
     :        g1(i),g2(i),eg1(i),eg2(i),cov,id,iflg
          nread=nread+1
          if (fw(i).lt.betapsf*2.3*1.1) goto 1
          if (iflg.ne.0) goto 1
          if (f.lt.fe*snmin) goto 1
      enddo
   99 ndat=i-1
      write(0,*) 'Including ',ndat,'/',nread,' shear estimators.'
      do i=1,ndat
          w(i)=g1(i)*g1(i)+g2(i)*g2(i)
      enddo
c take median g as the weight radius. Not far from sigma implied by fwhm
      rwt2=fmedian(ndat,w)/2.
      write(0,'('' Weight radius for g: '',f6.3)') sqrt(rwt2)
      rpt=.false.
      do i=1,ndat
          w(i)=exp(-0.5*w(i)/rwt2)
      enddo
c sort by flux, bin flux, and make averages of g1 and g2 per bin
      call sort4(ndat,fw,g1,g2,w)
      do ibin=0,9
          s=0
          s1=0
          s2=0
          i1=(ndat*ibin)/10+1
          i2=(ndat*(ibin+1))/10
          do i=i1,i2
              s=s+w(i)
              s1=s1+w(i)*g1(i)
              s2=s2+w(i)*g2(i)
          enddo
          write(*,'(2g15.5,4f10.5)') fw(i1),fw(i2),s1/s,s2/s,fmedian(i2-i1+1,g1(i1)),fmedian(i2-i1+1,g2(i1))
      enddo
      end
