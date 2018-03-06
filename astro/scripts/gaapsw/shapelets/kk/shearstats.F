
c from a shear estimator catalogue (from sh2g) estimate a weighted
c mean of the shear.
c weight is gaussian, depends on total g of each object
c iterate this to avoid an aperture bias

      parameter (maxcat=10*mmcat)
      real g1(maxcat),g2(maxcat),w(maxcat),eg1(maxcat),eg2(maxcat)
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
c read the leading coeffs of the PSF expansion, for reporting later
      open(10,file='psf.map',status='old')
      read(10,*) npsf,betapsf
      read(10,*) nfit,xmax,ymax
      read(10,*) psf00,(t,i=2,(nfit+1)*(nfit+2)/2)
      read(10,*,end=90) psf10,(t,i=2,(nfit+1)*(nfit+2)/2)
      read(10,*,end=90) psf01,(t,i=2,(nfit+1)*(nfit+2)/2)
      read(10,*,end=90) psf20,(t,i=2,(nfit+1)*(nfit+2)/2)
      read(10,*,end=90) psf11,(t,i=2,(nfit+1)*(nfit+2)/2)
      read(10,*,end=90) psf02,(t,i=2,(nfit+1)*(nfit+2)/2)
   90 close(10)
c only include unflagged objects with fwhm above 1.1 x that of the PSF 
c and S/N above the limit
      nread=0
      do i=1,maxcat
    1     read(*,*,err=1,end=99) x,y,f,fe,fw,a,b,pa,
     :        g1(i),g2(i),eg1(i),eg2(i),cov,id,iflg
          nread=nread+1
          if (fw.lt.betapsf*2.3*1.1 .or. iflg.ne.0 .or. f.lt.fe*snmin) then
              goto 1
          endif
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
      avg1=0
      avg2=0
c iterate the weighted mean for g1 and g2
    5 s1=0
      s2=0
      s=0
      do i=1,ndat
          s1=s1+g1(i)*w(i)
          s2=s2+g2(i)*w(i)
          s=s+w(i)
      enddo
      avg1new=s1/s
      avg2new=s2/s
      if (.not.rpt) then
          write(0,'('' g1, g2 mean:  '',2f10.4)') avg1new,avg2new
          rpt=.true.
      endif
      more=(abs(avg1new-avg1).gt.1e-5 .or. abs(avg2new-avg2).gt.1e-5)
      if (more) then
          avg1=avg1new
          avg2=avg2new
          do i=1,ndat
              w(i)=exp(-0.5/rwt2*((g1(i)-avg1)**2+(g2(i)-avg2)**2))
          enddo
          goto 5
      endif
      write(0,'('' g1, g2 mean:  '',2f10.4)') avg1new,avg2new

c now calculate the shear variance in the field as the weighted avg of 
c g1(a)g1(b) + g2(a)*g2(b), a!=b
c for this use a centrally symmetric weight function.

      s1=0
      s2=0
      s=0
      do i=1,ndat
          e2=eg1(i)**2+eg2(i)**2
          e2=max(1e-10,e2)
          w(i)=max(1.,0.04/e2)
      enddo
      do i=1,ndat
          do j=i+1,ndat
              s1=s1+g1(i)*g1(j)*w(i)*w(j)
              s2=s2+g2(i)*g2(j)*w(i)*w(j)
              s=s+w(i)*w(j)
          enddo
      enddo
      g1var=s1/s
      g2var=s2/s
      g1med=fmedian(ndat,g1)
      g2med=fmedian(ndat,g2)
      write(0,'('' g1, g2 median:'',2f10.4)') g1med,g2med
      write(*,'(4f10.5,2g15.5,f10.5,4f10.5,i8)') avg1,avg2,g1med,g2med,g1var,g2var,betapsf,psf00,psf20,psf11,psf02,ndat

      call pgbeg(0,'shears.ps/cps',1,1)
      call pgslw(3)
      call pgenv(-0.8,0.8,-0.8,0.8,1,1)
      call pgpt(ndat,g1,g2,5)
      write(lbl,'(i6,'' FILTERED OBJECTS S/N>'',f5.1)') ndat,snmin
      call pglab('G1','G2',lbl)
      call pgend
      write(0,*) 'Plot is in shears.ps'
      end
