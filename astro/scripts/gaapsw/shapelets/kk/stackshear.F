      real sh(msh),sumsh(msh),shpsf(msh),sumshpsf(msh),
     :    psfcoef(mmpoly2,msh-2),psfmat(msh-2,msh-2)
      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
      real coefs(mm/2+2),covar(mm/2+2,mm/2+2)
      real p(0:mm),frac(2:6)
      common /shearfns/ round,sheared1,sheared2,nshear


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

c set minimum S/N for inclusion
      snmin=10
      open(10,file='sn.dat',status='old',err=80)
      read(10,*) snmin
      close(10)
   80 write(0,*) 'Minimum S/N is ',snmin

c set cuts on maximum fraction of power that can be found at order 2-6
      shiftmax=9
      do i=2,6
          frac(i)=1
      enddo
      open(10,file='shapecuts.dat',status='old',err=3)
      read(10,*) shiftmax,(frac(i),i=2,6)
      close(10)
    3 write(0,'('' Shift, shape (2..6) cuts: '',6f6.3)') 
     :    shiftmax,(frac(i),i=2,6)

      nread=0
      beta=0
      nstacked=0
      avg1=0
      avg2=0
      ntot=0
      call pgbeg(0,'stackshear.ps/vps',4,5)
      call pgsch(2.)
      call pgslw(3)
      do i=1,mmcat
    1     read(*,*,err=1,end=99) x,y,f,fe,fw,a,b,th,id,iflg
          read(*,*) sum,dx,dy,bg,fitradius,sherr
          read(*,*) sh(1),n
          sh(2)=n
          read(*,*) (sh(k),k=3,(n+1)*(n+2)/2+2)
          read(*,*)
c apply filters to the catalogue - set flags if shape looks dubious
c penalty flag is 100 for each shape constraint violated; 
c penalty is 1000 if fwhm is smaller than 1.1*PSF fwhm
          k=2
          sump=0
          do nn=0,n
              p(nn)=0
              do j=0,nn
                  k=k+1
                  p(nn)=p(nn)+sh(k)**2
              enddo
              sump=sump+p(nn)
          enddo
          if (dx*dx+dy*dy .gt. shiftmax*shiftmax) iflg=iflg+100
          do nn=2,min(n,6)
              if (p(nn).gt.frac(nn)*sump) iflg=iflg+100
          enddo
          if (fw.lt.1.1*betapsf*2.3) iflg=iflg+1000
c only include unflagged objects with sufficiently high S/N
          if (iflg.ne.0) goto 1
          if (f.lt.fe*snmin) goto 1
          nread=nread+1
c          write(10,'(i8,8f10.2)') id,f/fe,fw/(betapsf*2.3),(dx*dx+dy*dy)/shiftmax**2,(p(nn)/frac(nn)/sump,nn=2,6)
          call interpolpsf(x,y,xmax,ymax,psfcoef,npsf,betapsf,nfit,shpsf)
c check whether the shapelet has the same beta as the stack.
          if (sh(1).ne.beta) then
              if (nstacked.gt.0) then
c do a shear fit
                  do k=3,(n+1)*(n+2)/2+2
                      sumshpsf(k)=sumshpsf(k)/nstacked
                  enddo
                  call setpsfmat(sumshpsf,beta,beta,psfmat)
                  nshear=n
                  call setroundshear(psfmat,round,sheared1,sheared2)
                  call fitshear(sumsh,g1,g2,g1sig,g2sig,coefs,covar,.false.)
                  write(*,'(f8.4,i8,2f8.4)') beta,nstacked,g1,g2
                  call plotshape(sumsh,21,'SUM')
                  call plotshape(sumshpsf,21,'PSF')
                  avg1=avg1+nstacked*g1
                  avg2=avg2+nstacked*g2
                  ntot=ntot+nstacked
              endif
              nstacked=1
              do k=1,(n+1)*(n+2)/2+2
                  sumsh(k)=sh(k)
                  sumshpsf(k)=shpsf(k)
              enddo
              beta=sumsh(1)
          else
              do k=1,(n+1)*(n+2)/2
                  sumsh(k+2)=sumsh(k+2)+sh(k+2)
                  sumshpsf(k+2)=sumshpsf(k+2)+shpsf(k+2)
              enddo
              nstacked=nstacked+1
          endif
      enddo
   99 write(0,*) 'Read ',nread,' objects.'
      if (ntot.gt.0) then
          avg1=avg1/ntot
          avg2=avg2/ntot
      endif
      write(*,'(f8.4,i8,2f8.4)') 0.,ntot,avg1,avg2
      call pgend
      end
