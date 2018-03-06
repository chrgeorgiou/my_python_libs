c define shape, make shapelet expansion, make plot of fit vs data.

      parameter (mm012=5*(mm/2)-1,mbox=256)
      real z(mbox,mbox),zerr(mbox,mbox),psf(mbox,mbox),xplot(mbox)
      integer nn(2)
      real xser(1000000),yser(1000000),xpsf(1000000),ypsf(1000000)
      character ch,ch3*3
      real shape(msh),shpsf(msh),sh012(mm012),psfmat(msh-2,msh-2)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :    ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2),coefs(mm/2+4),
     :    covar(mm/2+4,mm/2+4)
      logical plots

      common /shearshiftfns012/ round,sheared1,sheared2,shifted1,shifted2,nshear
c      common /shearfns/ round,sheared1,sheared2,nshear

c      write(*,*) 'n,nsersic,rpsf?'
c      read(*,*) n,nsersic,rpsf


      plots=.true.
      idum=-1
      write(0,*) 'Shapelet order n, PSF Moffat index?'
      read(*,*) n,bmoffat
      write(ch3,'(f3.1)') bmoffat
      open(15,file='moffat_'//ch3,status='old')
    1 read (15,*,err=1) (xpsf(i),ypsf(i),i=1,1000000)
      close(15)
      nsh=(n+1)*(n+2)/2
      if (plots) call pgbeg(0,'/xs',4,2)
      if (plots) call pgask(.false.)
      rpsf=6
      do i=1,mbox
          do j=1,mbox
              z(j,i)=0
              psf(j,i)=0
          enddo
      enddo
      ctr=mbox/2.
c     start MC building of PSF and source * PSF use 2 difft passes of
c     sources, starting at difft pt in list. make sure all pts are 8-pt
c     symmetric (no m=1,2)
      do i=1,8000000
c     write(0,*) i,mod(i+ioff,1000000)+1
          i8=(i-1)/8
          xxp=xpsf(mod(i8,1000000)+1)*rpsf
          yyp=ypsf(mod(i8,1000000)+1)*rpsf
          if (mod(i,8).eq.0) then
              xp= xxp
              yp= yyp
          elseif (mod(i,8).eq.1) then
              xp= yyp
              yp=-xxp
          elseif (mod(i,8).eq.2) then
              xp=-xxp
              yp=-yyp
          elseif (mod(i,8).eq.3) then
              xp=-yyp
              yp= xxp
          elseif (mod(i,8).eq.4) then
              xp=-xxp
              yp= yyp
          elseif (mod(i,8).eq.5) then
              xp= yyp
              yp= xxp
          elseif (mod(i,8).eq.6) then
              xp= xxp
              yp=-yyp
          elseif (mod(i,8).eq.7) then
              xp=-yyp
              yp=-xxp
          endif
c     pixellate
          ii=nint(ctr+xp)
          jj=nint(ctr+yp)
          if (ii.le.mbox .and.ii.ge.1.and.jj.le.mbox.and.jj.ge.1) psf(ii,jj)=psf(ii,jj)+1
      enddo
      write(0,*) 'MC-convolved with PSF'
      do i=2,-2,-1
          write(0,*) (z(nint(ctr)+i,nint(ctr)+j),j=-2,2)
      enddo

c measure FWHM of PSF; set beta accordingly
      zmx=psf(nint(ctr),nint(ctr))
      ic=nint(ctr)
      jc=nint(ctr)
      i=1
      do while (psf(ic+i,jc).gt.zmx*0.5)
c     write(0,*) ic+i,z(ic+i,jc)
          i=i+1
      enddo
      rhm1=i-1 + (psf(ic+i-1,jc)-zmx*0.5) / (psf(ic+i-1,jc)-psf(ic+i,jc))
      j=1
      do while (psf(ic,jc+j).gt.zmx*0.5)
c     write(0,*) jc+j,z(ic,jc+j)
          j=j+1
      enddo
      rhm2=j-1 + (psf(ic,jc+j-1)-zmx*0.5) / (psf(ic,jc+j-1)-psf(ic,jc+j))
      
      fwhm=2*sqrt(rhm1*rhm2)
      betapsf=fwhm/2.3
      write(0,*) 'Fitted PSF beta is ',betapsf
                  
      zmx=psf(nint(ctr),nint(ctr))
      if (plots) call showpatch(psf,mbox,mbox,mbox,mbox,mbox*0.5,mbox*0.5,mbox*0.1,0.,zmx)
      write(0,*) 'Displaying range ',0.,zmx
      
c fit shapelet exp with result of fitbeta
      bg=0.
      fitradius=max(10.,8*betapsf)
      call fitbeta(psf,mbox,mbox,mbox,ctr,ctr,fitradius,0.,betapsf)
      fitradius=max(10.,8*betapsf)
      call fitshapelets(psf,mbox,mbox,mbox,ctr,ctr,fitradius,bg,shpsf,n,betapsf)
c     betapsf=shpsf(1)
      call writeshapen('PSF',shpsf)
      if (plots) call plotshape(shpsf,mbox/8,'PSF SHAPELET FIT')

      do ix=0,0
          call decode(shpsf,z,mbox,99)
          do i=1,99
              write(12,*) i-50.-ix*0.1,z(i,50)
          enddo
      enddo

c now repeat fit with beta*1.3
      betapsf=betapsf*1.3
      fitradius=max(10.,8*betapsf)
      call fitshapelets(psf,mbox,mbox,mbox,ctr,ctr,fitradius,bg,shpsf,n,betapsf)
c     betapsf=shpsf(1)
      call writeshapen('PSF',shpsf)
      if (plots) call plotshape(shpsf,mbox/8,'PSF SHAPELET FIT')

      do ix=0,0
          call decode(shpsf,z,mbox,99)
          do i=1,99
              write(13,*) i-50.-ix*0.1,z(i,50)
          enddo
      enddo

      do i=mbox/2-49,mbox/2+49
          write(14,*) i-mbox/2,psf(i,mbox/2)
      enddo


      if (plots) call pgend
      end
