c define shape, make shapelet expansion, make plot of fit vs data.

      parameter (mm012=5*(mm/2)-1,mbox=256)
      real z(mbox,mbox),zerr(mbox,mbox),psf(mbox,mbox),xplot(mbox)
      integer nn(2)
      real xser(1000000),yser(1000000),xpsf(1000000),ypsf(1000000)
      character ch,ch3*3
      real shape(msh),shpsf(msh),sh012(mm012),psfmat(msh-2,msh-2)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :    ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2),coefs(mm/2+4),
     :     sheared11(mm012,0:mm/2),sheared12(mm012,0:mm/2),
     :     sheared22(mm012,0:mm/2),
     :     sheared111(mm012,0:mm/2),sheared112(mm012,0:mm/2),
     :     sheared122(mm012,0:mm/2),sheared222(mm012,0:mm/2),
     :    covar(mm/2+4,mm/2+4)
      logical plots

      common /shearshiftfns012_o3/ round,sheared1,sheared2,shifted1,shifted2,
     :     sheared11,sheared12,sheared22,
     :     sheared111,sheared112,sheared122,sheared222,nshear
c      common /shearfns/ round,sheared1,sheared2,nshear

c      write(*,*) 'n,nsersic,rpsf?'
c      read(*,*) n,nsersic,rpsf


      plots=.false.
      idum=-1
      write(0,*) 'Shapelet order n, PSF Moffat index, PSF type (0=rnd, 1=lop, 2=ell), # g-steps?'
      read(*,*) n,bmoffat,iasym, nstep
      write(ch3,'(f3.1)') bmoffat
      open(15,file='moffat_'//ch3,status='old')
    1 read (15,*,err=1) (xpsf(i),ypsf(i),i=1,1000000)
      nsh=(n+1)*(n+2)/2
      reff=4.
      if (plots) call pgbeg(0,'/xs',4,2)
      if (plots) call pgask(.false.)
      nsub=10
      do nsersic=0,4
          write(ch,'(i1)') nsersic
          open(15,file='sersic_'//ch,status='old')
    2     read (15,*,err=2) (xser(i),yser(i),i=1,1000000)
          write(0,*) 'Read sersic_'//ch
          do irpsf=2,6
              rpsf=irpsf
              nshear=n
              do ig1=0,nstep
                  g1in=0.05*ig1
                  g2in=0
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
                  do i=1,16000000
                      if (mod(i,8000000).eq.1) ioff=int(1000000*ran3(idum))
c                      write(0,*) i,mod(i+ioff,1000000)+1
                      i8=(i-1)/8
                      xxg=xser(mod(i8+ioff,1000000)+1)*reff
                      yyg=yser(mod(i8+ioff,1000000)+1)*reff
                      xxp=xpsf(mod(i8,1000000)+1)*rpsf
                      yyp=ypsf(mod(i8,1000000)+1)*rpsf
                      if (mod(i,8).eq.0) then
                          xg= xxg
                          yg= yyg
                          xp= xxp
                          yp= yyp
                      elseif (mod(i,8).eq.1) then
                          xg= yyg
                          yg=-xxg
                          xp= yyp
                          yp=-xxp
                      elseif (mod(i,8).eq.2) then
                          xg=-xxg
                          yg=-yyg
                          xp=-xxp
                          yp=-yyp
                      elseif (mod(i,8).eq.3) then
                          xg=-yyg
                          yg= xxg
                          xp=-yyp
                          yp= xxp
                      elseif (mod(i,8).eq.4) then
                          xg=-xxg
                          yg= yyg
                          xp=-xxp
                          yp= yyp
                      elseif (mod(i,8).eq.5) then
                          xg= yyg
                          yg= xxg
                          xp= yyp
                          yp= xxp
                      elseif (mod(i,8).eq.6) then
                          xg= xxg
                          yg=-yyg
                          xp= xxp
                          yp=-yyp
                      elseif (mod(i,8).eq.7) then
                          xg=-yyg
                          yg=-xxg
                          xp=-yyp
                          yp=-xxp
                      endif
                      if (iasym.eq.2) then
c     elliptical PSF:
                          xp=xp/1.1
                          yp=yp*1.1
                      elseif (iasym.eq.1) then
c     lopsided PSF:
                          if (ran3(idum).lt.0.25) then
                              xp=xp-rpsf
                          else
                              xp=xp+rpsf/3.
                          endif
                      endif
c     shear galaxy and add PSF 
                      xx=xg*(1+g1in)+yg*g2in  + xp
                      yy=xg*g2in+yg*(1-g1in) + yp
c     pixellate
                      ii=nint(ctr+xx)
                      jj=nint(ctr+yy)
                      if (ii.le.mbox .and.ii.ge.1.and.jj.le.mbox.and.jj.ge.1) z(ii,jj)=z(ii,jj)+1
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
                  
                  if (plots) call showpatch(psf,mbox,mbox,mbox,mbox,mbox*0.5,mbox*0.5,mbox*0.1,0.,zmx)
                  zmx=z(nint(ctr),nint(ctr))
                  write(0,*) 'Displaying range ',0.,zmx
                  if (plots) call showpatch(z,mbox,mbox,mbox,mbox,mbox*0.5,mbox*0.5,mbox*0.1,0.,zmx)
                  ic=nint(ctr)-1
                  jc=nint(ctr)-1
                  i=1
                  do while (z(ic+i,jc).gt.zmx*0.5)
c     write(0,*) ic+i,z(ic+i,jc)
                      i=i+1
                  enddo
                  rhm1=i-1 + (z(ic+i-1,jc)-zmx*0.5) / (z(ic+i-1,jc)-z(ic+i,jc))
                  j=1
                  do while (z(ic,jc+j).gt.zmx*0.5)
c     write(0,*) jc+j,z(ic,jc+j)
                      j=j+1
                  enddo
                  rhm2=j-1 + (z(ic,jc+j-1)-zmx*0.5) / (z(ic,jc+j-1)-z(ic,jc+j))
                  
                  fwhm=2*sqrt(rhm1*rhm2)
                  betafit=fwhm/2.3
                  write(0,*) 'Fitted beta is ',betafit
                  
                  do i=1,mbox
                      do j=1,mbox
                          zerr(i,j)=zmx/100.
                      enddo
                  enddo
                  
                  do ibeta=0,0
                      beta=betafit*1.15**(ibeta)
c                      call fitbg(z,mbox,mbox,mbox,ctr-1,ctr-1,10*beta,bg)
                      bg=0.
                      fitradius=max(10.,5*beta)
                      call fitbeta(z,mbox,mbox,mbox,ctr,ctr,fitradius,0.,beta)
                      beta=beta*1.3
                      call fitshapeletserr(z,zerr,mbox,mbox,mbox,ctr,ctr,5*beta,bg,shape,sherr,n,beta)
                      sum=seriessum(shape)
                      do i=3,nsh+2
                          shape(i)=shape(i)/sum
                      enddo
                      call writeshapen('CONVOLVED',shape)
                      if (plots) call plotshape(shape,mbox/4,'SHAPELET FIT')
                      bg=0.
                      fitradius=max(10.,8*betapsf)
                      call fitbeta(psf,mbox,mbox,mbox,ctr,ctr,fitradius,0.,betapsf)
                      betapsf=betapsf*1.3
                      fitradius=max(10.,8*betapsf)
                      call fitshapelets(psf,mbox,mbox,mbox,ctr,ctr,fitradius,bg,shpsf,n,betapsf)
c                      betapsf=shpsf(1)
                      sum=seriessum(shpsf)
                      do i=3,nsh+2
                          shpsf(i)=shpsf(i)/sum
                      enddo
                      call writeshapen('PSF',shpsf)
                      if (plots) call plotshape(shpsf,mbox/8,'PSF SHAPELET FIT')
c                      call to012(shape,sh012,n)
c                      write(10,'(40g13.5)') beta,(sh012(k),k=1,5*(n/2)-1)
                      betaintr=sqrt(max(0.21*betapsf**2,beta**2-betapsf**2))
                      betaintr=betapsf
c                      betaintr=beta
c     write(0,*) 'Beta, betapsf, betaintr :',beta,betapsf,betaintr
                      call setpsfmat(shpsf,beta,betaintr,psfmat)
                      nshear=n
                      call setroundshearshift_012_o3(psfmat,round,sheared1,sheared2,shifted1,shifted2,sheared11,sheared12,sheared22,sheared111,sheared112,sheared122,sheared222)
                      norder=n
                      call fitshearshift_err_012_o3(shape,sherr,g1,g2,g1sig,g2sig,coefs,covar,.false.)
                      write(*,'(2i4,2f8.4,10g15.5)') n,nsersic,rpsf,beta,bg,g1in,g1,(coefs(i),i=1,n/2),g2
                  enddo
              enddo
          enddo
      enddo
      if (plots) call pgend
      end
