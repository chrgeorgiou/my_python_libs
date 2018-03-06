c define shape, make shapelet expansion, make plot of fit vs data.

      parameter (mm012=5*(mm/2)-1,mbox=256)
      real z(mbox,mbox),zerr(mbox,mbox),psf(mbox,mbox),xplot(mbox)
      real z2(2*mbox,mbox),z3(2*mbox,mbox),z4(mbox,mbox)
      integer nn(2)
      real shape(msh),shpsf(msh),sh012(mm012),psfmat(msh-2,msh-2)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :    ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2),coefs(mm/2+4),
     :    covar(mm/2+4,mm/2+4)

      common /shearshiftfns012/ round,sheared1,sheared2,shifted1,shifted2,nshear
c      common /shearfns/ round,sheared1,sheared2,nshear

c      write(*,*) 'n,nsersic,rpsf?'
c      read(*,*) n,nsersic,rpsf

      write(0,*) 'Shapelet order n, PSF type ip?'
      read(*,*) n,ip
      rhmfac=sqrt(2**(2./5)-1)
      nsh=(n+1)*(n+2)/2
      rpct=20.
      reff=3.
      call pgbeg(0,'/xs',4,2)
      call pgask(.false.)
      nsub=10
      do irpsf=2,6
          rpsf=irpsf
          nshear=n
          do ig1=0,14
              g1in=0.05*ig1
              g2in=0
              do nsersic=0,4
                  do i=1,mbox
                      do j=1,mbox
                          zerr(j,i)=1e-3
                      enddo
                  enddo
                  ctr=mbox/2.
                  do i=1,mbox
                      do j=1,mbox
c PIXEL-SMEAR THE PSF IN THE SAME WAY AS THE SHEARED GALAXY
c without smearing the devauc profiles are too peaky.
                          z(i,j)=0
                          psf(i,j)=0
                          do isub=-2*nsub+1,2*nsub-1,2
                              dx=isub/real(4*nsub)
                              x=i-ctr+dx
                              do jsub=-2*nsub+1,2*nsub-1,2
                                  dy=jsub/real(4*nsub)
                                  y=j-ctr+dy
                                  xx=x*(1-g1in)-y*g2in
                                  yy=-x*g2in+y*(1+g1in)
                                  rsh=sqrt(xx*xx+yy*yy)
                                  if (rsh.le.rpct*3) then
                                      if (nsersic.eq.0) then
                                          z(i,j)=z(i,j)+exp(-0.5*(rsh/3.)**2)
                                      else
c                                          z(i,j)=z(i,j)+exp(-log(100.)*(rsh/rpct)**(1./nsersic))
                                          z(i,j)=z(i,j)+exp(-(2*nsersic-0.33)*(rsh/reff)**(1./nsersic))
                                      endif
                                  endif
                                  r=sqrt(x*x+y*y)
                                  re=sqrt((x*1.1)**2+(y/1.1)**2)
                                  if (r.gt.mbox/4) then
                                      psf(i,j)=0
                                  else
                                      if (ip.eq.1) then
                                          psf(i,j)=psf(i,j)+exp(-0.693*(r/rpsf)**2)
                                      elseif (ip.eq.2) then
                                          psf(i,j)=psf(i,j)+(1+(rhmfac*r/rpsf)**2)**(-2.5)
                                      elseif (ip.eq.3) then
                                          psf(i,j)=psf(i,j)+exp(-(0.93*r/rpsf)**2)+exp(-(0.93*r/rpsf/2)**2)/4.
                                      elseif (ip.eq.4) then
                                          psf(i,j)=psf(i,j)+exp(-0.693*r/rpsf)
                                      elseif (ip.eq.5) then
                                          psf(i,j)=psf(i,j)+exp(-0.693*(re/rpsf)**2)
                                      elseif (ip.eq.6) then
                                          psf(i,j)=psf(i,j)+(1+(rhmfac*re/rpsf)**2)**(-2.5)
                                      elseif (ip.eq.7) then
                                          psf(i,j)=psf(i,j)+exp(-(0.93*re/rpsf)**2)+exp(-(0.93*r/rpsf/2)**2)/4.
                                      elseif (ip.eq.8) then
                                          psf(i,j)=psf(i,j)+exp(-(0.93*r/rpsf)**2)+exp(-(0.93*re/rpsf/2)**2)/4.
                                      endif
                                  endif
                              enddo
                          enddo
                      enddo
                  enddo
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
                  
                  call showpatch(psf,mbox,mbox,mbox,mbox,mbox*0.5,mbox*0.5,mbox*0.1,0.,zmx)
                  zmx=z(nint(ctr),nint(ctr))
                  call showpatch(z,mbox,mbox,mbox,mbox,mbox*0.5,mbox*0.5,mbox*0.1,0.,zmx)
                  nn(1)=mbox
                  nn(2)=mbox

                  do j=1,mbox
                      do i=1,mbox
                          z3(2*i-1,j)=psf(i,j)
                          z3(2*i,j)=0
                      enddo
                  enddo
                  call fourn(z3,nn,2,1)
                  do i=1,mbox/2
                      do j=1,mbox/2
                          z4(i,j)=sqrt(z3(2*i-1,j)**2+z3(2*i,j)**2)
                      enddo
                  enddo
                  call showpatch(z4,mbox,mbox,mbox/2,mbox/2,mbox*0.125,mbox*0.125,mbox*0.25,0.,z4(2,2))
                  
                  do j=1,mbox
                      do i=1,mbox
                          z2(2*i-1,j)=z(i,j)
                          z2(2*i,j)=0
                      enddo
                  enddo
                  call fourn(z2,nn,2,1)
                  do i=1,mbox/2
                      do j=1,mbox/2
                          z4(i,j)=sqrt(z2(2*i-1,j)**2+z2(2*i,j)**2)
                      enddo
                  enddo
                  call showpatch(z4,mbox,mbox,mbox/2,mbox/2,mbox*0.125,mbox*0.125,mbox*0.25,0.,z4(2,2))
                  
                  isgn=1
                  do j=1,mbox
                      do i=1,mbox
c     convolution in F-space is multiplication.
c     must also shift PSF by a half-cycle: change the sign of all odd Fourier cpts.
                          prodr=z2(2*i-1,j)*z3(2*i-1,j)-z2(2*i,j)*z3(2*i,j)
                          prodi=z2(2*i-1,j)*z3(2*i,j)+z2(2*i,j)*z3(2*i-1,j)
                          z2(2*i-1,j)=prodr*isgn
                          z2(2*i,j)=prodi*isgn
                          isgn=-isgn
                      enddo
                      isgn=-isgn
                  enddo
                  do i=1,mbox/2
                      do j=1,mbox/2
                          z4(i,j)=sqrt(z2(2*i-1,j)**2+z2(2*i,j)**2)
                      enddo
                  enddo
                  call showpatch(z4,mbox,mbox,mbox/2,mbox/2,mbox*0.125,mbox*0.125,mbox*0.25,0.,z4(2,2))

                  call fourn(z2,nn,2,-1)

                  do j=1,mbox
                      do i=1,mbox
                          z(i,j)=z2(2*i-1,j)
                      enddo
                  enddo
                  
                  zmx=z(1,1)
                  zmn=zmx
                  do i=1,mbox
                      do j=1,mbox
                          zmx=max(zmx,z(i,j))
                          zmn=min(zmn,z(i,j))
                      enddo
                  enddo
                  write(0,*) 'Displaying range ',zmn,zmx
                  call showpatch(z,mbox,mbox,mbox,mbox,mbox*0.5,mbox*0.5,mbox*0.1,zmn,zmx)
                  
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
                      call fitshapeletserr(z,zerr,mbox,mbox,mbox,ctr-1,ctr-1,5*beta,bg,shape,sherr,n,beta)
                      sum=seriessum(shape)
                      do i=3,nsh+2
                          shape(i)=shape(i)/sum
                      enddo
                      call writeshapen('CONVOLVED',shape)
                      call plotshape(shape,mbox/4,'SHAPELET FIT')
                      bg=0.
                      fitradius=max(10.,8*betapsf)
                      call fitshapelets(psf,mbox,mbox,mbox,ctr,ctr,fitradius,bg,shpsf,n,betapsf)
c                      betapsf=shpsf(1)
                      sum=seriessum(shpsf)
                      do i=3,nsh+2
                          shpsf(i)=shpsf(i)/sum
                      enddo
                      call writeshapen('PSF',shpsf)
                      call plotshape(shpsf,mbox/8,'PSF SHAPELET FIT')
c                      call to012(shape,sh012,n)
c                      write(10,'(40g13.5)') beta,(sh012(k),k=1,5*(n/2)-1)
                      betaintr=sqrt(max(1.1**2*betapsf**2,beta**2-betapsf**2))
c     write(0,*) 'Beta, betapsf, betaintr :',beta,betapsf,betaintr
                      call setpsfmat(shpsf,beta,betaintr,psfmat)
                      nshear=n
                      call setroundshearshift_012(psfmat,round,sheared1,sheared2,shifted1,shifted2)
                      call fitshearshift_err_012(shape,sherr,g1,g2,g1sig,g2sig,coefs,covar,.false.)
                      write(*,*) n,nsersic,rpsf,beta,bg,g1in,g1,(coefs(i),i=1,n/2)
                  enddo
              enddo
          enddo
      enddo
      call pgend
      end
