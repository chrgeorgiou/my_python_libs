c test of the shear determination
c what is required?
c
c 1. A PSF decomposed into shapelets, and written as a matrix
c      call fitshapelets(zpsf,50,31,31,16.,16.,rpsf,0.,shpsf,npsf,betapsf)
c
c 2. An average shape to fit as a smeared, sheared cicular source
c      call fitshapelets(z,50,41,41,21.,21.,rfit,0.,shape,n,betaout)
c
c 3. A PSF matrix, adjusted to the proper scales of PSF and object
c      call setpsfmat(shpsf,betaout,beta,pmat)
c
c 4. A set of basis functions representing action of smear and shear
c      call setroundshear(pmat,round,sheared1,sheared2)
c
c 5. A call to fitshear to do all the fitting
c      call fitshear(shape,g1,g2,g1sig,g2sig,coefs,covar,.false.)
c
c That's it!
c Steps 3 and 4 only need to be done once per PSF and beta.


      real zpsf(50,50),shpsf(msh),cmat(msh-2,msh-2,msh-2),z(50,50)
      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
      real pmat(msh-2,msh-2)
      real shape0(msh),shape(msh),sh(msh),sh2(msh),sh3(msh),sh4(msh),
     :     coefs(mm/2+2),covar(mm/2+2,mm/2+2)
      double precision pp
      common /convcoef/ cmat,b1,b2,bout,norder
      common /shearfns/ round,sheared1,sheared2,nshear

c get PSF

      write(0,*) 'PSF: Half-ht radius, Moffat index n ?'
      read(*,*) prh,pn
      pcoef=(2.**(2./pn)-1)/prh**2
      do i=1,31
          x=i-16
          do j=1,31
              y=j-16
              zpsf(i,j)=(1+pcoef*(x*x+y*y))**(-0.5*pn)
          enddo
      enddo

      write(0,*) 'PSF shapelet parameters beta, order ?'
      read(*,*) betapsf,npsf
      rpsf=6*prh
      call fitshapelets(zpsf,50,31,31,16.,16.,rpsf,0.,shpsf,npsf,betapsf)
            
      call pgbeg(0,'/xs',2,2)
      write(0,*) 'Plotting PSF'
      call plotshape(shpsf,21,'PSF')
      call writeshapen('PSF',shpsf)
      
      do i=1,41
          do j=1,41
              z(i,j)=0
          enddo
      enddo

      write(0,*) 'Galaxy shear?'
      read(*,*) g1,g2

      write(*,*) 'Galaxy exponential scale length?'
      read(*,*) gr

      do i=1,41
         x=i-21
         do j=1,41
            y=j-21
            r=sqrt(((1-g1)*x-g2*y)**2+(-g2*x+(1+g1)*y)**2)/gr
            if (r.lt.4.) then
               g=exp(-r)
               do k=max(1,17-i),min(31,47-i)
                  do l=max(1,17-j),min(31,47-j)
                     z(i+k-16,j+l-16)=z(i+k-16,j+l-16)+g*zpsf(k,l)
                  enddo
               enddo
            endif
         enddo
      enddo

      write(0,*) 'Observed galaxy shapelet parameter beta ?'
      read(*,*) beta
      n=shpsf(2)
      betaout=beta
      rfit=6*beta
      call fitshapelets(z,50,41,41,21.,21.,rfit,0.,shape,n,betaout)
      write(0,*) 'Plotting convolved image to fit'
      call plotshape(shape,21,'Smeared')
      call writeshapen('Smeared',shape)

      call setpsfmat(shpsf,betaout,beta,pmat)

      nshear=norder
      call setroundshear(pmat,round,sheared1,sheared2)
      call fitshear(shape,g1,g2,g1sig,g2sig,coefs,covar,.true.)
      write(0,'("Shear g1, g2, errors:",4f10.4)') g1,g2,g1sig,g2sig

      sh(1)=shape(1)
      sh(2)=shape(2)
      n=sh(2)
      do k=3,(n+1)*(n+2)/2+2
         sh(k)=0
         do m=0,n-2,2
            sh(k)=sh(k)+coefs(m/2+1)*round(k-2,m/2)
         enddo
      enddo
      write(0,*) 'Plotting unsheared, convolved shape'
      call plotshape(sh,21,'Unsheared x PSF')
      call writeshapen('Unsheared x PSF',sh)

      call pgend
      end
      










c      write(0,*) 'testing funcs:'
      
c      do i=1,nshear/2+2
c         coefs(i)=1.4+i
c      enddo
c      call funcs(4.,coefs,y,der)
c      do i=1,nshear/2+2
c         coefs(i)=coefs(i)+0.1
c         call funcs(4.,coefs,y2,der2)
c         coefs(i)=coefs(i)-0.1
c         write(0,*) i,(y2-y)/0.1,der(i)
c      enddo
c      write(0,*)

