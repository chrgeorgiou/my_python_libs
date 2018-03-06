c define shape, make shapelet expansion, make plot of fit vs data.

      real z(41,41),z2(41,41)
      real shape(msh)
      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
     :    ,shifted1(msh-2,0:mm/2),shifted2(msh-2,0:mm/2),coefs(mm/2+4),
     :    covar(mm/2+4,mm/2+4)

      common /shearshiftfns/ round,sheared1,sheared2,shifted1,shifted2,nshear
c      common /shearfns/ round,sheared1,sheared2,nshear

      write(*,*) 'beta, n?'
      read(*,*) beta,n
c      write(*,*) 'Input shears g1,g2?'
c      read(*,*) g1in,g2in
      g1in=0.1
      g2in=0
      write(*,*) 'Noise level (fraction of peak)?'
      read(*,*) signoise
      write(*,*) 'Number of realizations?'
      read (*,*) niter
      idum=-1
      call pgbeg(0,'/xs',2,2)
      call pgask(.false.)
      do iter=1,niter
          do i=1,41
              x=(i-21.)/3.
              do j=1,41
                  y=(j-21.)/3.
                  xx=x*(1-g1in)-y*g2in
                  yy=-x*g2in+y*(1+g1in)
                  r=sqrt(.1+xx*xx+yy*yy)
                  if (r.gt.4.) then
                      z(i,j)=0
                  else 
                      z(i,j)=exp(-r)
                  endif
                  z(i,j)=z(i,j)+signoise*gasdev(idum)
              enddo
          enddo
          call fitshapelets(z,41,41,41,21.,21.,20.,0.,shape,n,beta)
          if (iter.le.4) call plotshape(shape,41,'SHAPELET FIT')
c          call writeshapen('Exp disk',shape)
c          call decode(shape,z2,41,41)
          
c          call pgenv(-0.5,1.5,-0.5,1.5,1,0)
c          call pgpt(41*41,z,z2,5)
c          call pgmove (-10.,-10.)
c          call pgdraw(10.,10.)
c          call pgenv(0.,30.,-0.1,1.1,0,0)
c          do i=1,41
c              do j=1,41
c                  r=sqrt((i-21.)**2+(j-21.)**2)
c                  call pgpt1(r,z(i,j),5)
c              enddo
c          enddo
c          call pgsci(3)
c          do i=1,41
c              do j=1,41
c                  r=sqrt((i-21.)**2+(j-21.)**2)
c                  call pgpt1(r,z2(i,j),5)
c              enddo
c          enddo
c          call pgsci(1)
          
          nshear=n
          call setroundshearshiftnopsf(beta,n,round,sheared1,sheared2,shifted1,shifted2)
          call fitshearshift_err(shape,0.1,g1,g2,g1sig,g2sig,coefs,covar,.false.)

c          call setroundshearnopsf(beta,n,round,sheared1,sheared2,shifted1,shifted2)
c          call fitshear_err(shape,0.1,g1,g2,g1sig,g2sig,coefs,covar,.false.)
          write(10,*) g1,g2,g1sig,g2sig,covar(n/2+1,n/2+2)/(g1sig*g2sig)
          sg1=sg1+g1
          sg2=sg2+g2
          ssg1=ssg1+g1**2
          ssg2=ssg2+g2**2
      enddo
      write(*,*) sg1/niter, sg2/niter, sqrt(ssg1/niter - (sg1/niter)**2), sqrt(ssg2/niter-(sg2/niter)**2)

      call pgend
      end
