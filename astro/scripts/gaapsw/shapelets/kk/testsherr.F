      real sh(mmcat,msh-2),sherr(mmcat),f(mmcat)

      fmin= 1e10
      fmax=-1e10
      do i=1,mmcat
    1     read(*,*,err=1,end=2) x,y,f(i),fe,fw,a,b,pa,id,iflg
          fmin=min(f(i),fmin)
          fmax=max(f(i),fmax)
          read(*,*) sum,dx,dy,bg,fitradius,sherr(i)
          read(*,*) beta,n
          read(*,*) (sh(i,k-2),k=3,(n+1)*(n+2)/2+2)
          read(*,*)
      enddo
    2 nobj=i-1
      write(0,*) 'Read ',nobj,' objects.'
      call pgbeg(0,'/xs',7,7)
      do k=1,(n+1)*(n+2)/2
          sum0=0
          sum1=0
          sum2=0
          do i=1,nobj
              sum0=sum0+1/sherr(i)**2
              sum1=sum1+sh(i,k)/sherr(i)**2
          enddo
          avg=sum1/sum0
          do i=1,nobj
              sum2=sum2+(sh(i,k)-avg)**2/sherr(i)**2
          enddo
          sig=sqrt(sum2/nobj)
          write(0,*) k,avg,sig
          call pgenv(fmin,fmax,-10.,10.,0,0)
          do i=1,nobj
              call pgpt1(f(i),max(-9.,min(9.,(sh(i,k)-avg)/sherr(i))),23)
          enddo
      enddo
      call pgend
      end
