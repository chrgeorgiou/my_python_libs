      real kercoef(mmpoly2,msh-2),sumkercoef(mmpoly2,msh-2)
      character kerfile*132
      logical readpsfmap

      do k=1,msh-2
          do kfit=1,mmpoly2
              sumkercoef(kfit,k)=0
          enddo
      enddo
c initialize min/max values in order to check for differences at the end
      xm0=1e30
      xm1=-xm0
      ym0=1e30
      ym1=-ym0
      nmax=0
      nfitmax=0
      beta0=1e30
      beta1=-beta0
      nsum=0

    1 read(*,'(a132)',err=1,end=2) kerfile
      if (readpsfmap(kerfile,kercoef,n,beta,nfit,xm,ym)) then
          nsum=nsum+1
          xm0=min(xm,xm0)
          xm1=max(xm,xm1)
          ym0=min(xm,xm0)
          ym1=max(xm,xm1)
          nmax=max(n,nmax)
          nfitmax=max(nfit,nfitmax)
          beta0=min(beta,beta0)
          beta1=max(beta,beta1)
          do k=1,(n+1)*(n+2)/2
              do kfit=1,(nfit+1)*(nfit+2)/2
                  sumkercoef(kfit,k)=sumkercoef(kfit,k)+kercoef(kfit,k)
              enddo
          enddo
          write(0,*) 'Added in map ',kerfile
      endif
      goto 1

    2 do k=1,(n+1)*(n+2)/2
          do kfit=1,(nfit+1)*(nfit+2)/2
              sumkercoef(kfit,k)=sumkercoef(kfit,k)/nsum
          enddo
      enddo

      if (xm0.ne.xm1) then
          write(0,*) '*** WARNING: nonzero range in x size of maps: ',xm0,xm1
      endif
      if (ym0.ne.ym1) then
          write(0,*) '*** WARNING: nonzero range in y size of maps: ',ym0,ym1
      endif
      if (beta0.ne.beta1) then
          write(0,*) '*** WARNING: nonzero range in beta of maps: ',beta0,beta1
      endif


      n=nmax
      nfit=nfitmax
      write(*,*) n,beta
      write(*,*) nfit,xm,ym
      write(*,*) 
      do k=1,(n+1)*(n+2)/2
         write(*,*) (sumkercoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         write(*,*)
      enddo
      write(0,*) 'Summed ',nsum,' maps.'
      end
