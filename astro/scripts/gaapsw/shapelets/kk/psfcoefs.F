      parameter (mm012=5*(mm/2)-1)
      real sh(msh),psfmean(msh),pk(mm),polarpsf(mm012)


      read(*,*) n,beta
      read(*,*) nfit,xmax,ymax
      read(*,*)
      nsh=(n+1)*(n+2)/2
c just read leading order PSF coeff
      psfmean(1)=beta
      psfmean(2)=n
      do k=1,nsh
         read(*,*) psfmean(k+2),(ttt,kfit=2,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo
      do k=nsh,1,-1
          psfmean(k+2)=psfmean(k+2)/psfmean(3)
      enddo
      call to012(psfmean,polarpsf,n)

c plot coeffs as function of order - linear scale
      call pgbeg(0,'psfcoefs.ps/ps',2,2)
      call pgslw(2)
      call pgsch(1.5)
      call pgenv(-1.,n+1.,-0.5,1.1,0,0)
      call pgmove(-0.5,0.)
      call pgdraw(n+0.5,0.)
      k=0
      do l=0,n
          do m=0,l
              k=k+1
              call pgpt(1,l*1.,psfmean(k+2),23)
          enddo
      enddo
      call pglab('ORDER','NORMALIZED COEFF','PSF COEFFS.')
c plot coeffs as function of order - log scale incl line with total power
      call pgenv(-1.,n+1.,-5.,0.1,0,0)
      k=0
      do l=0,n
          pk(l+1)=0
          do m=0,l
              k=k+1
              p=psfmean(k+2)
              if (p.gt.0.) then
                  call pgpt(1,l*1.,log10(p),23)
              else
                  call pgpt(1,l*1.,log10(-p),5)
              endif
              pk(l+1)=pk(l+1)+p*p
          enddo
          pk(l+1)=sqrt(pk(l+1))
      enddo
      call pgmove(0.,log10(pk(1)))
      do l=1,n
          call pgdraw(l*1.,log10(pk(l+1)))
      enddo
      call pglab('ORDER','log(NORMALIZED COEFF)','PSF COEFFS.')
          

c plot m=0 polar component as function of radial order - linear scale
      call pgenv(-1.,n+1.,-1.,1.,0,0)
      call pgmove(0.,polarpsf(1))
      do i=1,n/2
          call pgdraw(2.*i,polarpsf(1+i))
      enddo
      call pglab('ORDER','NORMALIZED COEFF','m=0')

c plot m=1,2 polar components as function of order - linear scale
      zmax=0
      do k=n/2+2,5*(n/2)-1
          zmax=max(zmax,abs(polarpsf(k)))
      enddo
      call pgenv(-1.,n+1.,-zmax*1.1,zmax*1.1,0,0)
      call pgsls(1)
      call pgmove(1.,polarpsf(n/2+2))
      do i=1,n/2-1
          call pgdraw(1.+2*i,polarpsf(n/2+2+i))
      enddo
      call pgsls(3)
      call pgmove(1.,polarpsf(n+2))
      do i=1,n/2-1
          call pgdraw(1.+2.*i,polarpsf(n+2+i))
      enddo
      call pgsls(1)
      call pgmove(2.,polarpsf(3*n/2+2))
      do i=1,n/2-2
          call pgdraw(2.+2.*i,(-1)**i*polarpsf(3*n/2+2+i))
      enddo
      call pgsls(3)
      call pgmove(2.,polarpsf(2*n+1))
      do i=1,n/2-2
          call pgdraw(2.+2.*i,(-1)**i*polarpsf(2*n+1+i))
      enddo
      call pglab('ORDER','NORMALIZED COEFFS','m=1,2')

      call pgend
      end

