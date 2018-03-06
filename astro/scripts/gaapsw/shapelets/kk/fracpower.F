      real sh(msh),f0(mmcat),f1(mmcat),f2(mmcat),f3(mmcat),f4(mmcat),f5(mmcat),
     :    f6(mmcat),f7(mmcat),f8(mmcat),flx(mmcat),fw(mmcat),sherr(mmcat),sn(mmcat)
      real p(0:mm),frac(2:6)
      integer iflg(mmcat)

      amagzero=35
      open(10,file='magzero.dat',status='old',err=2)
      read(10,*) amagzero
      close(10)
    2 write(0,'('' Magnitude zero point:''f6.2)') amagzero

      snmin=10
      open(10,file='sn.dat',status='old',err=80)
      read(10,*) snmin
      close(10)
   80 write(0,*) 'Minimum S/N is ',snmin

      shiftmax=9
      do i=2,6
          frac(i)=1
      enddo
      open(10,file='shapecuts.dat',status='old',err=3)
      read(10,*) shiftmax,(frac(i),i=2,6)
      close(10)
    3 write(0,'('' Shift, shape (2..6) cuts: '',6f6.3)') 
     :    shiftmax,(frac(i),i=2,6)

      do i=1,mmcat
    1     read(*,*,end=99) x,y,flx(i),fe,fw(i),a,b,pa,id,iflg(i)
          read(*,*) sum,dx,dy,bg,rfit,sherr(i)
          read(*,*) beta,n
          read(*,*) (sh(k),k=3,(n+1)*(n+2)/2+2)
          read(*,*)
          sn(i)=log10(flx(i)/fe)
          if (flx(i).lt.snmin*fe) iflg(i)=1000+iflg(i)
          sherr(i)=log10(sherr(i))
          f1(i)=dx*dx+dy*dy
          flx(i)=amagzero-2.5*log10(flx(i))
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
          f0(i)=log10(p(0)/sump)
          f1(i)=log10(f1(i))/2.
          f2(i)=log10(p(2)/sump)
          f3(i)=log10(p(3)/sump)
          f4(i)=log10(p(4)/sump)
          f5(i)=log10(p(5)/sump)
          f6(i)=log10(p(6)/sump)
          f7(i)=log10(p(7)/sump)
      enddo
   99 nobj=i-1
      call pgbeg(0,'fracpower-mag.ps/cps',3,3)
      call pgsch(2.)
      call pgslw(2)
      call pgenv(15.,29.,-4.,0.1,0,20)
      call pgpt(nobj,flx,f0,-1)
      call pglab('Mag','Fraction','0')
      call pgenv(15.,29.,-3.,1.,0,20)
      call pgpt(nobj,flx,f1,-1)
      call plothor(15.,29.,log10(shiftmax))
      call pglab('Mag','DR','1')
      call pgenv(15.,29.,-4.,0.1,0,20)
      call pgpt(nobj,flx,f2,-1)
      call plothor(15.,29.,log10(frac(2)))
      call pglab('Mag','Fraction','2')
      call pgenv(15.,29.,-4.,0.1,0,20)
      call pgpt(nobj,flx,f3,-1)
      call plothor(15.,29.,log10(frac(3)))
      call pglab('Mag','Fraction','3')
      call pgenv(15.,29.,-4.,0.1,0,20)
      call pgpt(nobj,flx,f4,-1)
      call plothor(15.,29.,log10(frac(4)))
      call pglab('Mag','Fraction','4')
      call pgenv(15.,29.,-4.,0.1,0,20)
      call pgpt(nobj,flx,f5,-1)
      call plothor(15.,29.,log10(frac(5)))
      call pglab('Mag','Fraction','5')
      call pgenv(15.,29.,-4.,0.1,0,20)
      call pgpt(nobj,flx,f6,-1)
      call plothor(15.,29.,log10(frac(6)))
      call pglab('Mag','Fraction','6')
      call pgenv(15.,29.,-4.,0.1,0,20)
      call pgpt(nobj,flx,f7,-1)
      call pglab('Mag','Fraction','7')
      call pgend

      call pgbeg(0,'fracpower-fwhm.ps/cps',3,3)
      call pgsch(2.)
      call pgslw(2)
      call pgenv(2.,12.,-4.,0.1,0,20)
      call pgpt(nobj,fw,f0,-1)
      call pglab('FWHM','Fraction','0')
      call pgenv(2.,12.,-3.,1.,0,20)
      call pgpt(nobj,fw,f1,-1)
      call plothor(2.,12.,log10(shiftmax))
      call pglab('FWHM','DR','1')
      call pgenv(2.,12.,-4.,0.1,0,20)
      call pgpt(nobj,fw,f2,-1)
      call plothor(2.,12.,log10(frac(2)))
      call pglab('FWHM','Fraction','2')
      call pgenv(2.,12.,-4.,0.1,0,20)
      call pgpt(nobj,fw,f3,-1)
      call plothor(2.,12.,log10(frac(3)))
      call pglab('FWHM','Fraction','3')
      call pgenv(2.,12.,-4.,0.1,0,20)
      call pgpt(nobj,fw,f4,-1)
      call plothor(2.,12.,log10(frac(4)))
      call pglab('FWHM','Fraction','4')
      call pgenv(2.,12.,-4.,0.1,0,20)
      call pgpt(nobj,fw,f5,-1)
      call plothor(2.,12.,log10(frac(5)))
      call pglab('FWHM','Fraction','5')
      call pgenv(2.,12.,-4.,0.1,0,20)
      call pgpt(nobj,fw,f6,-1)
      call plothor(2.,12.,log10(frac(6)))
      call pglab('FWHM','Fraction','6')
      call pgenv(2.,12.,-4.,0.1,0,20)
      call pgpt(nobj,fw,f7,-1)
      call pglab('FWHM','Fraction','7')
      call pgend

      call pgbeg(0,'fracpower-sherr.ps/cps',3,3)
      call pgsch(2.)
      call pgslw(2)
      call pgenv(-4.,0.,-4.,0.1,0,30)
      call pgpt(nobj,sherr,f0,-1)
      call pglab('SHERR','Fraction','0')
      call pgenv(-4.,0.,-3.,1.,0,30)
      call pgpt(nobj,sherr,f1,-1)
      call plothor(-4.,0.,log10(shiftmax))
      call pglab('SHERR','DR','1')
      call pgenv(-4.,0.,-4.,0.1,0,30)
      call pgpt(nobj,sherr,f2,-1)
      call plothor(-4.,0.,log10(frac(2)))
      call pglab('SHERR','Fraction','2')
      call pgenv(-4.,0.,-4.,0.1,0,30)
      call pgpt(nobj,sherr,f3,-1)
      call plothor(-4.,0.,log10(frac(3)))
      call pglab('SHERR','Fraction','3')
      call pgenv(-4.,0.,-4.,0.1,0,30)
      call pgpt(nobj,sherr,f4,-1)
      call plothor(-4.,0.,log10(frac(4)))
      call pglab('SHERR','Fraction','4')
      call pgenv(-4.,0.,-4.,0.1,0,30)
      call pgpt(nobj,sherr,f5,-1)
      call plothor(-4.,0.,log10(frac(5)))
      call pglab('SHERR','Fraction','5')
      call pgenv(-4.,0.,-4.,0.1,0,30)
      call pgpt(nobj,sherr,f6,-1)
      call plothor(-4.,0.,log10(frac(6)))
      call pglab('SHERR','Fraction','6')
      call pgenv(-4.,0.,-4.,0.1,0,30)
      call pgpt(nobj,sherr,f7,-1)
      call pglab('SHERR','Fraction','7')
      call pgend

      call pgbeg(0,'fracpower-sn.ps/cps',3,3)
      call pgsch(2.)
      call pgslw(2)
      call pgenv(0.45,3.,-4.,0.1,0,30)
      call pgpt(nobj,sn,f0,-1)
      call pglab('S/N','Fraction','0')
      call pgenv(0.45,3.,-3.,1.,0,30)
      call pgpt(nobj,sn,f1,-1)
      call plothor(0.45,3.,log10(shiftmax))
      call plotvert(log10(snmin),-3.,1.)
      call pglab('S/N','DR','1')
      call pgenv(0.45,3.,-4.,0.1,0,30)
      call pgpt(nobj,sn,f2,-1)
      call plothor(0.45,3.,log10(frac(2)))
      call plotvert(log10(snmin),-4.,0.)
      call pglab('S/N','Fraction','2')
      call pgenv(0.45,3.,-4.,0.1,0,30)
      call pgpt(nobj,sn,f3,-1)
      call plothor(0.45,3.,log10(frac(3)))
      call plotvert(log10(snmin),-4.,0.)
      call pglab('S/N','Fraction','3')
      call pgenv(0.45,3.,-4.,0.1,0,30)
      call pgpt(nobj,sn,f4,-1)
      call plothor(0.45,3.,log10(frac(4)))
      call plotvert(log10(snmin),-4.,0.)
      call pglab('S/N','Fraction','4')
      call pgenv(0.45,3.,-4.,0.1,0,30)
      call pgpt(nobj,sn,f5,-1)
      call plothor(0.45,3.,log10(frac(5)))
      call plotvert(log10(snmin),-4.,0.)
      call pglab('S/N','Fraction','5')
      call pgenv(0.45,3.,-4.,0.1,0,30)
      call pgpt(nobj,sn,f6,-1)
      call plothor(0.45,3.,log10(frac(6)))
      call plotvert(log10(snmin),-4.,0.)
      call pglab('S/N','Fraction','6')
      call pgenv(0.45,3.,-4.,0.1,0,30)
      call pgpt(nobj,sn,f7,-1)
      call plotvert(log10(snmin),-4.,0.)
      call pglab('S/N','Fraction','7')

      call pgend

      write(0,*) 'Plots are in fracpower-*.ps'
      end


      subroutine plotsegment(x1,y1,x2,y2,icol)
      call pgqci(icolsave)
      call pgsci(icol)
      call pgmove(x1,y1)
      call pgdraw(x2,y2)
      call pgsci(icolsave)
      return
      end

      subroutine plothor(x1,x2,y)
      call plotsegment(x1,y,x2,y,2)
      return
      end

      subroutine plotvert(x,y1,y2)
      call plotsegment(x,y1,x,y2,2)
      return
      end
