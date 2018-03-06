
c     fit the ndat data points (x,y,tofit) with polynomials in x and y
c     up to total order norder. Weights w are given as input; 
c     for clipped data w is set to zero at output. 
c     xm and ym are the factors by which to rescale x and y, 
c     niter the number of sigma-clipping iterations (at clip sigmas).  
c     output is the list of coefficients coefs, in order
c     00, 10, 01, 20, 11, 02, ...
c     and a list of best-fit values for the data points.


      subroutine fitpoly(ndat,xm,ym,x,y,tofit,w,norder,niter,clip,coefs,bestfit)
      real x(ndat),y(ndat),tofit(ndat),bestfit(ndat),w(ndat)
      real coefs(mmpoly2),mat(mmpoly2,mmpoly2),model(mmpoly2)
      double precision matd(mmpoly2,mmpoly2),coefsd(mmpoly2)
      real resid(ndat)
      logical more
      integer indx(mmpoly2)
      double precision sum

      nmod=(norder+1)*(norder+2)/2

      do iter=1,niter
          do kk=1,nmod
              do ll=1,nmod
                  matd(kk,ll)=0
              enddo
              coefsd(kk)=0
          enddo
          do i=1,ndat
              kk=0
              do k=0,norder
                  do m=0,k
                      l=k-m
                      kk=kk+1
                      model(kk)=polymodel(x(i),y(i),xm,ym,l,m)
                  enddo
              enddo
              do kk=1,nmod
                  do ll=1,nmod
                      matd(kk,ll)=matd(kk,ll)+model(kk)*model(ll)*w(i)
                  enddo
                  coefsd(kk)=coefsd(kk)+model(kk)*tofit(i)*w(i)
              enddo
          enddo
          do kk=1,nmod
              do ll=1,nmod
                  mat(kk,ll)=matd(kk,ll)
              enddo
              coefs(kk)=coefsd(kk)
          enddo
          call ludcmp(mat,nmod,mmpoly2,indx,det)
          call lubksb(mat,nmod,mmpoly2,indx,coefs)
          more=.false.
          do i=1,ndat
              kk=0
              do k=0,norder
                  do m=0,k
                      l=k-m
                      kk=kk+1
                      model(kk)=polymodel(x(i),y(i),xm,ym,l,m)
                  enddo
              enddo
              sum=0
              do kk=1,nmod
                  sum=sum+model(kk)*coefs(kk)
              enddo
              bestfit(i)=sum
              resid(i)=tofit(i)-bestfit(i)
          enddo
c calculate median and mean abs dev of the residuals
          residmed=fmedian(ndat,resid)
          sum0=0
          sum1=0
          do i=1,ndat
              sum0=sum0+w(i)
              sum1=sum1+abs(resid(i)-residmed)*w(i)
          enddo
          avgdev=sum1/sum0
          do i=1,ndat
              if (abs(resid(i)-residmed).gt.clip*avgdev .and. w(i).gt.0.) then 
                  w(i)=0
                  more=.true.
              endif
          enddo
          if (.not.more) goto 20
      enddo
   20 continue
      return
      end


      function polymodel(x,y,xm,ym,l,m)
c simple power of x and y, with range 0..xm remapped to -1..1
      if (l.eq.0) then
          px=1
      else
          px=(2*x/xm-1)**l
      endif
      if (m.eq.0) then
          py=1
      else
          py=(2*y/ym-1)**m
      endif
      polymodel=px * py
      return
      end


      subroutine interpolpsf(x,y,xm,ym,coef,n,beta,nfit,sh)
      real coef(mmpoly2,msh-2),sh(msh)
      sh(1)=beta
      sh(2)=n
      do k=1,(n+1)*(n+2)/2
          sum=0
          kfit=0
          do i=0,nfit
              do m=0,i
                  l=i-m
                  kfit=kfit+1
                  sum=sum+coef(kfit,k)*polymodel(x,y,xm,ym,l,m)
              enddo
          enddo
          sh(k+2)=sum
      enddo
      return
      end




c     fit the ndat data points (x,y,tofit) with polynomials in x and y
c     up to total order norder. Weights w are given as input; 
c     for clipped data w is set to zero at output. 
c     xm and ym are the factors by which to rescale x and y, 
c     niter the number of sigma-clipping iterations (at clip sigmas).  
c     output is the list of coefficients coefs, in order
c     00, 10, 01, 20, 11, 02, ...
c     and a list of best-fit values for the data points.


      subroutine fitcheby(ndat,xm,ym,x,y,tofit,w,norder,niter,clip,coefs,bestfit)
      real x(ndat),y(ndat),tofit(ndat),bestfit(ndat),w(ndat)
      real coefs(mmpoly2),mat(mmpoly2,mmpoly2),model(mmpoly2)
      double precision matd(mmpoly2,mmpoly2),coefsd(mmpoly2)
      real resid(ndat)
      logical more
      integer indx(mmpoly2)
      double precision sum

      nmod=(norder+1)*(norder+2)/2

      do iter=1,niter
          do kk=1,nmod
              do ll=1,nmod
                  matd(kk,ll)=0
              enddo
              coefsd(kk)=0
          enddo
          do i=1,ndat
              kk=0
              do k=0,norder
                  do m=0,k
                      l=k-m
                      kk=kk+1
                      model(kk)=chebmodel(x(i),y(i),xm,ym,l,m)
                  enddo
              enddo
              do kk=1,nmod
                  do ll=1,nmod
                      matd(kk,ll)=matd(kk,ll)+model(kk)*model(ll)*w(i)
                  enddo
                  coefsd(kk)=coefsd(kk)+model(kk)*tofit(i)*w(i)
              enddo
          enddo
          do kk=1,nmod
              do ll=1,nmod
                  mat(kk,ll)=matd(kk,ll)
              enddo
              coefs(kk)=coefsd(kk)
          enddo
          call ludcmp(mat,nmod,mmpoly2,indx,det)
          call lubksb(mat,nmod,mmpoly2,indx,coefs)
          more=.false.
          do i=1,ndat
              kk=0
              do k=0,norder
                  do m=0,k
                      l=k-m
                      kk=kk+1
                      model(kk)=chebmodel(x(i),y(i),xm,ym,l,m)
                  enddo
              enddo
              sum=0
              do kk=1,nmod
                  sum=sum+model(kk)*coefs(kk)
              enddo
              bestfit(i)=sum
              resid(i)=tofit(i)-bestfit(i)
          enddo
c calculate median and mean abs dev of the residuals
          residmed=fmedian(ndat,resid)
          sum0=0
          sum1=0
          do i=1,ndat
              sum0=sum0+w(i)
              sum1=sum1+abs(resid(i)-residmed)*w(i)
          enddo
          avgdev=sum1/sum0
          do i=1,ndat
              if (abs(resid(i)-residmed).gt.clip*avgdev .and. w(i).gt.0.) then 
                  w(i)=0
                  more=.true.
              endif
          enddo
          if (.not.more) goto 20
      enddo
   20 continue
      return
      end


      function chebmodel(x,y,xm,ym,l,m)
c cheby polynomials of x and y, with range 0..xm remapped to -1..1
      if (l.eq.0) then
          px=1
      else
         t=(2*x/xm-1)
         px=cheby(t,l)
      endif
      if (m.eq.0) then
          py=1
      else
         t=(2*y/ym-1)
         py=cheby(t,m)
      endif
      chebmodel=px * py
      return
      end


      function cheby(t,n)
      if (n.le.0) then
         cheby=1
      elseif (n.eq.1) then
         cheby=t
      else
         tnm2=1
         tnm1=t
         do k=0,n-2
            tn=2*t*tnm1-tnm2
            tnm2=tnm1
            tnm1=tn
         enddo
         cheby=tn
      endif
      return
      end





      subroutine interpolchebpsf(x,y,xm,ym,coef,n,beta,nfit,sh)
      real coef(mmpoly2,msh-2),sh(msh)
      sh(1)=beta
      sh(2)=n
      do k=1,(n+1)*(n+2)/2
          sum=0
          kfit=0
          do i=0,nfit
              do m=0,i
                  l=i-m
                  kfit=kfit+1
                  sum=sum+coef(kfit,k)*chebmodel(x,y,xm,ym,l,m)
              enddo
          enddo
          sh(k+2)=sum
      enddo
      return
      end




c     fit the ndat data points (x,y,tofit) with polynomials in x and y
c     up to total order norder. Weights w are given as input; 
c     for clipped data w is set to zero at output. 
c     xm and ym are the factors by which to rescale x and y, 
c     niter the number of sigma-clipping iterations (at clip sigmas).  
c     output is the list of coefficients coefs, in order
c     00, 10, 01, 20, 11, 02, ...
c     and a list of best-fit values for the data points.


      subroutine fitlegendre(ndat,xm,ym,x,y,tofit,w,norder,niter,clip,coefs,bestfit)
      real x(ndat),y(ndat),tofit(ndat),bestfit(ndat),w(ndat)
      real coefs(mmpoly2),mat(mmpoly2,mmpoly2),model(mmpoly2)
      double precision matd(mmpoly2,mmpoly2),coefsd(mmpoly2)
      real legendremodel
      real resid(ndat)
      logical more
      integer indx(mmpoly2)
      double precision sum

      nmod=(norder+1)*(norder+2)/2

      do iter=1,niter
          do kk=1,nmod
              do ll=1,nmod
                  matd(kk,ll)=0
              enddo
              coefsd(kk)=0
          enddo
          do i=1,ndat
              kk=0
              do k=0,norder
                  do m=0,k
                      l=k-m
                      kk=kk+1
                      model(kk)=legendremodel(x(i),y(i),xm,ym,l,m)
                  enddo
              enddo
              do kk=1,nmod
                  do ll=1,nmod
                      matd(kk,ll)=matd(kk,ll)+model(kk)*model(ll)*w(i)
                  enddo
                  coefsd(kk)=coefsd(kk)+model(kk)*tofit(i)*w(i)
              enddo
c              write(19,*) i,x(i),y(i),tofit(i),w(i),(model(k),k=1,nmod)
c              write(19,*) 'matrix (d.p.): ',((matd(ii,jj),jj=1,nmod),ii=1,nmod)
c              write(19,*) 'vector (d.p.): ',(coefsd(ii),ii=1,nmod)
          enddo
          do kk=1,nmod
              do ll=1,nmod
                  mat(kk,ll)=matd(kk,ll)
              enddo
              coefs(kk)=coefsd(kk)
          enddo
c          write(19,*) 'Solving system iteration ',iter
c          write(19,*) 'matrix: ',((mat(i,j),j=1,nmod),i=1,nmod)
c          write(19,*) 'vector: ',(coefs(i),i=1,nmod)
          call ludcmp(mat,nmod,mmpoly2,indx,det)
          call lubksb(mat,nmod,mmpoly2,indx,coefs)
c          write(19,*) '...done'
c          write(19,*)
          more=.false.
          do i=1,ndat
              kk=0
              do k=0,norder
                  do m=0,k
                      l=k-m
                      kk=kk+1
                      model(kk)=legendremodel(x(i),y(i),xm,ym,l,m)
                  enddo
              enddo
              sum=0
              do kk=1,nmod
                  sum=sum+model(kk)*coefs(kk)
              enddo
              bestfit(i)=sum
              resid(i)=tofit(i)-bestfit(i)
          enddo
c calculate median and mean abs dev of the residuals
          residmed=fmedian(ndat,resid)
          sum0=0
          sum1=0
          do i=1,ndat
              sum0=sum0+w(i)
              sum1=sum1+abs(resid(i)-residmed)*w(i)
          enddo
          avgdev=sum1/sum0
          do i=1,ndat
              if (abs(resid(i)-residmed).gt.clip*avgdev .and. w(i).gt.0.) then 
                  w(i)=0
                  more=.true.
              endif
          enddo
          if (.not.more) goto 20
      enddo
   20 continue
      return
      end


      real function legendremodel(x,y,xm,ym,l,m)
      real legendre
c cheby polynomials of x and y, with range 0..xm remapped to -1..1
      if (l.eq.0) then
          px=1
      else
         t=(2*x/xm-1)
         px=legendre(t,l)
      endif
      if (m.eq.0) then
          py=1
      else
         t=(2*y/ym-1)
         py=legendre(t,m)
      endif
      legendremodel=px * py
      return
      end


      real function legendre(t,n)
      if (n.le.0) then
         legendre=1
      elseif (n.eq.1) then
         legendre=t
      else
         tnm2=1
         tnm1=t
         do k=1,n-1
            tn=((2*k+1)*t*tnm1-k*tnm2)/(k+1)
            tnm2=tnm1
            tnm1=tn
         enddo
         legendre=tn
      endif
      return
      end





      subroutine interpollegendrepsf(x,y,xm,ym,coef,n,beta,nfit,sh)
      real coef(mmpoly2,msh-2),sh(msh)
      real legendremodel
      sh(1)=beta
      sh(2)=n
      do k=1,(n+1)*(n+2)/2
          sum=0
          kfit=0
          do i=0,nfit
              do m=0,i
                  l=i-m
                  kfit=kfit+1
                  sum=sum+coef(kfit,k)*legendremodel(x,y,xm,ym,l,m)
              enddo
          enddo
          sh(k+2)=sum
      enddo
      return
      end
