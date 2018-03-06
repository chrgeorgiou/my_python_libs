c     calculate the shear correlation function over quadrants, splitting
c     the field into nested square grids. Only use cells in which there
c     are at least 10 galaxies. 

      parameter (maxcat=10*mmcat)
      real g1(maxcat),g2(maxcat),w(maxcat),x(maxcat),y(maxcat),g1var(256)
     $    ,g2var(256),gg1(maxcat),gg2(maxcat),ww(maxcat)
      double precision s,s1,s2,s11,s22
      logical more,rpt
      character lbl*33

      nread=0
      do i=1,maxcat
    1     read(*,*,err=1,end=99) x(i),y(i),g1(i),g2(i),t,t,t,amag,e2
          if (e2.gt.0.7) goto 1
          if (amag.lt.17.) goto 1
          if (amag.gt.21.) goto 1
          nread=nread+1
          e2=e2*e2
          w(i)=1/e2
      enddo
   99 ndat=i-1
      write(0,*) 'Including ',ndat,'/',nread,' shear estimators.'

c now calculate the shear variance in the field as the weighted avg of 
c g1(a)g1(b) + g2(a)*g2(b), a!=b
c for this use a centrally symmetric weight function.

      nwinfull=2000
      do idiv=1,8
          nwin=nwinfull/idiv
          ncell=0
          do ix=0,2000/nwin-1
              do iy=0,2000/nwin-1
                  incell=0
                  do i=1,ndat
                      if (int(x(i)/nwin).eq.ix .and. int(y(i)/nwin).eq.iy) then
                          incell=incell+1
                          gg1(incell)=g1(i)
                          gg2(incell)=g2(i)
                          ww(incell)=w(i)
                      endif
                  enddo
                  if (incell.ge.10) then
                      ncell=ncell+1
                      s1=0
                      s2=0
                      s=0
                      do i=1,incell
                          do j=i+1,incell
                              s1=s1+gg1(i)*gg1(j)*ww(i)*ww(j)
                              s2=s2+gg2(i)*gg2(j)*ww(i)*ww(j)
                              s=s+ww(i)*ww(j)
                          enddo
                      enddo
                      g1var(ncell)=s1/s
                      g2var(ncell)=s2/s
                  endif
              enddo
          enddo
          if (ncell.gt.0) then
              s1=0
              s2=0
              s11=0
              s22=0
              do i=1,ncell
                  s1=s1+g1var(i)
                  s2=s2+g2var(i)
                  s11=s11+g1var(i)**2
                  s22=s22+g2var(i)**2
              enddo
              write(*,'(i3,4g15.5,i5)') idiv,s1/ncell,s2/ncell,s11/ncell,s22/ncell,ncell
          endif
      enddo
      end
