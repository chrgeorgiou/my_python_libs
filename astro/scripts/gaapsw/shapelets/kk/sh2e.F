c sh2e reads in a shapelets catalogue, and outputs id, e1, e2 values.

c ellipticities are KSB polarizations using beta as weight radius.


      real sh(msh)
      parameter (pi=3.1415926535)
      do i=1,mmcat
    1     read(*,*,err=1,end=2) x,y,f,fe,fw,a,b,th,id,iflg
          read(*,*) totsum,dx,dy,bg,tt,tt
          read(*,*) beta,n
          read(*,*) (sh(k),k=3,(n+1)*(n+2)/2+2)
          read(*,*)
          sh(1)=beta
          sh(2)=n
          a00=sh(3)
c     a10=sh(4)
c     a01=sh(5)
          a20=sh(6)
          a11=sh(7)
          a02=sh(8)
          txx=beta**3*sqrt(pi)*(0.5*a00+sqrt(0.5)*a20)
          txy=beta**3*sqrt(pi)*0.5*a11
          tyy=beta**3*sqrt(pi)*(0.5*a00+sqrt(0.5)*a02)
          e1=(txx-tyy)/(txx+tyy)
          e2=2*txy/(txx+tyy)
          write(*,*) id,e1,e2
      enddo

    2 end


      
