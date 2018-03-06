      dimension psfcoef(mmpoly2,msh-2),gkercoef(mmpoly2,msh-2)
      dimension sh(msh),z(64,64),zconv(64,64)

      open(10,file='psf.map',status='old')
      read(10,*) npsf,betapsf
      read(10,*) nfit,xmax,ymax
      read(10,*)
      do k=1,(npsf+1)*(npsf+2)/2
         read(10,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(10,*)
      enddo
      close(10)
      write(0,*) 'Read psf model from psf.map'
      write(0,*) 'Order for shape, spatial variation; beta:',npsf,nfit,betapsf

      open(10,file='gker.map',status='old')
      read(10,*) nker,betaker
      read(10,*) nfitker,xmaxker,ymaxker
      read(10,*)
      do k=1,(nker+1)*(nker+2)/2
         read(10,*) (gkercoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(10,*)
      enddo
      close(10)
      write(0,*) 'Read gaussianization kernal from gker.map'
      write(0,*) 'Order for shape, spatial variation; beta:',nker,nfitker,betaker

c just test the average components (kfit=1)

      call pgbeg(0,'/xw',2,1)

      sh(1)=betapsf
      sh(2)=npsf
      do k=1,(npsf+1)*(npsf+2)/2
         sh(k+2)=psfcoef(1,k)
      enddo 
      call plotshape(sh,64,'PSF')
      call decode(sh,z,64,-64)

      sh(1)=betaker
      sh(2)=nker
      do k=1,(nker+1)*(nker+2)/2
         sh(k+2)=gkercoef(1,k)
      enddo 
      call plotshape(sh,64,'KERNEL')
      call shxpix(sh,z,zconv,64)
      do i=1,64
         do j=1,64
            r2=(i-32.5)**2+(j-32.5)**2
            write(17,*) sqrt(r2),z(i,j),zconv(i,j),exp(-r2/2/betaker**2)/(2*3.14159265*betaker**2)
         enddo
      enddo
      call pgend
      end



      subroutine subbg(pix,mpix,npix,bg,pixout)
c subtract given b/g level from pixelized image.
c result is an array without padding.
      dimension pix(mpix,mpix),pixout(npix,npix)
      do iy=1,npix
          do ix=1,npix
              pixout(ix,iy)=pix(ix,iy)-bg
          enddo
      enddo
      return
      end

      subroutine shxpix(sh,cpix,cout,npix)

c convolve pixel array pix with shapelet object sh
c Use the first npix x npix pixels of the mpix x mpix pixel array pix
c pix can be called as a real npix x xpix array - inside this routine it is
c treated as a complex npix/2 x npix array

      parameter(twopi=2*3.1415926535)
      dimension sh(msh)
      dimension xtab(mm+1),ytab(mm+1)
      complex speq(1024),cpix(npix/2,npix),cout(npix/2,npix),filt


      beta=sh(1)
      n=sh(2)
c      write(0,*) 'beta, n, sum:',beta,n,seriessum(sh)
c Fourier transform pixel array
      do iy=1,npix
         do ix=1,npix/2
            cout(ix,iy)=cpix(ix,iy)
         enddo
      enddo
      call rlft3(cout,speq,npix,npix,1,1)
c      write(0,*) cout(1,1),cout(2,1),cout(1,2),speq(npix/2)
c construct shapelet function at frequency fxnorm, fynorm
c Hn(x/beta)exp(-0.5*x**2/beta**2) --> i^n Hn(f*beta) exp(-0.5*f**2*beta**2)
      fmax2=(sqrt(n*1.)+3)**2
      do iy=1,npix
          fynorm=(mod(iy-1+npix/2,npix)-npix/2)*twopi*beta/npix
          call htab(fynorm,n,ytab)
          do ix=1,npix/2
              fxnorm=(mod(ix-1+npix/2,npix)-npix/2)*twopi*beta/npix
              call htab(fxnorm,n,xtab)
              fnorm2=fxnorm*fxnorm+fynorm*fynorm
              filt=0
              if (fnorm2.gt.fmax2) goto 7 
              fac=exp(-0.5*fnorm2)/sqrt(3.1415926535)*beta*twopi/(npix*npix/2)
              k=0
              do nn=0,n
                  do j=0,nn
                      i=nn-j
                      k=k+1
                      if (j.ne.0) then
                          fac=fac/sqrt(j/real(i+1))
                      else 
                          if (i.ne.0) fac=fac/sqrt(real(i))
                      endif
                      filt=filt+(0,1)**nn*sh(k+2)*fac*xtab(i+1)*ytab(j+1)
                  enddo
                  fac=fac/sqrt(2.)
              enddo
c filter FT's pixel array
    7         cout(ix,iy)=cout(ix,iy)*filt
          enddo
c set the nyquist freq. component to zero. (if it was not zero there is
c a problem anyway!)
          speq(iy)=0
      enddo

c transform back
      call rlft3(cout,speq,npix,npix,1,-1)
      
      return
      end





      SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign)
      INTEGER isign,nn1,nn2,nn3
      COMPLEX data(nn1/2,nn2,nn3),speq(nn2,nn3)
CU    USES fourn
      INTEGER i1,i2,i3,j1,j2,j3,nn(3)
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      COMPLEX c1,c2,h1,h2,w
      c1=cmplx(0.5,0.0)
      c2=cmplx(0.0,-0.5*isign)
      theta=6.28318530717959d0/dble(isign*nn1)
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      nn(1)=nn1/2
      nn(2)=nn2
      nn(3)=nn3
      if(isign.eq.1)then
        call fourn(data,nn,3,isign)
        do 12 i3=1,nn3
          do 11 i2=1,nn2
            speq(i2,i3)=data(1,i2,i3)
11        continue
12      continue
      endif
      do 15 i3=1,nn3
        j3=1
        if (i3.ne.1) j3=nn3-i3+2
        wr=1.0d0
        wi=0.0d0
        do 14 i1=1,nn1/4+1
          j1=nn1/2-i1+2
          do 13 i2=1,nn2
            j2=1
            if (i2.ne.1) j2=nn2-i2+2
            if(i1.eq.1)then
              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
              data(1,i2,i3)=h1+h2
              speq(j2,j3)=conjg(h1-h2)
            else
              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
              data(i1,i2,i3)=h1+w*h2
              data(j1,j2,j3)=conjg(h1-w*h2)
            endif
13        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          w=cmplx(sngl(wr),sngl(wi))
14      continue
15    continue
      if(isign.eq.-1)then
        call fourn(data,nn,3,isign)
      endif
      return
      END

      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      REAL data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
     *k2,n,nprev,nrem,ntot
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END



      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END


      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END

      function fmedian(n,x)
      real x(n),xx(10*mmcat)
      do i=1,n
          xx(i)=x(i)
      enddo
      call sort(n,xx)
      if ((n/2)*2.eq.n) then
          fmedian=( xx(n/2)+xx(n/2+1) ) /2
      else
          fmedian=xx(n/2+1)
      endif
      return
      end


      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

