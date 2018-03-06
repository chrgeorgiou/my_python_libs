      parameter (mmo4=8,msho4=(mmo4+1)*(mmo4+2)/2+2)
      parameter (norder=4,norder2=(norder+1)*(norder+2)/2)

      real vec(norder2)
      real shearedcirc4(msho4-2,0:mmo4/2,norder2)
      
      do k=1,norder2
          do ic=0,mmo4
              do ish=1,msho4-2
                  shearedcirc4(ish,ic/2,k)=0
              enddo
          enddo
      enddo

      do kount=1,99999
    1     read (*,*,err=1,end=2) ic,k,l,m,(vec(i),i=1,nterm(l,m,ic))
          ish=k*(k+1)/2+m+1
          i=1
          n0=abs(k-ic)/2
          do k2=n0,norder,2
              do m2=mod(m,2),k2,2
                  l2=k2-m2
c icoef is the index of the 00,10,01,20,11,02, etc terms g1^l2 * g2^m2
                  icoef=k2*(k2+1)/2+m2+1
                  shearedcirc4(ish,ic/2,icoef)=vec(i)
                  i=i+1
              enddo
          enddo
      enddo
      
    2 write(*,1000)
 1000 format ('      data shearedcirc4 /')
      do ic=0,mmo4,2
          k=0
          do ntot=0,norder
              do ng2=0,ntot
                  ng1=ntot-ng2
                  k=k+1
                  write(*,2000) ic,ng1,ng2
 2000             format ('C       circ order ',i4,'   order g1,g2 ',2i3)
                  write(*,1001) (shearedcirc4(i,ic/2,k),i=1,msho4-2)
 1001             format (('     : ',6(f10.6,', ')))
              enddo
          enddo
      enddo
      write(*,1002)
 1002 format ('     :  /')
      end





      function nterm(l,m,ic)
      parameter (norder=4,norder2=(norder+1)*(norder+2)/2)

c     number of terms up to total order 4 in g1,g2 that are generated at
c     shapelet order l,m by basis function ic when sheared.

c k=total order of shapelet; n0=leading total order of terms generated
      k=l+m
      n0=abs(k-ic)/2
      nterm=0
      do n=n0,norder,2
          do m2=mod(m,2),n,2
              nterm=nterm+1
          enddo
      enddo
      return
      end
