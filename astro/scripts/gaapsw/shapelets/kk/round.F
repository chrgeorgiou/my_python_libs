
      function shround(x,y,n,beta)
      real t(mm/2+1)
      if (mod(n,2).eq.1) then
          shround=0
          return
      endif
      r2=(x*x+y*y)/beta**2
      call ltab(r2,n/2,t)
      shround=exp(-r2/2)*t(n/2+1)/sqrt(3.1415926535)/beta
      end

      subroutine ltab(x,nn,t)
c make table t of Laguerre poly evaluated at x, for order 0..n
      real t(nn+1)
      t(1)=1
      t(2)=1-x
c t(n+1) is poly at order n
      do n=1,nn-1
          t(n+2)= ( (2*n+1-x)*t(n+1)-n*t(n) ) / (n+1)
      enddo
c change convention to make leading-order x term positive
      do n=2,nn,2
          t(n)=-t(n)
      enddo
      return
      end
