      real legendre
      do i=-100,100
         x=0.01*i
         write(*,'(13f10.5)') x,(cheby(x,n),n=0,5),(legendre(x,n),n=0,5)
      enddo
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
