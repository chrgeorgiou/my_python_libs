
      real x(10),y(10),sig(10),coef(2),cov(10,10),alpha(10,10)
      external funcs

      a=2
      b=3
      do i=1,10
         x(i)=i
         y(i)=a+b*x(i)+a*b/x(i)
         sig(i)=1
      enddo
      alam=-1
      do k=1,20
         call mrqmin(x,y,sig,10,coef,2,cov,alpha,10,chisq,funcs,alam)
         write(0,*) alam,chisq,coef
      enddo
      alam=0
      call mrqmin(x,y,sig,10,coef,2,cov,alpha,10,chisq,funcs,alam)
      write(0,*) alam,chisq,coef

      end




      subroutine funcs(x,coef,y,dyda,ma)
      real coef(ma),dyda(ma)
      a=coef(1)
      b=coef(2)
      y=a+b*x+a*b/x
      dyda(1)=1+b/x
      dyda(2)=x+a/x
      return
      end
