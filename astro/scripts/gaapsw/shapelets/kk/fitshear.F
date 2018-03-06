
c     Fit a galaxy image as an intrinsically round galaxy, which, after
c     shear and translation by unknown amounts, is convolved with the
c     PSF.

      subroutine fitshearshift_err_012_o4(shape,sherr,g1,g2,g1sig,g2sig,coefs,covar,verbose)
      parameter(mm012=5*(mm/2)-1)
      real shearmat4(mm012,15,0:mm/2),
     :     shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2)
      real shape(*),shape012(mm012),xdata(mm012),sig(mm012),coefs(mm/2+4),
     :    covar(mm/2+4,mm/2+4),alpha(mm/2+4,mm/2+4)
      external funcshift_012_o4
      logical verbose
      common /shearshiftfns_012_o4/ shearmat4,shifted1,shifted2,norder


c      write(0,*) 'ENTERED fitshearshift_err_o4'

      beta=shape(1)
      norder=shape(2)
c fit terms up to order norder-2
      npolar=norder-2
      nsh=5*(norder/2)-6

      call to012(shape,shape012,npolar)
c      write(0,*) 'converted shape to shape012'
c      write(0,*) 'shape:',(i,shape(i),i=1,msh)
c      write(0,*) 'shape012:',(i,shape012(i),i=1,mm012)
c load the data to fit. xdata indexes the shapelet coefs. 
c ydata will be the polar m=0,+/-1,+/-2 shapelets
      do i=1,nsh
          xdata(i)=i
          sig(i)=sherr
c          write(0,*) i,xdata(i),shape012(i),sig(i)
      enddo

c the coefficients to fit are:
c  round shape coefs 0,2,4,...,norder-2
c   g1, g2
c   dx1,dx2
      coefs(1)=1.
      do i=2,norder-2,2
          coefs(i/2+1)=0.
      enddo
      coefs(norder/2+1)=0.
      coefs(norder/2+2)=0.
      coefs(norder/2+3)=0.
      coefs(norder/2+4)=0.
      if (verbose) 
     :    write(0,*) 'FITSHEAR iterations: chi^2, lambda, coefs'
      alam=-1.
      iter=0
      do while (alam.lt.10000. .and. iter.lt.20)
c          write(0,*) 'entering mrqmin:'
          call mrqmin(xdata,shape012,sig,nsh,coefs,norder/2+4,covar,alpha,
     :        mm/2+4,chisq,funcshift_012_o4,alam)
          if (verbose)  
     :        write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+4)
          iter=iter+1
      enddo
      alam=0
      call mrqmin(xdata,shape012,sig,nsh,coefs,norder/2+4,covar,alpha,
     :    mm/2+4,chisq,funcshift_012_o4,alam)
      if (verbose) then
          write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+4)
          write(0,*)
          write(0,*) 'Correlation matrix:'
          do i=1,norder/2+4
              write(0,'(12f7.3)') (covar(i,j)/sqrt(covar(i,i)*covar(j,j)),j=1,norder/2+4)
          enddo
          write(0,*)
          write(0,*) 'Best-fit coefs:'
          write(0,'(12g10.2)') (coefs(i),i=1,norder/2+4)
          write(0,*) 'Sigmas = sqrt (diag):'
          write(0,'(12g10.3)') (sqrt(covar(i,i)),i=1,norder/2+4)
          write(0,*)
      endif
      g1=coefs(norder/2+1)
      g2=coefs(norder/2+2)
      if (covar(norder/2+1,norder/2+1).gt.0. .and. covar(norder/2+2,norder/2+2).gt.0.) then
          g1sig=sqrt(covar(norder/2+1,norder/2+1))
          g2sig=sqrt(covar(norder/2+2,norder/2+2))
      else
          g1sig=-1
          g2sig=-1
      endif
      return
      end






c     Fit a galaxy image as an intrinsically round galaxy, which, after
c     shear and translation by unknown amounts, is convolved with the
c     PSF.

      subroutine fitshearshift_err_o4(shape,sherr,g1,g2,g1sig,g2sig,coefs,covar,verbose)
      
      real shearmat4(msh-2,15,0:mm/2),
     :     shifted1(msh-2,0:mm/2),shifted2(msh-2,0:mm/2)
      real shape(*),xdata(msh-2),sig(msh-2),coefs(mm/2+4),
     :    covar(mm/2+4,mm/2+4),alpha(mm/2+4,mm/2+4)
      external funcshift_o4
      logical verbose
      common /shearshiftfns_o4/ shearmat4,shifted1,shifted2,norder


c      write(0,*) 'ENTERED fitshearshift_err_o4'

      beta=shape(1)
      norder=shape(2)
c fit terms up to order norder-2
      nmin2sh=(norder-1)*(norder)/2

c load the data to fit. xdata indexes the shapelet coefs. 
c ydata will be the polar m=0,+/-1,+/-2 shapelets
      do i=1,nmin2sh
          xdata(i)=i
          sig(i)=sherr
      enddo

c the coefficients to fit are:
c  round shape coefs 0,2,4,...,norder-2
c   g1, g2
c   dx1,dx2
      coefs(1)=1.
      do i=2,norder-2,2
          coefs(i/2+1)=0.
      enddo
      coefs(norder/2+1)=0.
      coefs(norder/2+2)=0.
      coefs(norder/2+3)=0.
      coefs(norder/2+4)=0.

      if (verbose) 
     :    write(0,*) 'FITSHEAR iterations: chi^2, lambda, coefs'
      alam=-1.
      iter=0
      do while (alam.lt.10000. .and. iter.lt.20)
          call mrqmin(xdata,shape(3),sig,nmin2sh,coefs,norder/2+4,covar,alpha,
     :        mm/2+4,chisq,funcshift_o4,alam)
          if (verbose)  
     :        write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+4)
          iter=iter+1
      enddo
      alam=0
      call mrqmin(xdata,shape(3),sig,nmin2sh,coefs,norder/2+4,covar,alpha,
     :    mm/2+4,chisq,funcshift_o4,alam)
      if (verbose) then
          write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+4)
          write(0,*)
          write(0,*) 'Correlation matrix:'
          do i=1,norder/2+4
              write(0,'(12f7.3)') (covar(i,j)/sqrt(covar(i,i)*covar(j,j)),j=1,norder/2+4)
          enddo
          write(0,*)
          write(0,*) 'Best-fit coefs:'
          write(0,'(12g10.2)') (coefs(i),i=1,norder/2+4)
          write(0,*) 'Sigmas = sqrt (diag):'
          write(0,'(12g10.3)') (sqrt(covar(i,i)),i=1,norder/2+4)
          write(0,*)
      endif
      g1=coefs(norder/2+1)
      g2=coefs(norder/2+2)
      if (covar(norder/2+1,norder/2+1).gt.0. .and. covar(norder/2+2,norder/2+2).gt.0.) then
          g1sig=sqrt(covar(norder/2+1,norder/2+1))
          g2sig=sqrt(covar(norder/2+2,norder/2+2))
      else
          g1sig=-1
          g2sig=-1
      endif
      return
      end





c     Fit a galaxy image as an intrinsically round galaxy, which, after
c     shear and translation by unknown amounts, is convolved with the
c     PSF.








      subroutine fitshearshift_err_012_o3(shape,sherr,g1,g2,g1sig,g2sig,coefs,covar,verbose)

      parameter (mm012=5*(mm/2)-1)
      real shape(msh),shape012(mm012)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :     ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2)
     :     ,sheared11(mm012,0:mm/2),sheared12(mm012,0:mm/2),
     :     sheared22(mm012,0:mm/2),
     :     sheared111(mm012,0:mm/2),sheared112(mm012,0:mm/2),
     :     sheared122(mm012,0:mm/2),sheared222(mm012,0:mm/2)

      real xdata(mm012),sig(mm012),coefs(mm/2+4),
     :    covar(mm/2+4,mm/2+4),alpha(mm/2+4,mm/2+4)
      external funcshift012_o3
      logical verbose
      common /shearshiftfns012_o3/ round,sheared1,sheared2,shifted1,shifted2,
     :     sheared11,sheared12,sheared22,
     :     sheared111,sheared112,sheared122,sheared222,norder

      beta=shape(1)
      norder=shape(2)
      write(0,*) 'ENTERED FITSHEARSHIFT_ERR_012_O3'
      nsh=5*(norder/2)-6
      call to012(shape,shape012,norder-2)
c      call plotshape(shape,15,'INPUT SHAPE')

      if (verbose) then
          write(0,*) 'Converted shape to shape012'
          write(0,*) 'm=0:',(shape012(k),k=1,norder/2)
          write(0,*) 'cos:',(shape012(k),k=norder/2+1,2*(norder/2)-1)
          write(0,*) 'sin:',(shape012(k),k=2*(norder/2),3*(norder/2)-2)
          write(0,*) 'cos2:',(shape012(k),k=3*(norder/2)-1,4*(norder/2)-4)
          write(0,*) 'sin2:',(shape012(k),k=4*(norder/2)-3,5*(norder/2)-6)
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Round ',2*i
              write(0,*) 'm=0:',(round(k,i),k=1,norder/2)
              write(0,*) 'cos:',(round(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(round(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(round(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(round(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared1 ',2*i
              write(0,*) 'm=0:',(sheared1(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared1(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared1(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared1(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared1(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared2 ',2*i
              write(0,*) 'm=0:',(sheared2(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared2(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared2(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared2(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared2(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared11 ',2*i
              write(0,*) 'm=0:',(sheared11(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared11(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared11(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared11(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared11(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared12 ',2*i
              write(0,*) 'm=0:',(sheared12(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared12(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared12(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared12(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared12(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared12 ',2*i
              write(0,*) 'm=0:',(sheared22(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared22(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared22(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared22(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared22(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Shifted1 ',2*i
              write(0,*) 'm=0:',(shifted1(k,i),k=1,norder/2)
              write(0,*) 'cos:',(shifted1(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(shifted1(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(shifted1(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(shifted1(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Shifted2 ',2*i
              write(0,*) 'm=0:',(shifted2(k,i),k=1,norder/2)
              write(0,*) 'cos:',(shifted2(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(shifted2(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(shifted2(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(shifted2(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
      endif

c load the data to fit. xdata indexes the shapelet coefs. 
c ydata will be the polar m=0,+/-1,+/-2 shapelets
      do i=1,nsh
          xdata(i)=i
          sig(i)=sherr
      enddo

      if (verbose) then
          do i=1,nsh
              write(0,*) i,xdata(i),shape012(i),sig(i)
          enddo
      endif
c the coefficients to fit are:
c  round shape coefs 0,2,4,...,norder-2
c   g1, g2
c   dx1,dx2
      coefs(1)=1.
      do i=2,norder-2,2
          coefs(i/2+1)=0.
      enddo
      coefs(norder/2+1)=0.
      coefs(norder/2+2)=0.
      coefs(norder/2+3)=0.
      coefs(norder/2+4)=0.

      if (verbose) 
     :    write(0,*) 'FITSHEAR iterations: chi^2, lambda, coefs'
      alam=-1.
      iter=0
      do while (alam.lt.10000. .and. iter.lt.20)
          call mrqmin(xdata,shape012,sig,nsh,coefs,norder/2+4,covar,alpha,
     :        mm/2+4,chisq,funcshift012_o3,alam)
          if (verbose)  
     :        write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+4)
          iter=iter+1
      enddo
      alam=0
      call mrqmin(xdata,shape012,sig,nsh,coefs,norder/2+4,covar,alpha,
     :    mm/2+4,chisq,funcshift012_o3,alam)
      if (verbose) then
          write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+4)
          write(0,*)
          write(0,*) 'Correlation matrix:'
          do i=1,norder/2+4
              write(0,'(12f7.3)') (covar(i,j)/sqrt(covar(i,i)*covar(j,j)),j=1,norder/2+4)
          enddo
          write(0,*)
          write(0,*) 'Best-fit coefs:'
          write(0,'(12g10.2)') (coefs(i),i=1,norder/2+4)
          write(0,*) 'Sigmas = sqrt (diag):'
          write(0,'(12g10.3)') (sqrt(covar(i,i)),i=1,norder/2+4)
          write(0,*)
      endif
      g1=coefs(norder/2+1)
      g2=coefs(norder/2+2)
      if (covar(norder/2+1,norder/2+1).gt.0. .and. covar(norder/2+2,norder/2+2).gt.0.) then
          g1sig=sqrt(covar(norder/2+1,norder/2+1))
          g2sig=sqrt(covar(norder/2+2,norder/2+2))
      else
          g1sig=-1
          g2sig=-1
      endif
      return
      end


c     Fit a galaxy image as an intrinsically round galaxy, which, after
c     shear and translation by unknown amounts, is convolved with the
c     PSF.

      subroutine fitshearshift_err_012_o2(shape,sherr,g1,g2,g1sig,g2sig,coefs,covar,verbose)

      parameter (mm012=5*(mm/2)-1)
      real shape(msh),shape012(mm012)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :     ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2)
     :     ,sheared11(mm012,0:mm/2),sheared12(mm012,0:mm/2),sheared22(mm012,0:mm/2)

      real xdata(mm012),sig(mm012),coefs(mm/2+4),
     :    covar(mm/2+4,mm/2+4),alpha(mm/2+4,mm/2+4)
      external funcshift012_o2
      logical verbose
      common /shearshiftfns012_o2/ round,sheared1,sheared2,shifted1,shifted2,
     :     sheared11,sheared12,sheared22,norder
      

      beta=shape(1)
      if (norder.ne.nint(shape(2))) then
          write(0,*) 'FITSHEARSHIFT_ERR_012_O2: ORDERS DO NOT MATCH'
          stop
      endif
      write(0,*) 'ENTERED FITSHEARSHIFT_ERR_012_O2'
      nsh=5*(norder/2)-6
      call to012(shape,shape012,norder-2)
c      call plotshape(shape,15,'INPUT SHAPE')

      if (verbose) then
          write(0,*) 'Converted shape to shape012'
          write(0,*) 'm=0:',(shape012(k),k=1,norder/2)
          write(0,*) 'cos:',(shape012(k),k=norder/2+1,2*(norder/2)-1)
          write(0,*) 'sin:',(shape012(k),k=2*(norder/2),3*(norder/2)-2)
          write(0,*) 'cos2:',(shape012(k),k=3*(norder/2)-1,4*(norder/2)-4)
          write(0,*) 'sin2:',(shape012(k),k=4*(norder/2)-3,5*(norder/2)-6)
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Round ',2*i
              write(0,*) 'm=0:',(round(k,i),k=1,norder/2)
              write(0,*) 'cos:',(round(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(round(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(round(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(round(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared1 ',2*i
              write(0,*) 'm=0:',(sheared1(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared1(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared1(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared1(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared1(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared2 ',2*i
              write(0,*) 'm=0:',(sheared2(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared2(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared2(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared2(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared2(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared11 ',2*i
              write(0,*) 'm=0:',(sheared11(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared11(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared11(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared11(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared11(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared12 ',2*i
              write(0,*) 'm=0:',(sheared12(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared12(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared12(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared12(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared12(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared12 ',2*i
              write(0,*) 'm=0:',(sheared22(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared22(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared22(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared22(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared22(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Shifted1 ',2*i
              write(0,*) 'm=0:',(shifted1(k,i),k=1,norder/2)
              write(0,*) 'cos:',(shifted1(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(shifted1(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(shifted1(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(shifted1(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Shifted2 ',2*i
              write(0,*) 'm=0:',(shifted2(k,i),k=1,norder/2)
              write(0,*) 'cos:',(shifted2(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(shifted2(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(shifted2(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(shifted2(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
      endif

c load the data to fit. xdata indexes the shapelet coefs. 
c ydata will be the polar m=0,+/-1,+/-2 shapelets
      do i=1,nsh
          xdata(i)=i
          sig(i)=sherr
      enddo

      if (verbose) then
          do i=1,nsh
              write(0,*) i,xdata(i),shape012(i),sig(i)
          enddo
      endif
c the coefficients to fit are:
c  round shape coefs 0,2,4,...,norder-2
c   g1, g2
c   dx1,dx2
      coefs(1)=1.
      do i=2,norder-2,2
          coefs(i/2+1)=0.
      enddo
      coefs(norder/2+1)=0.
      coefs(norder/2+2)=0.
      coefs(norder/2+3)=0.
      coefs(norder/2+4)=0.

      if (verbose) 
     :    write(0,*) 'FITSHEAR iterations: chi^2, lambda, coefs'
      alam=-1.
      iter=0
      do while (alam.lt.10000. .and. iter.lt.20)
          call mrqmin(xdata,shape012,sig,nsh,coefs,norder/2+4,covar,alpha,
     :        mm/2+4,chisq,funcshift012_o2,alam)
          if (verbose)  
     :        write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+4)
          iter=iter+1
      enddo
      alam=0
      call mrqmin(xdata,shape012,sig,nsh,coefs,norder/2+4,covar,alpha,
     :    mm/2+4,chisq,funcshift012_o2,alam)
      if (verbose) then
          write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+4)
          write(0,*)
          write(0,*) 'Correlation matrix:'
          do i=1,norder/2+4
              write(0,'(12f7.3)') (covar(i,j)/sqrt(covar(i,i)*covar(j,j)),j=1,norder/2+4)
          enddo
          write(0,*)
          write(0,*) 'Best-fit coefs:'
          write(0,'(12g10.2)') (coefs(i),i=1,norder/2+4)
          write(0,*) 'Sigmas = sqrt (diag):'
          write(0,'(12g10.3)') (sqrt(covar(i,i)),i=1,norder/2+4)
          write(0,*)
      endif
      g1=coefs(norder/2+1)
      g2=coefs(norder/2+2)
      if (covar(norder/2+1,norder/2+1).gt.0. .and. covar(norder/2+2,norder/2+2).gt.0.) then
          g1sig=sqrt(covar(norder/2+1,norder/2+1))
          g2sig=sqrt(covar(norder/2+2,norder/2+2))
      else
          g1sig=-1
          g2sig=-1
      endif
      return
      end


      subroutine fitshearshift_err_012(shape,sherr,g1,g2,g1sig,g2sig,coefs,covar,verbose)

      parameter (mm012=5*(mm/2)-1)
      real shape(msh),shape012(mm012)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :    ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2)

      real xdata(mm012),sig(mm012),coefs(mm/2+4),
     :    covar(mm/2+4,mm/2+4),alpha(mm/2+4,mm/2+4)
      external funcshift012
      logical verbose
      common /shearshiftfns012/ round,sheared1,sheared2,shifted1,shifted2,norder
      

      beta=shape(1)
      norder=shape(2)

      nsh=5*(norder/2)-6
      call to012(shape,shape012,norder-2)
c      call plotshape(shape,15,'INPUT SHAPE')

      if (verbose) then
          write(0,*) 'Converted shape to shape012'
          write(0,*) 'm=0:',(shape012(k),k=1,norder/2)
          write(0,*) 'cos:',(shape012(k),k=norder/2+1,2*(norder/2)-1)
          write(0,*) 'sin:',(shape012(k),k=2*(norder/2),3*(norder/2)-2)
          write(0,*) 'cos2:',(shape012(k),k=3*(norder/2)-1,4*(norder/2)-4)
          write(0,*) 'sin2:',(shape012(k),k=4*(norder/2)-3,5*(norder/2)-6)
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Round ',2*i
              write(0,*) 'm=0:',(round(k,i),k=1,norder/2)
              write(0,*) 'cos:',(round(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(round(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(round(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(round(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared1 ',2*i
              write(0,*) 'm=0:',(sheared1(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared1(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared1(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared1(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared1(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Sheared2 ',2*i
              write(0,*) 'm=0:',(sheared2(k,i),k=1,norder/2)
              write(0,*) 'cos:',(sheared2(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(sheared2(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(sheared2(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(sheared2(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Shifted1 ',2*i
              write(0,*) 'm=0:',(shifted1(k,i),k=1,norder/2)
              write(0,*) 'cos:',(shifted1(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(shifted1(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(shifted1(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(shifted1(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
          
          do i=0,(norder-2)/2
              write(0,*)
              write(0,*) 'Shifted2 ',2*i
              write(0,*) 'm=0:',(shifted2(k,i),k=1,norder/2)
              write(0,*) 'cos:',(shifted2(k,i),k=norder/2+1,2*(norder/2)-1)
              write(0,*) 'sin:',(shifted2(k,i),k=2*(norder/2),3*(norder/2)-2)
              write(0,*) 'cos2:',(shifted2(k,i),k=3*(norder/2)-1,4*(norder/2)-4)
              write(0,*) 'sin2:',(shifted2(k,i),k=4*(norder/2)-3,5*(norder/2)-6)
          enddo
      endif

c load the data to fit. xdata indexes the shapelet coefs. 
c ydata will be the polar m=0,+/-1,+/-2 shapelets
      do i=1,nsh
          xdata(i)=i
          sig(i)=sherr
      enddo

      if (verbose) then
          do i=1,nsh
              write(0,*) i,xdata(i),shape012(i),sig(i)
          enddo
      endif
c the coefficients to fit are:
c  round shape coefs 0,2,4,...,norder-2
c   g1, g2
c   dx1,dx2
      coefs(1)=1.
      do i=2,norder-2,2
          coefs(i/2+1)=0.
      enddo
      coefs(norder/2+1)=0.
      coefs(norder/2+2)=0.
      coefs(norder/2+3)=0.
      coefs(norder/2+4)=0.

      if (verbose) 
     :    write(0,*) 'FITSHEAR iterations: chi^2, lambda, coefs'
      alam=-1.
      iter=0
      do while (alam.lt.10000. .and. iter.lt.20)
          call mrqmin(xdata,shape012,sig,nsh,coefs,norder/2+4,covar,alpha,
     :        mm/2+4,chisq,funcshift012,alam)
          if (verbose)  
     :        write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+4)
          iter=iter+1
      enddo
      alam=0
      call mrqmin(xdata,shape012,sig,nsh,coefs,norder/2+4,covar,alpha,
     :    mm/2+4,chisq,funcshift012,alam)
      if (verbose) then
          write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+4)
          write(0,*)
          write(0,*) 'Correlation matrix:'
          do i=1,norder/2+4
              write(0,'(12f7.3)') (covar(i,j)/sqrt(covar(i,i)*covar(j,j)),j=1,norder/2+4)
          enddo
          write(0,*)
          write(0,*) 'Best-fit coefs:'
          write(0,'(12g10.2)') (coefs(i),i=1,norder/2+4)
          write(0,*) 'Sigmas = sqrt (diag):'
          write(0,'(12g10.3)') (sqrt(covar(i,i)),i=1,norder/2+4)
          write(0,*)
      endif
      g1=coefs(norder/2+1)
      g2=coefs(norder/2+2)
      if (covar(norder/2+1,norder/2+1).gt.0. .and. covar(norder/2+2,norder/2+2).gt.0.) then
          g1sig=sqrt(covar(norder/2+1,norder/2+1))
          g2sig=sqrt(covar(norder/2+2,norder/2+2))
      else
          g1sig=-1
          g2sig=-1
      endif
      return
      end



      subroutine fitshearshift_err(shape,sherr,g1,g2,g1sig,g2sig,coefs,covar,verbose)

      real shape(msh)
      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
     :    ,shifted1(msh-2,0:mm/2),shifted2(msh-2,0:mm/2)

      real xdata(msh-2),sig(msh-2),coefs(mm/2+4),
     :    covar(mm/2+4,mm/2+4),alpha(mm/2+4,mm/2+4)
      external funcshift
      logical verbose
      common /shearshiftfns/ round,sheared1,sheared2,shifted1,shifted2,norder
      

      beta=shape(1)
      norder=nint(shape(2))-2

c MODIFIED TO FIT ONLY DOWN TO ORDER NORDER
      nsh=(norder+1)*(norder+2)/2

c load the data to fit. xdata indexes the shapelet coefs. 
c ydata will be shape(3..msh)
      do i=1,nsh
          xdata(i)=i
          sig(i)=sherr
      enddo
c the parameters of the fit are:
c  round shape coefs 0,2,4,...,norder-2
c   g1, g2
c   dx1,dx2
      coefs(1)=1.
      do i=2,norder-2,2
          coefs(i/2+1)=0.
      enddo
      coefs(norder/2+1)=0.
      coefs(norder/2+2)=0.
      coefs(norder/2+3)=0.
      coefs(norder/2+4)=0.

      if (verbose) 
     :    write(0,*) 'FITSHEAR iterations: chi^2, lambda, coefs'
      alam=-1.
      iter=0
      do while (alam.lt.10000. .and. iter.lt.20)
          call mrqmin(xdata,shape(3),sig,nsh,coefs,norder/2+4,covar,alpha,
     :        mm/2+4,chisq,funcshift,alam)
          if (verbose)  
     :        write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+4)
          iter=iter+1
      enddo
      alam=0
      call mrqmin(xdata,shape(3),sig,nsh,coefs,norder/2+4,covar,alpha,
     :    mm/2+4,chisq,funcshift,alam)
      if (verbose) then
          write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+4)
          write(0,*)
          write(0,*) 'Correlation matrix:'
          do i=1,norder/2+4
              write(0,'(12f7.3)') (covar(i,j)/sqrt(covar(i,i)*covar(j,j)),j=1,norder/2+4)
          enddo
          write(0,*)
          write(0,*) 'Best-fit coefs:'
          write(0,'(12g10.2)') (coefs(i),i=1,norder/2+4)
          write(0,*) 'Sigmas = sqrt (diag):'
          write(0,'(12g10.3)') (sqrt(covar(i,i)),i=1,norder/2+4)
          write(0,*)
      endif
      g1=coefs(norder/2+1)
      g2=coefs(norder/2+2)
      if (covar(norder/2+1,norder/2+1).gt.0. .and. covar(norder/2+2,norder/2+2).gt.0.) then
          g1sig=sqrt(covar(norder/2+1,norder/2+1))
          g2sig=sqrt(covar(norder/2+2,norder/2+2))
      else
          g1sig=-1
          g2sig=-1
      endif
      return
      end


      subroutine fitshear_err(shape,sherr,g1,g2,g1sig,g2sig,coefs,covar,verbose)

      real shape(msh)
      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)

      real xdata(msh-2),sig(msh-2),coefs(mm/2+2),
     :    covar(mm/2+2,mm/2+2),alpha(mm/2+2,mm/2+2)
      external funcs
      logical verbose
      common /shearfns/ round,sheared1,sheared2,norder
      

      beta=shape(1)
      norder=shape(2)
      nsh=(norder+1)*(norder+2)/2

c load the data to fit. xdata indexes the shapelet coefs. 
c ydata will be shape(3..msh)
      do i=1,nsh
          xdata(i)=i
          sig(i)=sherr
      enddo

      coefs(1)=1.
      do i=1,norder-2,2
          coefs(i/2+1)=0.1
      enddo
      coefs(norder/2+1)=0.
      coefs(norder/2+2)=0.

      if (verbose) 
     :    write(0,*) 'FITSHEAR iterations: chi^2, lambda, coefs'
      alam=-1.
      iter=0
      do while (alam.lt.10000. .and. iter.lt.20)
          call mrqmin(xdata,shape(3),sig,nsh,coefs,norder/2+2,covar,alpha,
     :        mm/2+2,chisq,funcs,alam)
          if (verbose)  
     :        write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+2)
          iter=iter+1
      enddo
      alam=0
      call mrqmin(xdata,shape(3),sig,nsh,coefs,norder/2+2,covar,alpha,
     :    mm/2+2,chisq,funcs,alam)
      if (verbose) then
          write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+2)
          write(0,*)
          write(0,*) 'Correlation matrix:'
          do i=1,norder/2+2
              write(0,'(12f7.3)') (covar(i,j)/sqrt(covar(i,i)*covar(j,j)),j=1,norder/2+2)
          enddo
          write(0,*)
          write(0,*) 'Best-fit coefs:'
          write(0,'(12g10.2)') (coefs(i),i=1,norder/2+2)
          write(0,*) 'Sigmas = sqrt (diag):'
          write(0,'(12g10.3)') (sqrt(covar(i,i)),i=1,norder/2+2)
          write(0,*)
      endif
      g1=coefs(norder/2+1)
      g2=coefs(norder/2+2)
      g1sig=sqrt(covar(norder/2+1,norder/2+1))
      g2sig=sqrt(covar(norder/2+2,norder/2+2))
      return
      end


      subroutine fitshear(shape,g1,g2,g1sig,g2sig,coefs,covar,verbose)

      real shape(msh)
      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)

      real xdata(msh-2),sig(msh-2),coefs(mm/2+2),
     :    covar(mm/2+2,mm/2+2),alpha(mm/2+2,mm/2+2)
      external funcs
      logical verbose
      common /shearfns/ round,sheared1,sheared2,norder
      

      beta=shape(1)
      norder=shape(2)
      nsh=(norder+1)*(norder+2)/2

c load the data to fit. xdata indexes the shapelet coefs. 
c ydata will be shape(3..msh)
      do i=1,nsh
          xdata(i)=i
          sig(i)=1.
      enddo

      coefs(1)=1.
      do i=1,norder-2,2
          coefs(i/2+1)=0.1
      enddo
      coefs(norder/2+1)=0.1
      coefs(norder/2+2)=0.1

      if (verbose) 
     :    write(0,*) 'FITSHEAR iterations: chi^2, lambda, coefs'
      alam=-1.
      iter=0
      do while (alam.lt.10000. .and. iter.lt.20)
          call mrqmin(xdata,shape(3),sig,nsh,coefs,norder/2+2,covar,alpha,
     :        mm/2+2,chisq,funcs,alam)
          if (verbose)  
     :        write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+2)
          iter=iter+1
      enddo
      alam=0
      call mrqmin(xdata,shape(3),sig,nsh,coefs,norder/2+2,covar,alpha,
     :    mm/2+2,chisq,funcs,alam)
      if (verbose) then
          write(0,'(20g11.3)') chisq,alam,(coefs(i),i=1,norder/2+2)
          write(0,*)
          write(0,*) 'Correlation matrix:'
          do i=1,norder/2+2
              write(0,'(12f7.3)') (covar(i,j)/sqrt(covar(i,i)*covar(j,j)),j=1,norder/2+2)
          enddo
          write(0,*)
          write(0,*) 'Best-fit coefs:'
          write(0,'(12g10.2)') (coefs(i),i=1,norder/2+2)
          write(0,*) 'Sigmas = sqrt (diag):'
          write(0,'(12g10.3)') (sqrt(covar(i,i)),i=1,norder/2+2)
          write(0,*)
      endif
      g1=coefs(norder/2+1)
      g2=coefs(norder/2+2)
      g1sig=sqrt(covar(norder/2+1,norder/2+1))
      g2sig=sqrt(covar(norder/2+2,norder/2+2))
      return
      end




c     The function to fit: Shapelet term `cpt' is built from the
c     corresponding terms in round, sheared1, sheared2, shifted1 and shifted2

      subroutine funcshift012_o3(cpt,a,y,dyda,ma)

      parameter (mm012=5*mm/2-1)
      real a(mm/2+4),dyda(mm/2+4)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :     ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2)
     :     ,sheared11(mm012,0:mm/2),sheared12(mm012,0:mm/2),
     :     sheared22(mm012,0:mm/2),
     :     sheared111(mm012,0:mm/2),sheared112(mm012,0:mm/2),
     :     sheared122(mm012,0:mm/2),sheared222(mm012,0:mm/2)
      double precision sum,sumg1,sumg2,sumt1,sumt2
      common /shearshiftfns012_o3/ round,sheared1,sheared2,shifted1,shifted2,
     :     sheared11,sheared12,sheared22,
     :     sheared111,sheared112,sheared122,sheared222,n

      if (ma.ne.n/2+4) then
          write(0,*) '****  FUNCSHIFT012_O3: # VARIABLES MISMATCHES ****'
          stop
      endif
      icpt=nint(cpt)
      g1=a(n/2+1)
      g2=a(n/2+2)
      dx=a(n/2+3)
      dy=a(n/2+4)
      sum=0
      sumg1=0
      sumg2=0
      sumt1=0
      sumt2=0
c do not include the highest order of cpt since this makes for degeneracies
      do i=0,n-2,2
         i2=i/2
         ci=a(i2+1)
         r=round(icpt,i2)
         t1=shifted1(icpt,i2)
         t2=shifted2(icpt,i2)
         s1=sheared1(icpt,i2)
         s2=sheared2(icpt,i2)
         s11=sheared11(icpt,i2)
         s12=sheared12(icpt,i2)
         s22=sheared22(icpt,i2)
         s111=sheared111(icpt,i2)
         s112=sheared112(icpt,i2)
         s122=sheared122(icpt,i2)
         s222=sheared222(icpt,i2)
         sumg1=sumg1+ci*(s1+g1*s11+g2*s12+0.5*g1*g1*s111+g1*g2*s112+0.5*g2*g2*s122)
         sumg2=sumg2+ci*(s2+g1*s12+g2*s22+0.5*g1*g1*s112+g1*g2*s122+0.5*g2*g2*s222)
         sumt1=sumt1+ci*t1
         sumt2=sumt2+ci*t2
         dyda(i2+1)=r+dx*t1+dy*t2+g1*s1+g2*s2+g1*g1/2*s11+g1*g2*s12+g2*g2/2*s22+g1*g1*g1/6*s111+g1*g1*g2/2*s112+g1*g2*g2/2*s122+g2*g2*g2/6*s222
         sum=sum+ci*dyda(i2+1)
      enddo
      y=sum
      dyda(n/2+1)=sumg1
      dyda(n/2+2)=sumg2
      dyda(n/2+3)=sumt1
      dyda(n/2+4)=sumt2

c      write(0,*) 'FUNCS:',cpt,(a(i),i=1,n/2+4)
c      write(0,*) y,(dyda(i),i=1,n/2+4)
c      write(0,*) (round(icpt,i),i=0,n/2)
c      write(0,*) (sheared1(icpt,i),i=0,n/2)
c      write(0,*) (sheared2(icpt,i),i=0,n/2)
c      write(0,*)
      return
      end



c     The function to fit: Shapelet term `cpt' is built from the
c     corresponding terms in round, sheared1, sheared2, shifted1 and shifted2

      subroutine funcshift012_o2(cpt,a,y,dyda,ma)

      parameter (mm012=5*mm/2-1)
      real a(mm/2+4),dyda(mm/2+4)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :     ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2)
     :     ,sheared11(mm012,0:mm/2),sheared12(mm012,0:mm/2),sheared22(mm012,0:mm/2)
      double precision sum,sumg1,sumg2,sumt1,sumt2
      common /shearshiftfns012_o2/ round,sheared1,sheared2,shifted1,shifted2,
     :     sheared11,sheared12,sheared22,n

      if (ma.ne.n/2+4) then
          write(0,*) '****  FUNCSHIFT012_O2: # VARIABLES MISMATCHES ****'
          stop
      endif
      icpt=nint(cpt)
      g1=a(n/2+1)
      g2=a(n/2+2)
      dx=a(n/2+3)
      dy=a(n/2+4)
      sum=0
      sumg1=0
      sumg2=0
      sumt1=0
      sumt2=0
c do not include the highest order of cpt since this makes for degeneracies
      do i=0,n-2,2
          sumg1=sumg1+a(i/2+1)*(sheared1(icpt,i/2)+g1*sheared11(icpt,i/2)+g2*sheared12(icpt,i/2))
          sumg2=sumg2+a(i/2+1)*(sheared2(icpt,i/2)+g1*sheared12(icpt,i/2)+g2*sheared22(icpt,i/2))
          sumt1=sumt1+a(i/2+1)*shifted1(icpt,i/2)
          sumt2=sumt2+a(i/2+1)*shifted2(icpt,i/2)
          dyda(i/2+1)=round(icpt,i/2)+g1*sheared1(icpt,i/2)
     :         +g2*sheared2(icpt,i/2)+dx*shifted1(icpt,i/2)+dy*shifted2(icpt,i/2)
     :         +g1*g1/2*sheared11(icpt,i/2)+g1*g2*sheared12(icpt,i/2)+g2*g2/2*sheared22(icpt,i/2)
          sum=sum+a(i/2+1)*dyda(i/2+1)
      enddo
      y=sum
      dyda(n/2+1)=sumg1
      dyda(n/2+2)=sumg2
      dyda(n/2+3)=sumt1
      dyda(n/2+4)=sumt2

c      write(0,*) 'FUNCS: Cpt',icpt
c      write(0,*) y,(dyda(i),i=1,n/2+4)
c      write(0,*) (round(icpt,i),i=0,n/2-1)
c      write(0,*) (sheared1(icpt,i),i=0,n/2-1)
c      write(0,*) (sheared2(icpt,i),i=0,n/2-1)
c      write(0,*) (sheared11(icpt,i),i=0,n/2-1)
c      write(0,*) (sheared12(icpt,i),i=0,n/2-1)
c      write(0,*) (sheared22(icpt,i),i=0,n/2-1)
c      write(0,*)
      return
      end


      subroutine funcshift012(cpt,a,y,dyda,ma)

      parameter (mm012=5*mm/2-1)
      real a(mm/2+4),dyda(mm/2+4)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :    ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2)
      double precision sum,sumg1,sumg2,sumt1,sumt2
      common /shearshiftfns012/ round,sheared1,sheared2,shifted1,shifted2,n

      if (ma.ne.n/2+4) then
          write(0,*) '****  FUNCSHIFT012: # VARIABLES MISMATCHES ****'
          stop
      endif
      icpt=nint(cpt)
      g1=a(n/2+1)
      g2=a(n/2+2)
      dx=a(n/2+3)
      dy=a(n/2+4)
      sum=0
      sumg1=0
      sumg2=0
      sumt1=0
      sumt2=0
c do not include the highest order of cpt since this makes for degeneracies
      do i=0,n-2,2
          sumg1=sumg1+a(i/2+1)*sheared1(icpt,i/2)
          sumg2=sumg2+a(i/2+1)*sheared2(icpt,i/2)
          sumt1=sumt1+a(i/2+1)*shifted1(icpt,i/2)
          sumt2=sumt2+a(i/2+1)*shifted2(icpt,i/2)
          dyda(i/2+1)=round(icpt,i/2)+g1*sheared1(icpt,i/2)
     :        +g2*sheared2(icpt,i/2)+dx*shifted1(icpt,i/2)+dy*shifted2(icpt,i/2)
          sum=sum+a(i/2+1)*dyda(i/2+1)
      enddo
      y=sum
      dyda(n/2+1)=sumg1
      dyda(n/2+2)=sumg2
      dyda(n/2+3)=sumt1
      dyda(n/2+4)=sumt2

c      write(0,*) 'FUNCS:',cpt,(a(i),i=1,n/2+4)
c      write(0,*) y,(dyda(i),i=1,n/2+4)
c      write(0,*) (round(icpt,i),i=0,n/2)
c      write(0,*) (sheared1(icpt,i),i=0,n/2)
c      write(0,*) (sheared2(icpt,i),i=0,n/2)
c      write(0,*)
      return
      end



      subroutine funcshift(cpt,a,y,dyda,ma)

      real a(mm/2+4),dyda(mm/2+4)
      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
     :    ,shifted1(msh-2,0:mm/2),shifted2(msh-2,0:mm/2)
      double precision sum,sumg1,sumg2,sumt1,sumt2
      common /shearshiftfns/ round,sheared1,sheared2,shifted1,shifted2,n

      if (ma.ne.n/2+4) then
          write(0,*) '****  FUNCSHIFT: # VARIABLES MISMATCHES ****'
          stop
      endif
      icpt=nint(cpt)
      g1=a(n/2+1)
      g2=a(n/2+2)
      dx=a(n/2+3)
      dy=a(n/2+4)
      sum=0
      sumg1=0
      sumg2=0
      sumt1=0
      sumt2=0
c do not include the highest order of cpt since this makes for degeneracies
      do i=0,n-2,2
          sumg1=sumg1+a(i/2+1)*sheared1(icpt,i/2)
          sumg2=sumg2+a(i/2+1)*sheared2(icpt,i/2)
          sumt1=sumt1+a(i/2+1)*shifted1(icpt,i/2)
          sumt2=sumt2+a(i/2+1)*shifted2(icpt,i/2)
          dyda(i/2+1)=round(icpt,i/2)+g1*sheared1(icpt,i/2)
     :        +g2*sheared2(icpt,i/2)+dx*shifted1(icpt,i/2)+dy*shifted2(icpt,i/2)
          sum=sum+a(i/2+1)*dyda(i/2+1)
      enddo
      y=sum
      dyda(n/2+1)=sumg1
      dyda(n/2+2)=sumg2
      dyda(n/2+3)=sumt1
      dyda(n/2+4)=sumt2

c      write(0,*) 'FUNCS:',cpt,(a(i),i=1,n/2+4)
c      write(0,*) y,(dyda(i),i=1,n/2+4)
c      write(0,*) (round(icpt,i),i=0,n/2)
c      write(0,*) (sheared1(icpt,i),i=0,n/2)
c      write(0,*) (sheared2(icpt,i),i=0,n/2)
c      write(0,*)
      return
      end



      subroutine funcshift_012_o4(cpt,a,y,dyda,ma)

      parameter(mm012=5*(mm/2)-1)
      real a(mm/2+4),dyda(mm/2+4)
      real shearmat4(mm012,15,0:mm/2),
     :     shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2)
      double precision sum,sumg1,sumg2,sumt1,sumt2
      common /shearshiftfns_012_o4/ shearmat4,shifted1,shifted2,n

c      write(0,*) 'ENTERED funcshift_o4',cpt,(a(kk),kk=1,ma)

      if (ma.ne.n/2+4) then
          write(0,*) '****  FUNCSHIFT: # VARIABLES MISMATCHES ****', ma,n
          stop
      endif
      icpt=nint(cpt)
      g1=a(n/2+1)
      g2=a(n/2+2)
      dx=a(n/2+3)
      dy=a(n/2+4)
      sum=0
      sumg1=0
      sumg2=0
      sumt1=0
      sumt2=0
c do not include the highest order of cpt since this makes for degeneracies
      do i=0,n-2,2
          coefi=a(i/2+1)
          dyda(i/2+1)=0
          sumt1=sumt1+coefi*shifted1(icpt,i/2)
          sumt2=sumt2+coefi*shifted2(icpt,i/2)
          k=0
c*************** loop up to desired order (max 4)
          do ig=0,4
              do ig2=0,ig
                  ig1=ig-ig2
                  k=k+1
                  if (ig1.gt.0) sumg1=sumg1+coefi*ig1*g1**(ig1-1)*g2**ig2*shearmat4(icpt,k,i/2)
                  if (ig2.gt.0) sumg2=sumg2+coefi*g1**ig1*ig2*g2**(ig2-1)*shearmat4(icpt,k,i/2)
                  dyda(i/2+1)=dyda(i/2+1)+g1**ig1*g2**ig2*shearmat4(icpt,k,i/2)
              enddo
          enddo
          dyda(i/2+1)=dyda(i/2+1)+dx*shifted1(icpt,i/2)+dy*shifted2(icpt,i/2)
          sum=sum+coefi*dyda(i/2+1)
      enddo
      y=sum
      dyda(n/2+1)=sumg1
      dyda(n/2+2)=sumg2
      dyda(n/2+3)=sumt1
      dyda(n/2+4)=sumt2

c      write(0,*) 'FUNCS:',cpt,(a(i),i=1,n/2+4)
c      write(0,*) y,(dyda(i),i=1,n/2+4)
c      write(0,*) (round(icpt,i),i=0,n/2)
c      write(0,*) (sheared1(icpt,i),i=0,n/2)
c      write(0,*) (sheared2(icpt,i),i=0,n/2)
c      write(0,*)
      return
      end


      subroutine funcshift_o4(cpt,a,y,dyda,ma)

      real a(mm/2+4),dyda(mm/2+4)
      real shearmat4(msh-2,15,0:mm/2),
     :     shifted1(msh-2,0:mm/2),shifted2(msh-2,0:mm/2)
      double precision sum,sumg1,sumg2,sumt1,sumt2
      common /shearshiftfns_o4/ shearmat4,shifted1,shifted2,n

c      write(0,*) 'ENTERED funcshift_o4',cpt,(a(kk),kk=1,ma)

      if (ma.ne.n/2+4) then
          write(0,*) '****  FUNCSHIFT: # VARIABLES MISMATCHES ****'
          stop
      endif
      icpt=nint(cpt)
      g1=a(n/2+1)
      g2=a(n/2+2)
      dx=a(n/2+3)
      dy=a(n/2+4)
      sum=0
      sumg1=0
      sumg2=0
      sumt1=0
      sumt2=0
c do not include the highest order of cpt since this makes for degeneracies
      do i=0,n-2,2
          coefi=a(i/2+1)
          dyda(i/2+1)=0
          sumt1=sumt1+coefi*shifted1(icpt,i/2)
          sumt2=sumt2+coefi*shifted2(icpt,i/2)
          k=0
c*************** loop up to desired order (max 4)
          do ig=0,4
              do ig2=0,ig
                  ig1=ig-ig2
                  k=k+1
                  if (ig1.gt.0) sumg1=sumg1+coefi*ig1*g1**(ig1-1)*g2**ig2*shearmat4(icpt,k,i/2)
                  if (ig2.gt.0) sumg2=sumg2+coefi*g1**ig1*ig2*g2**(ig2-1)*shearmat4(icpt,k,i/2)
                  dyda(i/2+1)=dyda(i/2+1)+g1**ig1*g2**ig2*shearmat4(icpt,k,i/2)  
              enddo
          enddo
          dyda(i/2+1)=dyda(i/2+1)+dx*shifted1(icpt,i/2)+dy*shifted2(icpt,i/2)
          sum=sum+coefi*dyda(i/2+1)
      enddo
      y=sum
      dyda(n/2+1)=sumg1
      dyda(n/2+2)=sumg2
      dyda(n/2+3)=sumt1
      dyda(n/2+4)=sumt2

c      write(0,*) 'FUNCS:',cpt,(a(i),i=1,n/2+4)
c      write(0,*) y,(dyda(i),i=1,n/2+4)
c      write(0,*) (round(icpt,i),i=0,n/2)
c      write(0,*) (sheared1(icpt,i),i=0,n/2)
c      write(0,*) (sheared2(icpt,i),i=0,n/2)
c      write(0,*)
      return
      end




c     The function to fit: Shapelet term `cpt' is built from the
c     corresponding terms in round, sheared1 and sheared2

      subroutine funcs(cpt,a,y,dyda,ma)

      real a(mm/2+2),dyda(mm/2+2)
      real round(msh-2,0:mm/2),sheared1(msh-2,0:mm/2),sheared2(msh-2,0:mm/2)
      double precision sum,sumg1,sumg2
      common /shearfns/ round,sheared1,sheared2,n

      icpt=nint(cpt)
      g1=a(n/2+1)
      g2=a(n/2+2)
      sum=0
      sumg1=0
      sumg2=0
c do not include the highest order of cpt since this makes for degeneracies
      do i=0,n-2,2
          sumg1=sumg1+a(i/2+1)*sheared1(icpt,i/2)
          sumg2=sumg2+a(i/2+1)*sheared2(icpt,i/2)
          dyda(i/2+1)=round(icpt,i/2)+g1*sheared1(icpt,i/2)+g2*sheared2(icpt,i/2)
          sum=sum+a(i/2+1)*dyda(i/2+1)
      enddo
      y=sum
      dyda(n/2+1)=sumg1
      dyda(n/2+2)=sumg2
c      write(0,*) 'FUNCS:',cpt,(a(i),i=1,n/2+3)
c      write(0,*) y,(dyda(i),i=1,n/2+3)
c      write(0,*) (round(icpt,i),i=0,n/2)
c      write(0,*) (sheared1(icpt,i),i=0,n/2)
c      write(0,*) (sheared2(icpt,i),i=0,n/2)
c      write(0,*)
      return
      end



