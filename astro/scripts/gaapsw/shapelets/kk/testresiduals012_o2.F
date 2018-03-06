c define shape, make shapelet expansion, make plot of fit vs data.

      parameter (mm012=5*(mm/2)-1,msize=99,mbox=2*msize+1)
      real z(mbox,mbox),zerr(mbox,mbox),xplot(mbox)
      real shape(msh),sh012(mm012)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :    ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2),coefs(mm/2+4),
     :     sheared11(mm012,0:mm/2),sheared12(mm012,0:mm/2),
     :    sheared22(mm012,0:mm/2),
     :    covar(mm/2+4,mm/2+4)

      common /shearshiftfns012_o2/ round,sheared1,sheared2,shifted1,shifted2,sheared11,sheared12,sheared22,nshear
c      common /shearfns/ round,sheared1,sheared2,nshear

      write(*,*) 'n,nsersic?'
      read(*,*) n,nsersic
      nsh=(n+1)*(n+2)/2
c      write(*,*) 'Input shears g1,g2?'
c      read(*,*) g1in,g2in

      do i=1,mbox
          do j=1,mbox
              zerr(j,i)=1e-3
          enddo
      enddo

      call pgbeg(0,'/xs',4,4)
      call pgask(.false.)
      do i=1,mbox
          xplot(i)=i-(msize+1.)
      enddo
      nshear=n
      rpct=20.
      do ig1=0,14
          g2in=0.05*ig1
          g1in=0
          do i=1,mbox
              do j=1,mbox
                  z(i,j)=0
                  do ix=0,0
                      x=(i-(msize+1.))+0.05*ix
                      do iy=0,0
                          y=(j-(msize+1.))+0.05*iy
                          xx=x*(1-g1in)-y*g2in
                          yy=-x*g2in+y*(1+g1in)
                          r=sqrt(xx*xx+yy*yy)
                          z(i,j)=z(i,j)+exp(-log(100.)*(r/rpct)**(1./nsersic))
                      enddo
                  enddo
              enddo
          enddo

          do ibeta=15,120,5
              beta=ibeta*0.1
              bg=0.
              call fitshapeletserr(z,zerr,mbox,mbox,mbox,(msize+1.),(msize+1.),4*beta,bg,shape,sherr,n,beta)
              sum=seriessum(shape)
              do i=3,nsh+2
                  shape(i)=shape(i)/sum
              enddo
              call plotshape(shape,mbox/4,'SHAPELET FIT')
c              call to012(shape,sh012,n)
c              write(10,'(40g13.5)') beta,(sh012(k),k=1,5*(n/2)-1)
              call setroundshearshift_012_o2_nopsf(beta,n,round,sheared1,sheared2,shifted1,shifted2,sheared11,sheared12,sheared22)

              call testfunc

              call fitshearshift_err_012_o2(shape,sherr,g1,g2,g1sig,g2sig,coefs,covar,.false.)
              write(*,*) n,nsersic,beta,bg,g1in,g1,(coefs(i),i=1,n/2),g2in,g2
          enddo
      enddo

      call pgend
      end



      subroutine testfunc
      parameter (mm012=5*(mm/2)-1,msize=99,mbox=2*msize+1)
      real c(mm/2+4),cgrad(mm/2+4),cgrad2(mm/2+4)
      real round(mm012,0:mm/2),sheared1(mm012,0:mm/2),sheared2(mm012,0:mm/2)
     :    ,shifted1(mm012,0:mm/2),shifted2(mm012,0:mm/2),coefs(mm/2+4),
     :     sheared11(mm012,0:mm/2),sheared12(mm012,0:mm/2),
     :    sheared22(mm012,0:mm/2)
      common /shearshiftfns012_o2/ round,sheared1,sheared2,shifted1,shifted2,sheared11,sheared12,sheared22,nshear

      np=nshear/2+4
      n=nshear
      do ii=1,3
          do k=1,np
              c(k)=5*ran3(idum)-0.5
          enddo
          do icpt=1,5*(nshear/2)-6
              write(0,*) 'FUNCS: Cpt',icpt
c     write(0,*) y,(dyda(i),i=1,n/2+4)
              write(0,*) (round(icpt,i),i=0,n/2-1)
              write(0,*) (sheared1(icpt,i),i=0,n/2-1)
              write(0,*) (sheared2(icpt,i),i=0,n/2-1)
              write(0,*) (sheared11(icpt,i),i=0,n/2-1)
              write(0,*) (sheared12(icpt,i),i=0,n/2-1)
              write(0,*) (sheared22(icpt,i),i=0,n/2-1)
              write(0,*)
              call funcshift012_o2(real(icpt),c,y,cgrad,np)
              do k=1,np
                  c(k)=c(k)+0.01
                  call funcshift012_o2(real(icpt),c,y2,cgrad2,np)
                  c(k)=c(k)-0.01
                  write(0,*) icpt,k,c(k),(y2-y)/.01,cgrad(k)
              enddo
          enddo
          write(0,*)
      enddo
      return
      end
