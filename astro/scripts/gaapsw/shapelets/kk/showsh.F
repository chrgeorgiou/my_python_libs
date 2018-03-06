      real sh(msh),shx(msh),shy(msh),pix(mmpix,mmpix)
      character*11 lbl

      call get2('inimage.fits',nx,ny,pix,mmpix,mmpix)

      call pgbeg(0,'showsh.ps/ps',4,3)
      call pgsch(1.5)
      do i=1,mmcat
          read(*,*,end=99) x,y,f,fe,fw,a,b,th,id,iflg
          read(*,*) sum,dx,dy,bg,fitradius
          read(*,*) beta,n
          read(*,*) (sh(k),k=3,(n+1)*(n+2)/2+2)
          read(*,*)
          sh(1)=beta
          sh(2)=n
          do k=3,(n+1)*(n+2)/2
              sh(k)=sh(k)*sum
          enddo
          call trans1(sh,shx)
          call trans2(sh,shy)
          do k=3,(n+1)*(n+2)/2+2
              sh(k)=sh(k)-dx/beta*shx(k)-dy/beta*shy(k)
          enddo
          write(lbl,'("ID=",i8)') id
          rplot=nint(4*beta)
          call plotshape(sh,2*int(rplot)+1,lbl)
          call compare(pix,mmpix,nx,ny,x,y,bg,rplot,sh)
      enddo
   99 write(0,*) 'Plot written to showps.ps'
      call pgend
      end
