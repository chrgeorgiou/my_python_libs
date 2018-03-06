      real dx(mmcat),dy(mmcat),x(mmcat),y(mmcat),flx(mmcat)
      real vx(40,40),vy(40,40),tr(6)
      double precision sum,sumx,sumy

c     read a stellar shapelet catalogue, plot smoothed dx,dy
c     field. Weight by inverse flux.

      do i=1,mmcat
    1     read(*,*,err=1,end=2) x(i),y(i),t,t,t,t,t,t,ii,iflg
          read(*,*) flx(i),dx(i),dy(i),t,t,t
          read(*,*) beta,n
          read(*,*) (t,k=1,(n+1)*(n+2)/2)
          read(*,*)
          if (iflg.ne.0) goto 1
      enddo
    2 nstar=i-1
      write(*,*) 'Read ',nstar,' star models.'

c determine size of the star field
      xmin=x(1)
      xmax=x(1)
      ymin=y(1)
      ymax=y(1)
      do i=1,nstar
          xmin=min(xmin,x(i))
          xmax=max(xmax,x(i))
          ymin=min(ymin,y(i))
          ymax=max(ymax,y(i))
      enddo
c plot smoothed field at 40x40 resolution over full extent. Smooth radius 0.01
      r=min(xmax-xmin,ymax-ymin)/100
      call pgbeg(0,'dxdy.ps/ps',1,1)
      call pgslw(4)
      call pgenv(1.05*xmin-0.05*xmax,1.05*xmax-0.05*xmin,
     :           1.05*ymin-0.05*ymax,1.05*ymax-0.05*ymin,1,0)
      call pglab('X [pix]','Y [pix]','Smoothed dx, dy')
      call pgslw(1)
      call pgpt(nstar,x,y,5)
      efac=0.5/r/r
      do ix=1,40
          xx=xmin+(ix-0.5)*(xmax-xmin)/40.
          do iy=1,40
              yy=ymin+(iy-0.5)*(ymax-ymin)/40.
              sum=0
              sumx=0
              sumy=0
              do i=1,nstar
                  r2=(x(i)-xx)**2+(y(i)-yy)**2
                  w=exp(-r2*efac)*flx(i)
                  sum=sum+w
                  sumx=sumx+dx(i)*w
                  sumy=sumy+dy(i)*w
              enddo
              vx(ix,iy)=sumx/sum
              vy(ix,iy)=sumy/sum
          enddo
      enddo
      d=(xmax-xmin)/40
      tr(1)=xmin-0.5*d
      tr(2)=d
      tr(3)=0
      d=(ymax-ymin)/40
      tr(4)=ymin-0.5*d
      tr(5)=0
      tr(6)=d
      call pgsch(0.1)
      call pgslw(4)
      call pgvect(vx,vy,40,40,1,40,1,40,(xmax-xmin)/40.,0,tr,99999.)
      call pgend
      end

