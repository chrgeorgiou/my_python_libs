      parameter (mmplot=10000000)
      real xplot(mmplot),yplot(mmplot)
      character*20 lab1,lab2,lab3*50,filename*80
      character ch
      logical diag

      iaxes=0
      just=0

      read(*,*) jcol1
      read(*,'(a)') lab1
      read(*,*) jcol2
      read(*,'(a)') lab2
      read(*,'(a)') lab3
      icol1=abs(jcol1)
      icol2=abs(jcol2)
      if (jcol1.lt.0) iaxes=iaxes+10
      if (jcol2.lt.0) iaxes=iaxes+20
      xmax=-1e30
      ymax=xmax
      xmin=-xmax
      ymin=xmin
      do i=1,mmplot
   99     if (icol1.lt.icol2) then
              read(*,*,end=1,err=99) (t,k=1,icol1-1),xplot(i),
     :                     (t,k=icol1+1,icol2-1),yplot(i)
          elseif (icol2.lt.icol1) then
              read(*,*,end=1,err=99) (t,k=1,icol2-1),yplot(i),
     :                     (t,k=icol2+1,icol1-1),xplot(i)
          else
              write(0,*) 'You asked for the same column twice.'
              stop
          endif
          if (jcol1.lt.0) xplot(i)=log10(xplot(i))
          if (jcol2.lt.0) yplot(i)=log10(yplot(i))
          xmin=min(xmin,xplot(i))
          xmax=max(xmax,xplot(i))
          ymin=min(ymin,yplot(i))
          ymax=max(ymax,yplot(i))
      enddo
      write(0,*) '***Warning: max no of input lines reached: ',mmplot

    1 if (xmin*xmax.le.0. .and. ymin*ymax.le.0.) iaxes=iaxes+1
      if (((xmax-xmin).gt.0.5*(ymax-ymin)) 
     :    .and.((xmax-xmin).lt.2*(ymax-ymin))) just=1

      nplot=i-1
      dx=(xmax-xmin)*0.05
      if (dx.eq.0.) dx=0.1
      xmax=xmax+dx
      xmin=xmin-dx
      dy=(ymax-ymin)*0.05
      if (dy.eq.0.) dy=0.1
      ymax=ymax+dy
      ymin=ymin-dy
      xminsave=xmin
      xmaxsave=xmax
      yminsave=ymin
      ymaxsave=ymax

      call pgbeg(0,'/xs',1,1)
      call pgask(.false.)
      isym=23
      if (nplot.gt.2000) isym=-1
      diag=.false.

    2 continue
c      write(*,*) 'iaxes, just: ', iaxes, just
      call pgenv(xmin,xmax,ymin,ymax,just,iaxes)
      call pglab(lab1,lab2,lab3)
      call pgpt(nplot,xplot,yplot,isym)
      if (diag) then
          call pgmove(min(xmin,ymin),min(xmin,ymin))
          call pgdraw(max(xmax,ymax),max(xmax,ymax))
      endif

      xc=0.5*(xmin+xmax)
      yc=0.5*(ymin+ymax)
      
    3 call pgband(7,1,xc,yc,xc,yc,ch)
      if (ch.eq.' ') then
          xmin=xc
          ymin=yc
          xmax=(xmin+xmax)/2
          ymax=(ymin+ymax)/2
          call pgband(2,1,xmin,ymin,xmax,ymax,ch)
          goto 2
      endif
      if (ch.eq.'r') then
          xmax=xmaxsave
          ymax=ymaxsave
          xmin=xminsave
          ymin=yminsave
          goto 2
      endif
      if (ch.eq.'j') then
          just=1-just
          goto 2
      endif
      if (ch.eq.'a') then
          isaxis=mod(iaxes,10)
          iaxes=10*(iaxes/10) + (1-isaxis)
          goto 2
      endif
      if (ch.eq.'m') then
          isym=22-isym
          goto 2
      endif
      if (ch.eq.'i') then
          dx=(xmax-xmin)*0.1
          xmax=xmax-dx
          xmin=xmin+dx
          dy=(ymax-ymin)*0.1
          ymax=ymax-dy
          ymin=ymin+dy
          goto 2
      endif
      if (ch.eq.'o') then
          dx=(xmax-xmin)*0.1
          xmax=xmax+dx
          xmin=xmin-dx
          dy=(ymax-ymin)*0.1
          ymax=ymax+dy
          ymin=ymin-dy
          goto 2
      endif
      if (ch.eq.'d') then
          diag=.not.diag
          goto 2
      endif
      if (ch.eq.'s') then
          call pgend
          call pgbeg(0,'colplot.ps/ps',1,1)
          call pgslw(3)
          call pgenv(xmin,xmax,ymin,ymax,just,iaxes)
          call pglab(lab1,lab2,lab3)
          call pgpt(nplot,xplot,yplot,isym)
          if (diag) then
              call pgmove(min(xmin,ymin),min(xmin,ymin))
              call pgdraw(max(xmax,ymax),max(xmax,ymax))
          endif
          call pgend
          write(0,*) 'Wrote PS copy of the plot to colplot.ps'
          call pgbeg(0,'/xs',1,1)
          goto 2
      endif
      if (ch.eq.'?') then
          write(0,*) 'To zoom in, select bottom left with cursor,'
          write(0,*) '  hit SPACE, then select top-right, hit SPACE'
          write(0,*) 'To zoom in or out 20%, hit i or o'
          write(0,*) 'To reset plot window, hit r'
          write(0,*) 'To toggle justification or axes, hit j or a'
          write(0,*) 'To toggle marker type circle/dot, hit m'
          write(0,*) 'To toggle a diagonal line x=y, hit d'
          write(0,*) 'To save as PS file, hit s'
          write(0,*) 'To show this list, hit ?'
          write(0,*) 'To quit, hit q'
          goto 3
      endif
      if (ch.ne.'q') goto 3
      call pgend
      end
