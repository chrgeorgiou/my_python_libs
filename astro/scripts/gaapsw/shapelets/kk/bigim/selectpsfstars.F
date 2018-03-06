
c     reads a sextractor catalogue from stdin, displays it, and asks for
c     selection of the stellar locus. A filtered list s sent to stdout.

      parameter (m=mmcat)
      dimension x(m),y(m),f(m),fe(m),fw(m),a(m),b(m),th(m),id(m),iflg(m),amag(m)
      logical isstar(m)
      character ch

      amagzero=35
      open(10,file='magzero.dat',status='old',err=2)
      read(10,*) amagzero
      close(10)
    2 write(0,'('' Magnitude zero point:''f6.2)') amagzero

      npsf=0
      xmax=-1e30
      ymax=-1e30
      do i=1,m
 1       read(*,*,err=1,end=99) x(i),y(i),f(i),fe(i),fw(i),a(i),b(i),th(i),id(i),iflg(i)
         amag(i)=amagzero-2.5*log10(max(1.e-10,f(i)))
         xmax=max(xmax,x(i))
         ymax=max(ymax,y(i))
      enddo
 99   nobj=i-1
      write(0,*) 'Read ',nobj,' sources.'
      call pgbeg(0,'/xs',1,1)
      call pgask(.false.)
    6 call pgenv(0.,8.,30.,15.,0,0)
      call pgpt(nobj,fw,amag,5)
      call pgsci(3)
      ch=' '
      do while (ch.ne.'s')
          write(0,*) 'Select bottom left, top right of the PSF sequence'
          fw0=4
          am0=25
          call pgband(7,1,fw0,am0,fw0,am0,ch)
          fw1=fw0+1
          am1=am0-3
          call pgband(2, 1, fw0,am0, fw1,am1, ch)
          call pgmove(fw0,am0)
          call pgdraw(fw0,am1)
          call pgdraw(fw1,am1)
          call pgdraw(fw1,am0)
          call pgdraw(fw0,am0)
          write(0,*) 'Hit s to select these stars, other key to repeat.'
          call pgcurs(x,y,ch)
      enddo
      call pgsci(1)
      call pgenv(1.,xmax,1.,ymax,1,0)
      call pgpt(nobj,x,y,-1)
      call pgsci(3)
      do i=1,nobj
          if ((amag(i)-am1)*(amag(i)-am0).lt.0. .and. (fw(i)-fw1)*(fw(i)-fw0).lt.0.) then
              call pgpt(1,x(i),y(i),23)
          endif
      enddo
      call pgsci(1)
      write(0,*) 'Hit r to return to selection, other key to end.'
      call pgcurs(xmax/2,ymax/2,ch)
      if (ch.eq.'r') goto 6
      call pgend

c now display all selected objects
      nstars=0
      do i=1,nobj
          if ((amag(i)-am1)*(amag(i)-am0).lt.0. .and. (fw(i)-fw1)*(fw(i)-fw0).lt.0.) then
              isstar(i)=.true.
              nstars=nstars+1
          else
              isstar(i)=.false.
          endif
      enddo

c Now make PS plots of the final selection.
      call pgbeg(0,'psfsel.ps/ps',1,1)
      call pgslw(2)
      call pgenv(0.,8.,30.,15.,0,0)
      call pgpt(nobj,fw,amag,5)
      call pglab('FWHM','INSTR MAG','PSF OBJECTS')
      call pgmove(fw0,am0)
      call pgdraw(fw0,am1)
      call pgdraw(fw1,am1)
      call pgdraw(fw1,am0)
      call pgdraw(fw0,am0)
      call pgenv(1.,xmax,1.,ymax,1,0)
      call pgpt(nobj,x,y,-1)
      call pgsci(3)
      do i=1,nobj
          if (isstar(i))  call pgpt(1,x(i),y(i),23)
      enddo
      call pglab('X (pix)','Y (pix)','PSF OBJECTS')
      call pgend



      k=0
      do i=1,nobj
          if (isstar(i)) then
              write(*,'(f10.3,f11.3,g14.6,g12.5,f9.2,2f10.3,f6.1,i11,i6)')
     :            x(i),y(i),f(i),fe(i),fw(i),a(i),b(i),th(i),id(i),iflg(i)
              k=k+1
          endif
      enddo
      write(0,*) 'Selected ',k,' PSF stars.'
      end
