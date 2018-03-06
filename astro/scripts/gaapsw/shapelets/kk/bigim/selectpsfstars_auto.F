
c     reads a sextractor catalogue from stdin, displays it, and asks for
c     selection of the stellar locus. A filtered list s sent to stdout.
c
c     MODIFIED VERSION TO SELECT THE STARS AUTOMATICALLY
c
      dimension x(mmcat),y(mmcat),f(mmcat),fe(mmcat),fw(mmcat),
     :    a(mmcat),b(mmcat),th(mmcat),id(mmcat),iflg(mmcat),
     :    temp(mmcat),amag(mmcat),temp2(mmcat)
      logical isstar(mmcat)

c      call get2('inimage.fits',nx,ny,pix,mmpix,mmpix)

      amagzero=35
      open(10,file='magzero.dat',status='old',err=2)
      read(10,*) amagzero
      close(10)
    2 write(0,'('' Magnitude zero point:''f6.2)') amagzero

      npsf=0
      xmax=-1e30
      ymax=-1e30
      do i=1,mmcat
 1       read(*,*,err=1,end=99) x(i),y(i),f(i),fe(i),fw(i),a(i),b(i),th(i),id(i),iflg(i)
         amag(i)=amagzero-2.5*log10(max(1.e-10,f(i)))
         xmax=max(xmax,x(i))
         ymax=max(ymax,y(i))
      enddo
 99   nobj=i-1
      write(0,*) 'Read ',nobj,' sources.'

c     To find the stars, do the following:
c     first find the brightest 10% of sources, at least 100
c     of those find the smallest 20%, at least 20
c     anything at most 0.5 pixel wider than the median of these is taken.
c     of those, take away top and bottom 1/2 magnitude
c     parameters can be re-defined in starsel.par file

      topbright=0.1
      topsmall=0.2
      striptop=1
      snmin=25
      open(10,file='starsel.par',status='old',err=3)
      read(10,*) topbright,topsmall,striptop,snmin
      close(10)
    3 continue

      do i=1,nobj
          temp(i)=amag(i)
          temp2(i)=fw(i)
      enddo
      call sort2(nobj,temp,temp2)
      n1=nint(nobj*topbright)
      n1=min(nobj,max(n1,100))
      am0=temp(n1)
c am0 is now the bottom of the bright, compact sources
      call sort2(n1,temp2,temp)
      n2=max(3,nint(n1*topsmall))
c      n2=min(n1,max(n2,20))
      fw1=fmedian(n2,temp2)+0.5
c fw1 is now the median fwhm of the compact, bright sources
c      write(0,*) 'Nobj, n1, n2:',nobj,n1,n2
c      write(0,'('' PSF stars have fwhm below'',f6.2)') fw1
      fw0=(fw1-0.5)/3.
c fw0 is now a sensible lower limit to source size (1/3 of median)
      k=0
c now search for the brightest compact sources
      am1=am0
      do i=1,nobj
          if (fw(i).lt.fw1 .and. amag(i).lt.am1) then
              am1=amag(i)
          endif
      enddo
c strip off the top 1 magnitude to get rid of saturated objects
      am1=am1+striptop

      write(0,'('' PSF star mags above '',f6.2)')am1
      write(0,'('' PSF S/N at least '',f6.1)')snmin
      write(0,'('' PSF star fwhm between '',f5.2,'' and '',f5.2)')fw0,fw1


c now identify all selected objects
      nstars=0
      topmag=-1e10
      botmag=1e10
      do i=1,nobj
c          if ((amag(i)-am1)*(amag(i)-am0).lt.0. .and. (fw(i)-fw1)*(fw(i)-fw0).lt.0.) then
          if (amag(i).gt.am1 .and. (fw(i)-fw1)*(fw(i)-fw0).lt.0. .and. f(i).gt.snmin*fe(i)) then
              isstar(i)=.true.
              nstars=nstars+1
          else
              isstar(i)=.false.
          endif
          topmag=max(topmag,amag(i))
          botmag=min(botmag,amag(i))
      enddo

c Now make PS plots of the final selection.
      call pgbeg(0,'psfsel.ps/ps',1,1)
      call pgslw(2)
      call pgenv(0.,8.,topmag+0.5,botmag-0.5,0,0)
      call pgpt(nobj,fw,amag,5)
      call pglab('FWHM','INSTR MAG','PSF OBJECTS')
      call pgmove(fw0,am1+100)
      call pgdraw(fw0,am1)
      call pgdraw(fw1,am1)
      call pgdraw(fw1,am1+100)
      do i=1,nobj
          if (isstar(i))  call pgpt(1,fw(i),amag(i),23)
      enddo

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
