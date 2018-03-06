c read salient shapelet expansion numbers, write out as single line per source
      write(*,*) '#   x     y    flx   flxerr   fwhm    ID   FLG Shflx  dx    dy    bg    Sherr   beta'
      do i=1,mmcat
 1       read(*,*,err=1,end=2) x,y,f,fe,fw,a,b,pa,id,iflg,
     :        shum,dx,dy,bg,rfit,sherr,beta,nsh,
     :        (sh,k=1,(nsh+1)*(nsh+2)/2)
         write(*,1001) x,y,f,fe,fw,id,iflg,shum,dx,dy,bg,sherr,beta
      enddo
 2    n=i-1
 1001 format (2f12.3,2g15.5,f11.3,i10,i6,g15.5,2f10.3,g15.5,g15.5,f11.3)
      write(0,*) 'Wrote out data for ',n,' sources.'
      end
