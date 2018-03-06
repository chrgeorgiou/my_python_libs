c     go through a raw object catalogue, flag all sources which have pixels of
c     zero weight within radius 1*FWHM.
c     weights are obtained from weight image "weights.fits"

      parameter (maxstamp=800,maxstamp2=maxstamp*maxstamp)
      dimension wtsec(maxstamp2)
      logical haszerowt

      call openfits('weights.fits',51,nx,ny)

      do k=1,maxstamp2
          wtsec(k)=0
      enddo

      nflagged=0
      nread=0

c     go through shapelet catalogue from std input
    2 read(*,*,err=2,end=1) x,y,f,fe,fw,a,b,th,id,iflg
      nread=nread+1
      fitradius=max(10.,fw)
c first catch those within fw of the edge
      if ((x-fw-1)*(nx-fw-x).lt.0. .or. (y-fw-1)*(ny-fw-y).lt.0.) then
         iflg=iflg+4000
         nflagged=nflagged+1
         goto 3
      endif
      npix=int(2*fitradius)+3
      ix1=int(x-npix/2)
      ix1=max(1,min(ix1,nx-npix+1))
      iy1=int(y-npix/2)
      iy1=max(1,min(iy1,ny-npix+1))
      isok=igetsec(51,nx,ny,ix1,iy1,npix,npix,wtsec,maxstamp2)
      if (isok.eq.0) then
          if (haszerowt(wtsec,0.,npix,x-ix1+1,y-iy1+1,fitradius)) then
              iflg=iflg+4000
c              write(0,*) 'Flagged object ID,x,y:',id,x,y
              nflagged=nflagged+1
              goto 3
          endif
      else
          iflg=iflg+4000
c          write(0,*) 'Flagged object ID,x,y:',id,x,y
          nflagged=nflagged+1
      endif
    3 continue
      if (mod(nread,100).eq.0) write(0,*) 'Flagged ',nflagged,' / ',nread
      write(*,1001) x,y,f,fe,fw,a,b,th,id,iflg
 1001 format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i10,i5)
      goto 2

    1 continue
      write(0,*) 'Flagged ',nflagged,' / ',nread
      end


      logical function haszerowt(wt,thresh,n,xc,yc,r)
      dimension wt(n,n)
c scan disk of radius r around xc,yc in wt(n,n) array for low wt
c any wt below thresh results in flag
      haszerowt=.false.
      r2=r*r
      wmin=1e32
      wmax=-wmin
      do iy=max(1,int(yc-r)),min(n,int(yc+r))
          dy=iy-yc
          do ix=max(1,int(xc-r)),min(n,int(xc+r))
              dx=ix-xc
              dr2=dx*dx+dy*dy
              if (dr2.le.r2 .and. wt(ix,iy).le.thresh) then
                  haszerowt=.true.
                  return
              endif
          enddo
      enddo
      return
      end
