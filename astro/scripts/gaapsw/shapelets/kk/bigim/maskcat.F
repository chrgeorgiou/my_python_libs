c     go through a shapelet catalogue and flag all sources which have pixels of
c     zero weight within the fitting radius somewhere.
c     weights are obtained from weight image "weights.fits"

      parameter (maxstamp=800,maxstamp2=maxstamp*maxstamp)
      dimension wtsec(maxstamp2)
      dimension sh(msh)
      logical haszerowt

      call openfits('weights.fits',51,nx,ny)
      do k=1,maxstamp2
          wtsec(k)=0
      enddo
c     go through shapelet catalogue from std input
    2 read(*,*,err=1,end=1) x,y,f,fe,fw,a,b,th,id,iflg
      read(*,*) sum,dx,dy,bg,fitradius,sherr
      read(*,*) beta,n
      read(*,*) (sh(k),k=1,(n+1)*(n+2)/2)
      npix=int(2*fitradius)+3
      ix1=int(x-npix/2)
      ix1=max(1,min(ix1,nx-npix+1))
      iy1=int(y-npix/2)
      iy1=max(1,min(iy1,ny-npix+1))
      isok=igetsec(51,nx,ny,ix1,iy1,npix,npix,wtsec,maxstamp2)
      if (isok.eq.0) then
          if (haszerowt(wtsec,npix,x-ix1+1,y-iy1+1,fitradius)) then
              iflg=iflg+4000
              write(0,*) 'Flagged object ID,x,y:',id,x,y
              goto 3
          endif
      else
          iflg=iflg+4000
          write(0,*) 'Flagged object ID,x,y:',id,x,y
      endif


    3 continue
      write(*,1001) x,y,f,fe,fw,a,b,th,id,iflg
 1001 format(2f10.3,g14.6,g12.5,f8.2,2f8.3,f6.1,i10,i5)
      write(*,1002) sum,dx,dy,bg,fitradius,sherr
 1002 format (g14.7,2f10.5,g15.7,f10.5,g15.6)
      write(*,*) beta,n
      write(*,*) (sh(k),k=1,(n+1)*(n+2)/2)
      write(*,*)
      goto 2

    1 continue
      write(*,*) 
      end


      logical function haszerowt(wt,n,xc,yc,r)
      dimension wt(n,n)
c scan disk of radius r around xc,yc in wt(n,n) array for nonpositive values
      haszerowt=.false.
      r2=r*r
      do iy=1,n
          dy=iy-yc
          do ix=1,n
              dx=ix-xc
              dr2=dx*dx+dy*dy
              if (dr2.le.r2 .and. wt(ix,iy).le.0.) then
                  haszerowt=.true.
                  return
              endif
          enddo
      enddo
      return
      end
