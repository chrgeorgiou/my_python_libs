
      parameter (m=512,maxnx=30000)
      dimension zbufin(m*maxnx),zbufout(m*maxnx)
      logical zap(maxnx),zap2(maxnx)
      dimension naxlen(7)
      data naxlen /7*0/ 

c read in a weight frame and find all pixels where the weight is zero.
c then set all pixels in a ngrow x ngrow box around these pixels to zero.
c
c useful after convolution of a data frame by a compact kernel.
c

      read(*,*) ngrow
      write(0,*) 'Growing zero-weight pixels by: ',ngrow

      nzeromax=(2*ngrow+1)**2/2
c open input image 
      call openfits('weights.fits',51,nx,ny)
      if (nx.gt.maxnx) then
         write(0,*) ' *** IMAGE X dimension exceeds ',maxnx,' *** '
         stop
      endif
      naxlen(1)=nx
      naxlen(2)=ny

c open output image
      call copyfits_nodata('weights.fits',61,'trimmedweights.fits')

c write out first ngrow lines - these are all zero      
      do i=1,nx*ngrow
          zbufout(i)=0
      enddo
      call putsec(61,zbufout,nx,ngrow,1,1,naxlen)
      write(0,*) 'Writing lines',1,'-',ngrow
c read in the image in strips of m pixels wide.
      linein1=1
      do while (linein1.lt.ny-ngrow)
         linein2=min(ny,linein1+m-1)
         nlinein=linein2-linein1+1
         istat=igetsec(51,nx,ny,1,linein1,nx,nlinein,zbufin,maxnx*m)
c         write(0,*) 'Reading lines ',linein1,'-',linein2,'/',ny
c prepare the outbut buffer, which excludes first and last ngrow lines. 
         do i=1,nx*(nlinein-2*ngrow)
            zbufout(i)=zbufin(i+nx*ngrow)
         enddo
c go through all lines in the buffer: flag the zero-weight pixels
         do line=1,nlinein
            do i=1,nx
               zap(i)= (zbufin(i+(line-1)*nx).eq.0)
            enddo
c zero the pixels in the output buffer within ngrow in x and y of those
c that have been flagged
c first run a boxcar (in x) over zap to make zap2
            do i=1,ngrow
                zap2(i)=.true.
            enddo
c nz is running number of zeroes between pixel i-ngrow and i+ngrow-1
c so add 1 if zap(i+ngrow)=T; set zap2(i) if nz>0; subtract 1 if zap(i-ngrow)=T
            nz=0
            do i=1,2*ngrow
                if (zap(i)) nz=nz+1
            enddo
            do i=ngrow+1,nx-ngrow
                if (zap(i+ngrow)) nz=nz+1
                zap2(i)= (nz.gt.0) 
                if (zap(i-ngrow)) nz=nz-1
            enddo
            do i=nx-ngrow+1,nx
                zap2(i)=.true.
            enddo
c then zero all pixels with zap2=TRUE within ngrow in y - avoid edges of the image
            lout1=max(1,line-2*ngrow)
            lout2=min(nlinein-2*ngrow,line)
            do i=1,nx
               if (zap2(i)) then
                  do lout=lout1,lout2
                     zbufout(i+(lout-1)*nx)=0
                  enddo
               endif
            enddo
c ready to move on to the next line
         enddo
c write out the completed section and get ready to read the next one
         call putsec(61,zbufout,nx,nlinein-2*ngrow,1,linein1+ngrow,naxlen)
         write(0,*) 'Writing lines',linein1+ngrow,'-',linein1+nlinein-ngrow-1,'/',ny
         linein1=linein1+m-2*ngrow
      enddo
      do i=1,nx*ngrow
          zbufout(i)=0
      enddo
      call putsec(61,zbufout,nx,ngrow,1,ny-ngrow+1,naxlen)
      write(0,*) 'Writing lines',ny-ngrow+1,'-',ny,'/',ny
c Done!
      call closefits(51)
      call closefits(61)
      end


      subroutine addzeroes(z,line,n,nz)
c add 1 to every nz(i) counter whenever z(i,line) = 0
      dimension z(n,line),nz(n)
      do i=1,n
         if (z(i,line).eq.0) nz(i)=nz(i)+1
      enddo
      return
      end

      subroutine subzeroes(z,line,n,nz)
c add 1 to every nz(i) counter whenever z(i,line) = 0
      dimension z(n,line),nz(n)
      do i=1,n
         if (z(i,line).eq.0) nz(i)=nz(i)-1
      enddo
      return
      end


