
      parameter (m=512,maxnx=30000)
      dimension zbufin(m*maxnx),zbufout(m*maxnx)
      dimension nzero(maxnx)
      logical zap(maxnx),zap2(maxnx)
      dimension naxlen(7)
      data naxlen /7*0/ 

c read in a weight frame and find all pixels where at least half of the surrounding pixels within ntrim are zero.
c then set all pixels in a ntrim x ntrim box around these pixels to zero.
c this has the effect of growing the edges of the empty regions of an image by a certain width 
c useful after convolution by a compact kernel.
c

      read(*,*) ntrim
      write(0,*) 'Trimming width: ',ntrim

      nzeromax=(2*ntrim+1)**2/2
c open image to be trimmed
      call openfits('weights.fits',51,nx,ny)
      if (nx.gt.maxnx) then
         write(0,*) ' *** IMAGE X dimension exceeds ',maxnx,' *** '
         stop
      endif
      naxlen(1)=nx
      naxlen(2)=ny

c open output image
      call copyfits_nodata('weights.fits',61,'trimmedweights.fits')

c read in the image in strips of m pixels wide.
c start ntrim below the edge of the image with all zeroes, to make sure edge effects counted properly.

      linein1=1-2*ntrim
      do while (linein1.lt.ny-2*ntrim)
         linein2=min(ny+2*ntrim,linein1+m-1)
         nlinein=linein2-linein1+1
         istat=igetsecpad(51,nx,ny,1,linein1,nx,nlinein,zbufin,maxnx*m,0.)
         write(0,*) 'Reading lines ',linein1,'-',linein2,'/',ny
c prepare the outbut buffer, which excludes first and last ntrim lines. In the end will exclude 2*ntrim lines top and bottom
         do i=1,nx*(nlinein-2*ntrim)
            zbufout(i)=zbufin(i+nx*ntrim)
         enddo
c calculate the number of zeroes within |dy|<=ntrim for each pixel of rom ntrim+1 (this is the first line to write out later)
c but do not add row ntrim*2 + 1 yet, this happens inside the loop
         do i=1,nx
            nzero(i)=0
         enddo
         do line=1,2*ntrim
            call addzeroes(zbufin,line,nx,nzero)
         enddo
c here nzero countains summed zero line counts per pixel for lines 1-ntrim to ntrim. So need to add ntrim+1 to get full set
         do line=ntrim+1,nlinein-ntrim
            call addzeroes(zbufin,line+ntrim,nx,nzero)
c figure out which pixels to use as centers for zap: those with more than half their neighbours equal to zero, AND themselves zero
c so set zap(i) to FALSE at all non-zero pixels in the line. These will stay false.
            do i=1,nx
               zap(i)= (zbufin(i+(line-1)*nx).eq.0)
            enddo
c nzero(i) is always the number of zero pixels within y +/- ntrim in the same column (i)
c nz is runing boxcar in x over +/- ntrim to see which pixels have more than 50% zero neighbours
c initialize at first pixel: count ntrim columns to left of column 1 as zeroes, and first ntrim+1 columns of image
            nz=(2*ntrim+1)*ntrim
            do i=1,ntrim+1
               nz=nz+nzero(i)
            enddo
            zap(1)=(nz.gt.nzeromax)
            do i=2,ntrim+1
               nz=nz-(2*ntrim+1)+nzero(i+ntrim)
               zap(i)=zap(i) .and. (nz.gt.nzeromax)
            enddo
            do i=ntrim+2,nx-ntrim
               nz=nz-nzero(i-ntrim-1)+nzero(i+ntrim)
               zap(i)=zap(i) .and. (nz.gt.nzeromax)
            enddo
            do i=nx-ntrim+1,nx
               nz=nz-nzero(i-ntrim-1)+(2*ntrim+1)
               zap(i)=zap(i) .and. (nz.gt.nzeromax)
            enddo
c zero the pixels in the output buffer within ntrim in x and y of those that have been flagged as edges
c first run a boxcar (in x) over zap to make zap2
            do i=1,ntrim
               zap2(i)=.true.
            enddo
            do i=ntrim+1,nx-ntrim
               zap2(i)=.false.
               do j=-ntrim,ntrim
                  zap2(i)=zap2(i) .or. zap(i-j)
               enddo
            enddo
            do i=nx-ntrim+1,nx
               zap2(i)=.true.
            enddo
c then zero all pixels with zap2=TRUE within ntrim in y
            do i=1,nx
               if (zap2(i)) then
                  do lout=max(1,line-2*ntrim),min(nlinein-2*ntrim,line)
                     zbufout(i+(lout-1)*nx)=0
                  enddo
               endif
            enddo
c update the counts of nearby wt=zero pixels per column by subtracting the lowest line used (line-ntrim).
c will need to add new line (line+ntrim+1) in the next iteration (done in beginning of the loop)
            call subzeroes(zbufin,line-ntrim,nx,nzero)
c ready to move on to the next line
         enddo
c write out the completed section and get ready to read the next one
         call putsec(61,zbufout(1+ntrim*nx),nx,nlinein-4*ntrim,1,linein1+2*ntrim,naxlen)
         write(0,*) 'Writing lines',linein1+2*ntrim,'-',linein1+nlinein-2*ntrim-1
         linein1=linein1+m-4*ntrim
      enddo
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


