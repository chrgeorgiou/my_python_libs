
      subroutine readorders(n,nfit)
      open(20,file='orders.par',status='old')
      read(20,*) n
      read(20,*) nfit
      close(20)
c      write(0,*) 'Order for shapelets, spatial interpolation: ',n,nfit
      return
      end



      logical function readpsfmap(filename,psf,nsh,beta,nfit,xmax,ymax)
c read in a shapelet psfmap called 'filename'.
c argument is TRUE on success, FALSE on failure
      character *(*) filename
      real psf(mmpoly2,msh-2)
      readpsfmap=.false.
      open(23,file=filename,status='old',err=2)
 1    read(23,*,err=1) nsh,beta
      read(23,*) nfit,xmax,ymax
      read(23,*) 
      do k=1,(nsh+1)*(nsh+2)/2
         read(23,*,err=2) (psf(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(23,*)
      enddo
      readpsfmap=.true.
 2    return
      end
