      
      parameter (mm012=5*(mm/2)-1)
      real psfcoef(mmpoly2,msh-2),roundpsfcoef(mmpoly2,msh-2)
      real wts(mm+1)
      
c     circularize a PSF map: azimuthally average all coefficients.  this
c     makes all odd orders zero, and for orders that are even in x orer
c     a and y order b they are averaged over all even coefs with the
c     same N=a+b according to weights
c
c     W_a = (N/2;a/2) sqrt (N;a)   where (X;Y) is the combinatorial coeff.
c
c     i.e., a set C_a --> sum_b C_b W_a/W_b

c read input PSF map from std input

      read(*,*) npsf,betapsf
      read(*,*) nfit,xmax,ymax
      read(*,*)
      do k=1,(npsf+1)*(npsf+2)/2
         read(*,*) (psfcoef(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(*,*)
      enddo
      write(0,*) 'Read psf model from psf.map'
      write(0,*) 'Order for shape, spatial variation; beta:',npsf,nfit,betapsf

c number of poly coefs per shapelet term; number of shapelet terms
      npoly=(nfit+1)*(nfit+2)/2
      nsh=(npsf+1)*(npsf+2)/2

      do kfit=1,npoly
          k=0
          do kk=1,nsh
              roundpsfcoef(kfit,kk)=0
          enddo
          do n=0,npsf
c if even then set wts for order n...                      
              do m=0,n
                  k=k+1
                  if (mod(n,2).eq.0 .and. mod(m,2).eq.0) then
                      roundpsfcoef(kfit,k)=roundpsfcoef(kfit,k)+
                      

              
