c        **********************  cn  **************************
c        *                       --                           *
c        * this function calculates the absolute value of
c        * the wave refractive index                          *
c        *                                                    *
c        ******************************************************
c
c-----------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      r - r-coordinate of the point where the wave refractive     !
c          index is calculated.        		         	   !
c      cnz, cnr, cm - n_z, n_r and r*n_phi                         !
c------------------------------------------------------------------
	real*8 function cn(r,cnz,cnr,cm)
      implicit none !integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      real*8 r,cnz,cnr,cm !INPUT
      real*8 r2 ! local

       r2=r*r
	cn=dsqrt(cnz*cnz+cnr*cnr+cm*cm/r2)

      return
      end function cn
