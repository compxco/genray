c        ********************** dydz *************************** 
c        *                      ----                           *
c        * this function calculates  the  derivative  of  y    *
c        * (being omega_cyclotron_i/omega) with  respect  to  z*
c        * (i - type of particles)                             *
c        *******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r phi - coordinates of the point where the derivative is !
c                 calculated.        		         	   !
c      i = 1 - electrons,                                          !
c        > 1 - ions.                                               !
c      rho is from common one.i
c------------------------------------------------------------------
	double precision
     1function dydz(z,r,phi,i)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include'one.i'

	dydz=w(i)*dbmdz
      return
      end
