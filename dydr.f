c        ********************** dydr ***************************
c        *                      ----                           *
c        * this function calculates  the  derivative  of  y    *
c        * (being omega_cyclotron_i/omega) with  respect  to  r*
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
     1function dydr(z,r,phi,i)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include'one.i'

	dydr=w(i)*dbmdr

      return
      end
