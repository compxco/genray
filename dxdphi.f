c        ********************** dxdphi **********************
c        *                      ------                      *
c        * this function calculates  the  derivative  of  x *
c        * (being (omega_pl_i/omega)**2 ) with  respect  to *
c        * phi (i - type of particles)                      *
c        ****************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r phi - coordinates of the point where the derivative is !
c                 calculated.        		         	   !
c      i = 1 - electrons,                                          !
c        > 1 - ions.                                               !
c       rho is in common /one/
c------------------------------------------------------------------
      real*8 function dxdphi(z,r,phi,i)
c      implicit double precision (a-h,o-z)

      implicit none 

      include 'param.i'
      include 'one.i'
c-----input
      real*8 z,r,phi
      integer i     !number of plasma specie
c-----external
      real*8 d_density_r_z_i_d_phi,dense,dvarddph

c-----locals
      real*8 den

cSAP090206
c      den=densrho(rho,i)

c      write(*,*)'dxdphi before dense'
      if (rho.gt.1.d0-1.d-10) then   
c----------------------------------------------------------
c        the point is outside LCFS
c---------------------------------------------------------------
cSAP090403
         dxdphi=0.d0
         if(n_wall.gt.1) then
c-----------------------------------------------------------------
c          calculate density derivative using spline at RZ mesh
c-----------------------------------------------------
           dxdphi=v(i)*d_density_r_z_i_d_phi(z,r,phi,i) !derivative from RZ spline
         endif
      else
c----------------------------------------------------------
c        the point is inside LCFS
c---------------------------------------------------------------
         den=dense(z,r,phi,i)


         if(den.lt.0.d0)then
            den=0.d0
	    write(*,*)'in dxdphi den.lt.0.0 rho,i',rho,i
         endif
         dxdphi=0.d0
         dxdphi=v(i)*den*dvarddph(z,r,phi)

c         write(*,*)'in dxdphi v(i),den,dvarddph(z,r,phi),dxdphi',
c     &                        v(i),den,dvarddph(z,r,phi),dxdphi
      endif

      return
      end
