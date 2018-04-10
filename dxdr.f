c/*
c        ********************** dxdr ************************
c        *                      ----                        *
c        * this function calculates  the  derivative  of  x *
c        * (being (omega_pl_i/omega)**2 ) with  respect  to *
c        * r (i - type of particles)                        *
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
c      rho is from common/one/                                     !
c------------------------------------------------------------------
c         uses
c          v,dense0,denseb,rn1de,rn2de,rho,idens   from common 'one'
c           psilim,psimag            from common 'three'
c           functions psif, drhopsi, ddnsdrho
c----------------------------------------------------------------------
      double precision FUNCTION dxdr_old(z,r,phi,i)
      IMPLICIT double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'

c      write(*,*)' dxdr.f z,r,phi,i,iboundb',z,r,phi,i,iboundb
c      write(*,*)'rho',rho
      if(iboundb.ge.1)then
c       the point is outside the plasma
c	if(iboundb.eq.2) dxdr=0.d0
c   	if(iboundb.eq.1) then
c	  dn_drho=ddnsrho(rho,i)
c          delr=0.5d0*(rmax-rmin)
c          if (r.le.rzmin) drhodr=-1.d0/delr
c          if (r.ge.rzmin) drhodr=1.d0/delr
c          dxdr=v(i)*dn_drho*drhodr
c	endif
cSm070325
        dn_drho=ddnsrho(rho,i)
        dxdr=v(i)*dn_drho*drhodr(z,r,phi)
      else
        psi_xr=psif(z,r)
        dro_dpsi=drhopsi(psi_xr)
c------------------------------------------------------------
c       the spline form only  !
c------------------------------------------------------------
c        write(*,*)'dxdr.f rho,i',rho,i

	dn_drho=ddnsrho(rho,i)

c        write(*,*)'dxdr.f dn_drho,dro_dpsi,dpdrd',
c     &                    dn_drho,dro_dpsi,dpdrd

        dxdr=v(i)*dn_drho*dro_dpsi*dpdrd*(1.d0+vardens(z,r,phi))
        den=densrho(rho,i)
        if(den.lt.0.d0)then
           den=0.d0
           write(*,*)'in dxdr den.lt.0.0 rho,i',rho,i
        endif
        dxdr=dxdr+v(i)*den*dvarddr(z,r,phi)
      endif

      return
      END


cSAP090209
c        ********************** dxdz ************************
c        *                      ----                        *
c        * this function calculates  the  derivative  of  x *
c        * (being (omega_pl_i/omega)**2 ) with  respect  to *
c        * z (i - type of particles)                        *
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
c------------------------------------------------------------------
c         uses
c          v,dense0,denseb,rn1de,rn2de,rho,idens   from common 'one'
c          psilim, psimag            from common 'three'
c          functions psif, drhopsi, ddnsdrho
c          d_density_r_z_i_d_r !density derivative from RZ spline
c----------------------------------------------------------------------
      double precision FUNCTION dxdr(z,r,phi,i)
c      IMPLICIT double precision (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
      include 'five.i'
c-----input
      real*8 z,r,phi !space coordinates
      integer i      !the plasma species number
c-----locals
      real*8  
     &  theta_pol,      !poloidal angle [radians] pi< theta_pol<2pi 
     &  dens_rho_theta,d_dens_rho_theta_d_rho,
     &  d_dens_rho_theta_d_theta,
     &  den,dn_drho,dro_dpsi,psi_xr,drhopsi
c-----externals
      real*8 thetapol,drhodr,dthetadr,psif,dvarddr,vardens,
     & ddnsrho,densrho,
     & d_density_r_z_i_d_r !density derivative from spline at RZ mesh
cSAP090228
c      if(iboundb.ge.1)then
      if (rho.gt.1.d0-1.d-10) then
c----------------------------------------------------------
c        the point is outside LCFS
c---------------------------------------------------------------
cSAP090403
         if(n_wall.gt.1) then
c-----------------------------------------------------------------
c          calculate density derivative using spline at RZ mesh
c----------------------------------------------------------------
           dxdr=v(i)*d_density_r_z_i_d_r(z,r,phi,i) !derivative from RZ spline
         else
cSAP090209
c----------------------------------------------------------------
c          calculate density derivative using formula versus small radius
c          and poloidal angle
c----------------------------------------------------------------
           theta_pol=thetapol(z,r) ! -pi <thetapol =<pi
           if (theta_pol.lt.0d0) then
              theta_pol=theta_pol+2*pi !pi< theta_pol<2pi
           endif

c          write(*,*)'dxdr z,r,phi,i,rho',z,r,phi,i,rho

           call dens_rho_theta_LCFS(rho,theta_pol,i,
     &        dens_rho_theta,d_dens_rho_theta_d_rho,
     &        d_dens_rho_theta_d_theta)
cSm070325
c        dn_drho=ddnsrho(rho,i)
c        dxdr=v(i)*dn_drho*drhodr(z,r,phi)
           dxdr=v(i)*(d_dens_rho_theta_d_rho*drhodr(z,r,phi)
     &              +d_dens_rho_theta_d_theta*dthetadr(z,r))
        
c         write(*,*)'d_dens_rho_theta_d_rho,drhodr(z,r,phi)',
c     &              d_dens_rho_theta_d_rho,drhodr(z,r,phi)
c         write(*,*)'d_dens_rho_theta_d_rho*drhodr(z,r,phi)',
c     &              d_dens_rho_theta_d_rho*drhodr(z,r,phi)

c         write(*,*)'d_dens_rho_theta_d_theta,dthetadr(z,r)',
c     &              d_dens_rho_theta_d_theta,dthetadr(z,r)
c          write(*,*)'d_dens_rho_theta_d_theta*dthetadr(z,r)',
c     &               d_dens_rho_theta_d_theta*dthetadr(z,r)

c         write(*,*)'dxdr',dxdr
         endif !n_wall.ge.1
      else
c-----------------------------------------------------------
c       the point inside LCFS
c-----------------------------------------------------------
        psi_xr=psif(z,r)
        dro_dpsi=drhopsi(psi_xr)
c------------------------------------------------------------
c       spline form
c------------------------------------------------------------
        dn_drho=ddnsrho(rho,i)
        dxdr=v(i)*dn_drho*dro_dpsi*dpdrd*(1.d0+vardens(z,r,phi))
        den=densrho(rho,i)

cSAP091017
c        if(den.lt.0.d0)then
c          den=0.d0
c	  write(*,*)'in dxdz den.lt.0.0 rho,i',rho,i
c        endif
c        dxdr=dxdr+v(i)*den*dvarddr(z,r,phi)

        if(den.lt.0.d0)then
           den=0.d0
	   write(*,*)'in dxdz den.lt.0.0 rho,i',rho,i
           dxdr=0.d0
        else
           dxdr=dxdr+v(i)*den*dvarddr(z,r,phi)
        endif


c       write(*,*)'in dxdr v(i),dn_drho,dro_dpsi,dpdrd,dxdr',
c     1                    v(i),dn_drho,dro_dpsi,dpdrd,dxdr
      endif
      return
      END
