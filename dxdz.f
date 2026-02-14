c/*
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
c          d_density_r_z_i_d_z      !density derivative from RZ spline
c----------------------------------------------------------------------
      real*8 FUNCTION dxdz(z,r,phi,i)
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
      integer isp,itype !local
      real*8 val,dvalz,dvalr,dvalphi !local
c-----externals
      real*8 thetapol,drhodz,dthetadz,psif,dvarddz,vardens,
     & ddnsrho,densrho,
     & d_density_r_z_i_d_z !spline art RZ mesh
     
     
![called by many...  - dddrz1, rside, transmit_coef_ox,
!find_rho_X, reflect, bound, boundc, plasmray, hotdervs]
      if((model_rho_dens.eq.7).and.(dens_read.eq.'enabled')) then 
         ![2024-08-14] Interpolate from dengrid_zrp(iz,ir,iphi) points
         !that are adjacent to the given point (z,r,phi).
         !Also - interpolate derivatives from grid points to (z,r,phi)
         itype=1 !=1 for density(and derivs); =2 for T; =3 for tpop
         isp=i !species number
         call interp_zrp(z,r,phi, val,dvalz,dvalr,dvalphi, isp,itype) !dn/d*
         !dxdr=   v(i)*dvalr
         !dxdphi= v(i)*dvalphi
         dxdz=   v(i)*dvalz !Here: in func dxdz(z,r,phi,i)
         return !Done here
      endif !((model_rho_dens.eq.7).and.(dens_read.eq.'enabled'))
     
cSAP090228
c      if(iboundb.ge.1)then
      if (rho.gt.1.d0-1.d-10) then
c----------------------------------------------------------
c       the point is outside LCFS
c---------------------------------------------------------------
cSAP090403 
         if(n_wall.gt.1) then
c-----------------------------------------------------------------
c          calculate density derivative using spline at RZ mesh
c-----------------------------------------------------
c          write(*,*)'in dense before density_r_z_i'
           dxdz=v(i)*d_density_r_z_i_d_z(z,r,phi,i) !derivative from RZ spline
         else
c---------------------------------------------------------------
cSAP090209
c----------------------------------------------------------------
c        calculate density derivative using formula versus small radius
c        and poloidal angle
c----------------------------------------------------------------
           theta_pol=thetapol(z,r) ! -pi <thetapol =<pi
           if (theta_pol.lt.0d0) then
             theta_pol=theta_pol+2*pi !pi< theta_pol<2pi
           endif

c          write(*,*)'dxdz z,r,phi,i,rho',rho,z,r,phi,rho

           call dens_rho_theta_LCFS(rho,theta_pol,i,
     &        dens_rho_theta,d_dens_rho_theta_d_rho,
     &        d_dens_rho_theta_d_theta)
cSm070325
c        dn_drho=ddnsrho(rho,i)
c        dxdz=v(i)*dn_drho*drhodz(z,r,phi)
           dxdz=v(i)*(d_dens_rho_theta_d_rho*drhodz(z,r,phi)
     &              +d_dens_rho_theta_d_theta*dthetadz(z,r))

c         write(*,*)'d_dens_rho_theta_d_rho,drhodz(z,r,phi)',
c     &              d_dens_rho_theta_d_rho,drhodz(z,r,phi)
c          write(*,*)'d_dens_rho_theta_d_rho*drhodz(z,r,phi)',
c     &               d_dens_rho_theta_d_rho*drhodz(z,r,phi)
c         write(*,*)'d_dens_rho_theta_d_theta,dthetadz(z,r)',
c     &              d_dens_rho_theta_d_theta,dthetadz(z,r)
c         write(*,*)'d_dens_rho_theta_d_theta*dthetadz(z,r)',
c     &              d_dens_rho_theta_d_theta*dthetadz(z,r)
c         write(*,*)'dxdz',dxdz
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
        dxdz=v(i)*dn_drho*dro_dpsi*dpdzd*(1.d0+vardens(z,r,phi))
        den=densrho(rho,i)

cSAP091017
c        if(den.lt.0.d0)then
c          den=0.d0
c	  write(*,*)'in dxdz den.lt.0.0 rho,i',rho,i
c        endif
c        dxdz=dxdz+v(i)*den*dvarddz(z,r,phi)
        if(den.lt.0.d0)then
          den=0.d0
	  write(*,*)'in dxdz den.lt.0.0 rho,i',rho,i
          dxdz=0.d0
        else
          dxdz=dxdz+v(i)*den*dvarddz(z,r,phi)
        endif

c       write(*,*)'in dxdz v(i),dn_drho,dro_dpsi,dpdzd,dxdz',
c     1                    v(i),dn_drho,dro_dpsi,dpdzd,dxdz

        

      endif
      return
      END
