      double precision function psi0_minus_psi_zr(r)
c--------------------------------------------------
c     calculates function: psi0_minus_psi_zr=psi_0-psi(z_0,r)    
c--------------------------------------------------
      implicit none
c-----input
       double precision
     & r          ! major radius

       double precision
     & psi_0,     !poloidal flux from /common_psi_m_psi_zr/ 
     & z_0        !vertical coordinate from /common_psi_m_psi_zr/ 
      common /common_psi_m_psi_zr/ psi_0,z_0 ! input data
c-----externals
      double precision psif
     
      psi0_minus_psi_zr=psi_0-psif(z_0,r)

      return
      end

      subroutine r_min_r_max_from_psi_z(psi0,z0,accuracy,
     & rmin_psi_z,rmax_psi_z)
c----------------------------------------------------------
c     calculates two major radus points: rmin_psi_z,rmax_psi_z
c     at the flux surface with the given poloidal flux: psi
c     at the given vertical coordinate: z
c----------------------------------------------------------
      implicit none
c-----input
      double precision
     & psi0,     ! poloidal flux
     & z0,       ! vertical coordinate
     & accuracy  ! accuracy of root calciulations
     

      include 'param.i'
      include 'three.i'       
      include 'five.i'
c-----output
      double precision rmin_psi_z,rmax_psi_z
c-----local
      double precision X1,X2
        
c-----externals
      double precision  rtbis,psi_m_psi_zr
      external psi_m_psi_zr
      double precision
     & psi_0,     !poloidal flux in /common_psi_m_psi_zr/ 
     & z_0        !vertical coordinate in /common_psi_m_psi_zr/ 
      common /common_psi_m_psi_zr/ psi_0,z_0 ! input data

c----------------------------------------------------------
c     put parameters to common / common_psi_m_psi_zr/
c----------------------------------------------------------
      psi_0=psi0
      z_0=z_0
c----------------------------------------------------------

      X1=xma
      X2=rmax*(1.d0+1.d-4)
      rmax_psi_z=rtbis(psi_m_psi_zr,x1,x2,accuracy)

      X1=rmin*(1.d0-1.d-4)
      X2=xma
      rmin_psi_z=rtbis(psi_m_psi_zr,x1,x2,accuracy)

      return
      end

