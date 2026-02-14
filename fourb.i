
c the arrays from the input eqdskin (default="equilib.dat") file
      real*8 peqd,feqd,pres,ffpeqd,
     1 ppeqd,qpsi,rlimit,zlimit,req,xeq,yeq,zeq,
     2 rves,zves,rr,zz,flux

c	Set up some variables for re-gridding of input f(1:nveqd), etc.,
c       arrays onto a equispaced mesh of dimension nxeqd.
c       BobH020822.

      integer iworkka,itabl,i1p,n_wall_add,n_limiter_add
      parameter(iworkka=3*nxeqda+1)
      common/fourb_int/ n_limiter_add, itabl(3), i1p(2),
     &                  n_wall_add(0:max_limiters_a)
      !YuP[2020-01] Collected all integer variables into a separate block,
      !for better alignment

      
      
       
      real*8 psiar,d2feqd,d2pres,d2ffpeqd,d2ppeqd,d2qpsi,workk,r8temp,
     &tabl,
     & thetapol_wall,thetapol_limiter,
     & rho_wall,rho_limiter,
     &r_wall_add,z_wall_add,                  
     &theta_pol_limit,        
     &rho_geom_limit    

      real*8, pointer ::
     &dens_r_z_in(:,:,:),         !(1:nxeqd_add,1:nyeqd_add,1:nbulk)
     &temperature_r_z_in(:,:,:),  !(1:nxeqd_add,1:nyeqd_add,1:nbulk)
     &zeff_r_z_in(:,:),           !(1:nxeqd_add,1:nyeqd_add)
cSAP111012
     &rr_add(:),                  !(1:nxeqd_add)
     &zz_add(:),                  !(1:nyeqd_add)
     &distance_to_wall(:,:,:),    !(1:nxeqd_add,1:nyeqd_add,0:max_limiters)
                                  !density at rz mesh
     &density_r_z(:,:,:,:),       !(1:nxeqd_add,1:nyeqd_add,1:nbulk,
                                  !0:max_limiters)
                                  !second derivatives
     &density_r_z_rr(:,:,:,:),    !(1:nxeqd_add,1:nyeqd_add,1:nbulk,
                                  !0:max_limiters)
                                  !second derivatives
     &density_r_z_zz(:,:,:,:),    !(1:nxeqd_add,1:nyeqd_add,1:nbulk,!
                                  !0:max_limiters)
                                  !4th derivatives 
     &density_r_z_rrzz(:,:,:,:),  !(1:nxeqd_add,1:nyeqd_add,1:nbulk,
                                  !0:max_limiters) 
                                  !temperature at rz mesh
     &temperature_r_z(:,:,:),     !(1:nxeqd_add,1:nyeqd_add,1:nbulk)
                                  !second derivatives  
     &temperature_r_z_rr(:,:,:),  !(1:nxeqd_add,1:nyeqd_add,1:nbulk)
                                  !second derivatives
     &temperature_r_z_zz(:,:,:),  !(1:nxeqd_add,1:nyeqd_add,1:nbulk)
                                  !4th derivatives   
     &temperature_r_z_rrzz(:,:,:), !(1:nxeqd_add,nyeqd_add,1:nbulk)
     &zeff_r_z(:,:),              !(1:nxeqd_add,1:nyeqd_add) zeff at rz mesh
     &zeff_r_z_rr(:,:),           !(1:nxeqd_add,1:nyeqd_add) second derivatives
     &zeff_r_z_zz(:,:),           !(1:nxeqd_add,1:nyeqd_add) second derivatives
     &zeff_r_z_rrzz(:,:)          !(1:nxeqd_add,1:nyeqd_add) 4th derivatives
      common/fourb/ feqd(nxeqda),pres(nxeqda),ffpeqd(nxeqda),
     1             ppeqd(nxeqda),peqd(nxeqda,nyeqda),
     2             qpsi(nxeqda),
     3             rlimit(nlimit),zlimit(nlimit),
     4             rves(nves),zves(nves),
     5             rr(nxeqda),zz(nyeqda),
     5             req(nxeqda),xeq(2*nxeqda),yeq(2*nxeqda),zeq(nyeqda),
     6             flux(nxeqda),
     &             theta_pol_limit(nlimit),        
     &             rho_geom_limit(nlimit),           
c------BobH02022
     &             psiar(nxeqda),d2feqd(nxeqda),d2pres(nxeqda),
     1             d2ffpeqd(nxeqda),d2ppeqd(nxeqda),d2qpsi(nxeqda),
     2            r8temp(nxeqda),workk(iworkka),tabl(3),
c-----for wall and limiter
     &             thetapol_wall(n_wall_a),
     &             thetapol_limiter(n_limiter_a,max_limiters_a),
     &             rho_wall(n_wall_a),
     &             rho_limiter(n_limiter_a,max_limiters_a),
     &             r_wall_add( n_wall_add_a,0:max_limiters_a),
     &             z_wall_add( n_wall_add_a,0:max_limiters_a),
c-----for density fall near the wall
c-------pointers at rz mesh 
     &             rr_add,
     &             zz_add,
     &             distance_to_wall,
     &             density_r_z,
     &             density_r_z_rr,
     &             density_r_z_zz,
     &             density_r_z_rrzz,
     &             temperature_r_z,  
     &             temperature_r_z_rr,
     &             temperature_r_z_zz,
     &             temperature_r_z_rrzz,  
     &             zeff_r_z,  
     &             zeff_r_z_rr,
     &             zeff_r_z_zz,
     &             zeff_r_z_rrzz,  
c-------pointers for density and zeff at rz mesh 
     &dens_r_z_in,  
     &temperature_r_z_in,  
     &zeff_r_z_in

c
c  nxeqda must be >= nxeqd
c  nyeqda must be >= nyeqd
c
c----------------------------------------------------------------
c thetapol_wall,       ! poloidal angles [radians] of
c thetapol_limiter     ! wall and limiter points     
c rho_wall,rho_limiter !small radius
c
!=======================================================================
!YuP[2024-08-14] added, for new option of density(z,r,phi) profile
!and temperature(z,r,phi) (for each species).
c Input data from dendsk file (model_rho_dens.eq.7):
      integer nzden,nrden,nphiden !Grid sizes, read from file
      !Note: Assume same (z,r,phi)-grid for density and temperature!
      real*4 denmin,denmax
      real*4 tempmin,tempmax
      
      !Note for zden,rden,phiden: Same grid for density and temperature!
      real*4 rdenmin,rdenmax, phidenmin,phidenmax, zdenmin,zdenmax
      real*4 drden,dphiden,dzden  !grid spacings
      real*4,pointer :: zden(:)   !(nzden) !grid in Z
      real*4,pointer :: rden(:)   !(nrden) !grid in R
      real*4,pointer :: phiden(:) !(nphiden) !grid in phi
      
      !Density and derivatives of density:
      !save in single precision - could be a huge array:
      real*4,pointer :: dengrid_zrp(:,:,:,:) !(nzden,nrden,nphiden,nbulk) 
      real*4,pointer :: dnzrp_dz(:,:,:,:)    !(nzden,nrden,nphiden,nbulk)
      real*4,pointer :: dnzrp_dr(:,:,:,:)    !(nzden,nrden,nphiden,nbulk)
      real*4,pointer :: dnzrp_dphi(:,:,:,:)  !(nzden,nrden,nphiden,nbulk)     
      common/dendsk_zrp/
     1         nzden,nrden,nphiden,
     2  zdenmin,zdenmax, rdenmin,rdenmax, phidenmin,phidenmax,
     2               dzden,drden,dphiden,
     3               zden, ! uniform-grid
     3               rden, ! uniform-grid
     3               phiden, ! uniform-grid
     1      denmin,denmax, !To be found: density min/max 
     4      dengrid_zrp, dnzrp_dz,dnzrp_dr,dnzrp_dphi
     
      !Temperature and derivatives of T:  
      !save in single precision - could be huge arrays:
      real*4,pointer :: tempgrid_zrp(:,:,:,:) !(nzden,nrden,nphiden,nbulk) 
      real*4,pointer :: dtzrp_dz(:,:,:,:)     !(nzden,nrden,nphiden,nbulk)
      real*4,pointer :: dtzrp_dr(:,:,:,:)     !(nzden,nrden,nphiden,nbulk)
      real*4,pointer :: dtzrp_dphi(:,:,:,:)   !(nzden,nrden,nphiden,nbulk)
      common/tempdsk_zrp/
     1      tempmin,tempmax, !To be found: T min/max [keV]
     4      tempgrid_zrp,dtzrp_dz,dtzrp_dr,dtzrp_dphi
     
      ![2024-08-14] Added, for Tpop=T_perp/T_parallel and derivatives of Tpop
      real*4 tpopmin,tpopmax
      !save in single precision - could be huge arrays:
      real*4,pointer :: tpopgrid_zrp(:,:,:,:)  !(nzden,nrden,nphiden,nbulk) 
      real*4,pointer :: dtpopzrp_dz(:,:,:,:)   !(nzden,nrden,nphiden,nbulk)
      real*4,pointer :: dtpopzrp_dr(:,:,:,:)   !(nzden,nrden,nphiden,nbulk)
      real*4,pointer :: dtpopzrp_dphi(:,:,:,:) !(nzden,nrden,nphiden,nbulk)
      common/tpopdsk_zrp/
     1      tpopmin,tpopmax,    !To be found: Tpop min/max
     4      tpopgrid_zrp,dtpopzrp_dz,dtpopzrp_dr,dtpopzrp_dphi