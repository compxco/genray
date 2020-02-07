
c the arrays from the input eqdskin (default="equilib.dat") file
      real*8 peqd,feqd,pres,ffpeqd,
     1ppeqd,qpsi,rlimit,zlimit,
     2rves,zves,rr,zz,flux

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
