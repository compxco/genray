
c the arrays from the input eqdskin (default="equilib.dat") file
      real*8 peqd,feqd,pres,ffpeqd,
     1ppeqd,qpsi,rlimit,zlimit,
     2rves,zves,rr,zz,flux

c	Set up some variables for re-gridding of input f(1:nveqd), etc.,
c       arrays onto a equispaced mesh of dimension nxeqd.
c       BobH020822.
      integer iworkka,itabl,i1p,n_wall_add,n_limiter_add

      parameter(iworkka=3*nxeqda+1)
      real*8 psiar,d2feqd,d2pres,d2ffpeqd,d2ppeqd,d2qpsi,workk,r8temp,
     &tabl,
     & thetapol_wall,thetapol_limiter,
     & rho_wall,rho_limiter,
     &r_wall_add,z_wall_add,
     &rr_add,zz_add,distance_to_wall,
     &density_r_z,                   !density at rz mesh
     &density_r_z_rr,density_r_z_zz, !second derivatives
     &density_r_z_rrzz,               !4th derivatives 
     &temperature_r_z,                   !density at rz mesh
     &dtemperature_r_z_rr,temperature_r_z_zz, !second derivatives
     &temperature_r_z_rrzz,               !4th derivatives    
     &theta_pol_limit,        
     &rho_geom_limit    

      real*8, pointer ::
     &dens_rz_in(:,:,:),      !(1:nxeqd_add,1:nyeqd_add,1:nbulk)
     &zeff_rz_in(:,:)         !(1:nxeqd_add,1:nyeqd_add)

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
     2            r8temp(nxeqda),i1p(2),workk(iworkka),tabl(3),itabl(3),
c-----for wall and limiter
     &             thetapol_wall(n_wall_a),
     &             thetapol_limiter(n_limiter_a,max_limiters_a),
     &             rho_wall(n_wall_a),
     &             rho_limiter(n_limiter_a,max_limiters_a),
     &             r_wall_add( n_wall_add_a,0:max_limiters_a),
     &             z_wall_add( n_wall_add_a,0:max_limiters_a),
c-----for density fall near the wall
     &             rr_add(nxeqd_add_a),zz_add(nyeqd_add_a),
     &distance_to_wall(nxeqd_add_a,nyeqd_add_a,0:max_limiters_a),
     &density_r_z(nxeqd_add_a,nyeqd_add_a,nbulka,0:max_limiters_a),
     &density_r_z_rr(nxeqd_add_a,nyeqd_add_a,nbulka,0:max_limiters_a),
     &density_r_z_zz(nxeqd_add_a,nyeqd_add_a,nbulka,0:max_limiters_a),
     &density_r_z_rrzz(nxeqd_add_a,nyeqd_add_a,nbulka,0:max_limiters_a),
     &temperature_r_z(nxeqd_add_a,nyeqd_add_a,nbulka),  
     &dtemperature_r_z_r(nxeqd_add_a,nyeqd_add_a,nbulka),
     &temperature_r_z_zz(nxeqd_add_a,nyeqd_add_a,nbulka),
     &temperature_r_z_rrzz(nxeqd_add_a,nyeqd_add_a,nbulka),  
     &   n_wall_add(0:max_limiters_a),
c-------pointers for density and zeff at rz mesh 
     &dens_rz_in,      
     &zeff_rz_in
c
c  nxeqda must be >= nxeqd
c  nyeqda must be >= nyeqd
c
c----------------------------------------------------------------
c thetapol_wall,       ! poloidal angles [radians] of
c thetapol_limiter     ! wall and limiter points     
c rho_wall,rho_limiter !small radius
c
