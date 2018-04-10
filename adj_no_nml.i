      integer 
     & iout3,iout4,iout5,iout6,iout7,iout8

c      real*8
c     & psis,
c     & cxc,
c     & czc,
c     & bbar, 
c     & th0a,ua,
c     & dene_chi,teme_chi,zi_chi,     
c     & chi_3d, 
c     & bmin_chi,bmax_chi,umax_chi,
c     & th0max_chi,th0a_chi,
c----second derivatibes
c     &chi_tt,
c     &chi_uu, 
c-----fourth derivatives                         
c     &chi_ttuu,
c------------------------------
c     &thetae_1,
c     &ba_2d,
c     &d_chi_d_theta,d_chi_d_u,
c--------------------------------
c     &Q_safety_adj,sigma_E_adj,dens_averaged_ar,rho_adj

      real*8, pointer ::
     & psis(:),                   !(1:npsi0_a),
     & cxc(:,:),                  !(1:nthp0_a+1,1:npsi0_a),
     & czc(:,:),                  !(1:nthp0_a+1,1:npsi0_a),
     & bbar(:),                   !(1:npsi0_a),
c     & th0a(:),                   !(1:imax_chi_a),     
     & ua(:),                     !(1:npsi0_a),
     & dene_chi(:),               !(1:npsi0_a),
     & teme_chi(:),               !(1:npsi0_a),
     & zi_chi(:),                 !(1:npsi0_a),
     & chi_3d(:,:,:),             !(1:imax_chi_a,1:nmax_chi_a,1:npsi0_a)
     & bmin_chi(:),               !(1:npsi0_a),
     & bmax_chi(:),               !(1:npsi0_a),
c     & umax_chi(:),               !(1:npsi0_a), 
     & chi_tt(:,:,:),             !(1:imax_chi_a,1:nmax_chi_a,1:npsi0_a)
     & chi_uu(:,:,:),             !(1:imax_chi_a,1:nmax_chi_a,1:npsi0_a) 
     & chi_ttuu(:,:,:),           !(1:imax_chi_a,1:nmax_chi_a,1:npsi0_a)
c     & th0max_chi(:),             !(1:npsi0_a),
     & th0a_chi(:,:),             !(1:imax_chi_a,1:npsi0_a),
     & thetae_1(:),               !(0:nthp0_a-1)
     & ba_2d(:,:),                !(0:nthp0_a-1,1:npsi0_a),
     & d_chi_d_theta(:,:,:),      !(1:imax_chi_a,1:nmax_chi_a,1:npsi0_a),
     & d_chi_d_u(:,:,:),          !(1:imax_chi_a,1:nmax_chi_a,1:npsi0_a), 
     & Q_safety_adj(:),           !(1:npsi0_a),
     & sigma_E_adj(:),            !(1:npsi0_a),  
     & dens_averaged_ar(:),       !(1:npsi0_a),
     & rho_adj(:)                 !(1:npsi0_a),

      common /adj_no_nml/ 
     & psis,
     & cxc,
     & czc,
     & bbar,
c     & th0a,         !   pitch angle
     & ua,           !   momentum
     & dene_chi,     !   electron density in 10**13 cm**-3
     & teme_chi,     !   electron temperature in KeV
     & zi_chi,       !   zeff 
     & chi_3d,       !   3D chi function
     & bmin_chi,
     & bmax_chi,
c     & umax_chi,     !   maximal velocity [cm/sec]
c-----second derivativers
     & chi_tt,       !   d^2chi/dpitch
     & chi_uu,       !   d^2chi/dx  x=u0
c-----fourth derivatives                         
     & chi_ttuu,     !   d^4chi/(d/dpitch)**2*(d/dx)**2
c-----
c     & th0max_chi,
     & th0a_chi,
c-------------------------
     & thetae_1,     !   poloidal angle mesh for b line integration 
     & ba_2d,        !   b/b0  along b line
c---------------------------------------------------------------------------
     & d_chi_d_theta,
     & d_chi_d_u,
     & Q_safety_adj, !   safety factor used in adj
     & sigma_E_adj,  !   conductivity used in adj_adj
     & dens_averaged_ar,  !bounce averaged density from
                          !the relativistic maxwellian distribution
     & rho_adj,           !small radii where adj function is calculated 
     & iout3,iout4,iout5,iout6,iout7,iout8
 
