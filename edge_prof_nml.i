c---------------------------------------------------------------------------
c data to set density profile ouside LCFS
c---------------------------------------------------------------------------
      integer
     &i_edge_dens_anal,  ! =0
                         ! =1 analytical formula for sigmedgn(theta_pol)
                         ! =2 table data for  sigmedgn(theta_pol)
     &i_edge_dens_rz_mesh, ! =0 do not use density_rz 
                           ! =1 to read density_rz from namelist /dens_rz/ 
     &n_pol_edge_dens,   ! number of poloidal points of 
                         ! the poloidal mesh
     &nxeqd_add,nyeqd_add ! number of poits ar RZ mesh used for creation 
                         ! density_r_z  
      real*8
cSAP090304
     &sigmedgn,          !normalized exponential density and temperature fall 
     &sigmedgt,          !off dist utside LCFS starting at rho=1
     &theta_pol_edge_dens_ar_degree, !poloidal angle mesh
                                     !to set
                                     !sigmedgn_ar
     &sigmedgn_ar,                   !exponential density fall off distance
                                     !outside LCFS starting at rho=1 density 
     &dens_min_edge,                 !minimal edge density 
                                     !10**13/cm**3 for all plasma species
     &temp_min_edge,                 !minimal edge temperature
                                     ![KeV] for all plasma species

     &theta_pol_edge_1_degree,       !for analitical formula of sigma_edge_n
     &theta_pol_edge_2_degree,       ! 
     &sigma_theta_pol_edge_1_degree,
     &sigma_theta_pol_edge_2_degree,   
     &sigma_edgen_0,   
     &sigma_edgen_1,
     &sigma_edgen_2,
     &sigma_wall_n,
     &sigma_lim_toroidal_degree 

      common /edge_prof_nml/
c---- real*8  
c-----from namelist /name_edge_prof_nml/ 
     &sigmedgn,          
     &sigmedgt,   
     &theta_pol_edge_dens_ar_degree(1:n_pol_edge_dens_a),
     &sigmedgn_ar(1:n_pol_edge_dens_a),
     &dens_min_edge,temp_min_edge,
     &theta_pol_edge_1_degree,  !for analitical formula of sigma_edge_n
     &theta_pol_edge_2_degree,             
     &sigma_theta_pol_edge_1_degree,
     &sigma_theta_pol_edge_2_degree,
     &sigma_edgen_0,     
     &sigma_edgen_1,
     &sigma_edgen_2,
     &sigma_wall_n,
     &sigma_lim_toroidal_degree, 
c-----integer
     &i_edge_dens_anal, 
     &i_edge_dens_rz_mesh,
     &n_pol_edge_dens,
     &nxeqd_add,nyeqd_add  

