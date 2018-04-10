c------------------------------------------------------------------
c namelist to set the table of exponential density fall off distance
c outside LCFS starting at rho=1 density 
c at different polidal mesh  
c--------------------------------------------------------------------
      namelist /edge_prof_nml/ 
     &i_edge_dens_anal,
     &i_edge_dens_rz_mesh,
     &n_pol_edge_dens,
     &sigmedgn,sigmedgt,
     &theta_pol_edge_dens_ar_degree,
     &sigmedgn_ar,
     &dens_min_edge,temp_min_edge,  
     &theta_pol_edge_1_degree,  !for analitical formula of sigma_edge_n
     &theta_pol_edge_2_degree,            ! 
     &sigma_theta_pol_edge_1_degree,
     &sigma_theta_pol_edge_2_degree,
     &sigma_edgen_0,
     &sigma_edgen_1,
     &sigma_edgen_2,
     &sigma_wall_n,
     &sigma_lim_toroidal_degree, 
     &nxeqd_add,nyeqd_add
