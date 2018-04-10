c----------------------------------------------------------------
c data for density profile ouside LCFS which are not in the input
c namelist /edge_prof_nml/
c----------------------------------------------------------------- 
      real*8
     &theta_pol_edge_dens_ar_radian, !polidal angle mesh at radian 
     &sigmedgn_deriv,!second derivatives from sigmedgn by poloidal angle
                     !for the spline approximation
     &theta_pol_edge_1_radian,  !for analitical formula of sigma_edge_n
     &theta_pol_edge_2_radian,             
     &sigma_theta_pol_edge_1_radian,
     &sigma_theta_pol_edge_2_radian
   
      common /edge_prof_no_nml/
     &theta_pol_edge_dens_ar_radian(n_pol_edge_dens_a),
     &sigmedgn_deriv(n_pol_edge_dens_a),
     &theta_pol_edge_1_radian,           
     &theta_pol_edge_2_radian,             
     &sigma_theta_pol_edge_1_radian,
     &sigma_theta_pol_edge_2_radian 
    
