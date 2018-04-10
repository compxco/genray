c data for LSC approach
c      velocity mesh is normalized at v_e_thermal different
c     at different radial points (1:n_psi_TSC)
     
      real*8, pointer :: 
     &rho_bin_lsc(:),                    !(1:n_psi_TSC+1)
     &rho_bin_center_lsc(:),             !(1:n_psi_TSC)
     &binvol_lsc(:),                     !(1:n_psi_TSC)
     &binarea_lsc(:),                    !(1:n_psi_TSC
     &v_te_bin_center_lsc(:),            !(1:n_psi_TSC)
     &tau_n_bin_center_lsc(:),           !(1:n_psi_TSC) 
     &del_s_pol_bin_lsc(:),              !(1:n_psi_TSC)  
     &v_par_mesh_lsc(:),                 !(-nv_lsc_ql_lh : nv_lsc_ql_lh)
     &d_ql_lsc_lh_ar(:,:),               !(1:n_psi_TSC,-nv_lsc_ql_lh : nv_lsc_ql_lh)
     &integral_fe_lsc(:,:),              !(1:n_psi_TSC,-nv_lsc_ql_lh : nv_lsc_ql_lh)
     &fe_lsc(:,:),                       !(1:n_psi_TSC,-nv_lsc_ql_lh : nv_lsc_ql_lh)
     &d_fe_dv_lsc(:,:),                  !(1:n_psi_TSC,-nv_lsc_ql_lh : nv_lsc_ql_lh)
     &d_ql_lsc_lh_ar_old(:,:),           !(1:n_psi_TSC,-nv_lsc_ql_lh : nv_lsc_ql_lh)
     &fe_lsc_old(:,:),                   !(1:n_psi_TSC,-nv_lsc_ql_lh : nv_lsc_ql_lh)
     &delpwr_nc_old_lsc(:,:),            !(nrelta,nraya)   
     &j_rf_TSC_1D(:),                    !(1:n_psi_TSC)
     &d_j_rf_d_Edc_TSC_1D(:),            !(1:n_psi_TSC)
     &d_ln_j_rf_d_ln_Edc_TSC_1D(:),      !(1:n_psi_TSC)
     &power_dens_watt_m3_TSC_1D(:),      !(1:n_psi_TSC)
     &CD_dens_no_Edc_a_m2_TSC_1D(:),     !(1:n_psi_TSC)
     &CD_dens_small_Edc_a_m2_TSC_1D(:),  !(1:n_psi_TSC)
     &d_ql_lsc_lh_ar_loc(:)              !(-nv_lsc_ql_lh : nv_lsc_ql_lh)

      common /lsc_approach_no_nml/
     &rho_bin_lsc,
     &rho_bin_center_lsc,
     &binvol_lsc,                        ![cm**2]
     &binarea_lsc,                       ![cm**3]
     &v_te_bin_center_lsc,               ![cm/sec]
     &tau_n_bin_center_lsc,              ![sec]
     &del_s_pol_bin_lsc,     
     &v_par_mesh_lsc,
     &d_ql_lsc_lh_ar,     
     &integral_fe_lsc,
     &fe_lsc,   
     &d_fe_dv_lsc,  
     &d_ql_lsc_lh_ar_old,
     &fe_lsc_old,               
     &delpwr_nc_old_lsc,   
     &j_rf_TSC_1D,                     
     &d_j_rf_d_Edc_TSC_1D,              
     &d_ln_j_rf_d_ln_Edc_TSC_1D,
     &power_dens_watt_m3_TSC_1D,          ![watt/m**3]
     &CD_dens_no_Edc_a_m2_TSC_1D,         !A/m**2  CD density without E_DC
     &CD_dens_small_Edc_a_m2_TSC_1D,      !A/m**2  CD density with small E_DC
     &d_ql_lsc_lh_ar_loc        
