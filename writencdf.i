

c arrays for the ray data for the FP (CQL3D) code
c as a netcdf file
c     nrelta is max value of nrelt
c     nraya  is max value of nray

      real*8, pointer ::
     &ws_nc(:,:),              !(nrelta,nraya),
     &seikon_nc(:,:),          !(nrelta,nraya),
     &spsi_nc(:,:),            !(nrelta,nraya),
     &wr_nc(:,:),              !(nrelta,nraya),
     &wphi_nc(:,:),            !(nrelta,nraya),
     &wz_nc(:,:),              !(nrelta,nraya),
     &wnpar_nc(:,:),           !(nrelta,nraya),
     &wnper_nc(:,:),           !(nrelta,nraya),
     &delpwr_nc(:,:),          !(nrelta,nraya),
     &sdpwr_nc(:,:),           !(nrelta,nraya),
     &wdnpar_nc(:,:),          !(nrelta,nraya),
     &fluxn_nc(:,:),           !(nrelta,nraya),
     &sbtot_nc(:,:),           !(nrelta,nraya),
     &sb_z_nc(:,:),            !(nrelta,nraya),
     &sb_r_nc(:,:),            !(nrelta,nraya),
     &sb_phi_nc(:,:),          !(nrelta,nraya),
     &sene_nc(:,:),            !(nrelta,nraya),
     &ste_nc(:,:),             !(nrelta,nraya),
     &szeff_nc(:,:),           !(nrelta,nraya),
     &salphac_nc(:,:),         !(nrelta,nraya),
     &salphal_nc(:,:),         !(nrelta,nraya),
     &salphas_nc(:,:,:),       !(nrelta,nraya,nbulka),
     &vgr_z_nc(:,:),           !(nrelta,nraya),
     &vgr_r_nc(:,:),           !(nrelta,nraya),
     &vgr_phi_nc(:,:),         !(nrelta,nraya),
     &flux_z_nc(:,:),          !(nrelta,nraya),
     &flux_r_nc(:,:),          !(nrelta,nraya),
     &flux_phi_nc(:,:),        !(nrelta,nraya),
     &wn_r_nc(:,:),            !(nrelta,nraya),
     &wn_z_nc(:,:),            !(nrelta,nraya),
     &wn_phi_nc(:,:),          !(nrelta,nraya),
     &transm_ox_nc(:),               !(nraya),
     &cn_par_optimal_nc(:),          !(nraya),
     &cnpar_ox_nc(:),                !(nraya),
     &cn_b_gradpsi_nc(:),            !(nraya),
c-----emission data--------------------------      
     &wsn_nc(:,:),             !(nrelta,nfreqa),  
     &wcnpar_em_nc(:,:),       !(nrelta,nfreqa),  
     &wcnper_em_nc(:,:),       !(nrelta,nfreqa),  
     &wz_em_nc(:,:),           !(nrelta,nfreqa),  
     &wr_em_nc(:,:),           !(nrelta,nfreqa),  
     &wphi_em_nc(:,:),         !(nrelta,nfreqa),  
     &wal_emis_nc(:,:),        !(nrelta,nfreqa),  
     &wj_emis_nc(:,:),         !(nrelta,nfreqa),            
     &wnray_nc(:,:),           !(nrelta,nfreqa),   
     &win_sn_nc(:,:),          !(nrelta,nfreqa),  
     &win_0_nc(:,:),           !(nrelta,nfreqa),  
     &w_specific_intensity_nc(:,:),             !(nrelta,nfreqa),  
     &wi_0_nc(:,:),                  !(nraya,nfreqa),
     &wi_0sn_nc(:,:),          !(nrelta,nfreqa),
     &wtemp_em_nc(:,:),        !(nrelta,nfreqa),
     &wtemp_rad_em_nc(:,:),    !(nrelta,nfreqa),
     &wtemp_rad_fr_wall_nc(:,:),     !(nraya,nfreqa),
     &waveraged_temp_rad_fr_nc(:),                 !(nfreqa),
     &waveraged_temp_rad_fr_wall_nc(:),            !(nfreqa),
     &wtemp_rad_fr_nc(:,:),          !(nraya,nfreqa),
     &wtemp_pl_fr_nc(:,:),           !(nraya,nfreqa),
     &wr0_em_nc(:,:),                !(nraya,nfreqa),
     &wz0_em_nc(:),                                !nfreqa),
     &wrho0_em_nc(:,:),              !(nraya,nfreqa),
     &wtaun_em_nc(:,:),        !(nrelta,nfreqa),
     &wfreq_nc(:),                                 !(nfreqa),
     &wtau_em_nc(:,:),               !(nraya,nfreqa),
     &wi_0t_nc(:,:),                 !(nraya,nfreqa),
     &wr_2nd_harm_nc(:,:),           !(nraya,nfreqa),
     &wtemp_2nd_harm_nc(:,:),        !(nraya,nfreqa),
     &wr_emis_initial_mesh_nc(:,:),         !(nrelta,nfreqa),
     &wp_perpmax_dmc_nc(:,:,:),    !(nrelta,n_relt_harm1a:n_relt_harm2a,nfreqa),
     &wp_parmin_dmc_nc(:,:,:),     !(nrelta,n_relt_harm1a:n_relt_harm2a,nfreqa)
     &wp_parmax_dmc_nc(:,:,:),     !(nrelta,n_relt_harm1a:n_relt_harm2a,nfreqa)
     &wj_emis_x_nc(:,:,:),                !(nrelta,jx_kin_a,nfreqa),
     &win_sn_x_nc(:,:,:),                 !(nrelta,jx_kin_a,nfreqa),
     &win_0_x_nc(:,:,:),                  !(nrelta,jx_kin_a,nfreqa),
     &wi_0sn_x_nc(:,:,:),                 !(nrelta,jx_kin_a,nfreqa),
     &wi_0_x_nc(:,:),                              !(jx_kin_a,nfreqa),
     &w_specific_intensity_x_nc(:,:,:),   !(nrelta,jx_kin_a,nfreqa),
     &wdye_nc(:,:),                !(nrelta,nfreqa),
c------------------plasma profiles VS r at z=0  ------------------    
     &w_dens_vs_r_nc(:,:),             !(NRA,nbulka),
     &w_temp_vs_r_nc(:,:),             !(NRA,nbulka),
     &w_zeff_vs_r_nc(:),                  !(NRA),
     &w_r_densprof_nc(:),                 !(NRA),
C---------------------------------------------------------------
     &w_eff_nc(:,:),               !(nrelta,nraya),
     &w_theta_pol_nc(:,:)          !(nrelta,nraya),


      complex*16, pointer ::
     &cwexde_nc(:,:),          !(nrelta,nraya)
     &cweyde_nc(:,:),          !(nrelta,nraya)
     &cwezde_nc(:,:),          !(nrelta,nraya)
     &cweps11_nc(:,:),         !(nrelta,nraya)
     &cweps12_nc(:,:),         !(nrelta,nraya)
     &cweps13_nc(:,:),         !(nrelta,nraya)
     &cweps21_nc(:,:),         !(nrelta,nraya)
     &cweps22_nc(:,:),         !(nrelta,nraya)
     &cweps23_nc(:,:),         !(nrelta,nraya)
     &cweps31_nc(:,:),         !(nrelta,nraya)
     &cweps32_nc(:,:),         !(nrelta,nraya)
     &cweps33_nc(:,:)          !(nrelta,nraya)

    
      integer, pointer :: 
     &i_ox_conversion_nc(:),                   !(nraya),
     &nrayelt_nc(:),                           !(nraya),
c-----emission data-------------------------- 
     &nrayelt_emis_nc(:),                      !(nfreqa),
     &nrayelt_emis_initial_mesh_nc(:),          !(nfreqa),
c-------------------------------------------
cSAP100514
     &iray_status_nc(:,:)                          !(nray,nfreq)


      real*8    
     &freqcy,
c-----emission data--------------------------    
     &freqncy0_nc,
c-----data for power absorprion at reflections 
     &w_tot_pow_absorb_at_refl_nc   !scalar

      integer
     &nharm,nrayl,irayl

  
      common/writencdf/
     &cwexde_nc,
     &cweyde_nc,
     &cwezde_nc,
     &cweps11_nc,
     &cweps12_nc,
     &cweps13_nc,
     &cweps21_nc,
     &cweps22_nc,
     &cweps23_nc,
     &cweps31_nc,
     &cweps32_nc,
     &cweps33_nc,
     &ws_nc,
     &seikon_nc,
     &spsi_nc,
     &wr_nc,
     &wphi_nc,
     &wz_nc,
     &wnpar_nc,
     &wnper_nc,
     &delpwr_nc,
     &sdpwr_nc,
     &wdnpar_nc,
     &fluxn_nc,
     &sbtot_nc,
     &sb_z_nc,
     &sb_r_nc,
     &sb_phi_nc,
     &sene_nc,
     &ste_nc,
     &szeff_nc,
     &salphac_nc,
     &salphal_nc,
     &salphas_nc,
     &vgr_z_nc,
     &vgr_r_nc,
     &vgr_phi_nc,
     &flux_z_nc,
     &flux_r_nc,
     &flux_phi_nc,
     &wn_r_nc,
     &wn_z_nc,
     &wn_phi_nc,
     &transm_ox_nc,
     &cn_par_optimal_nc,
     &cnpar_ox_nc,
     &cn_b_gradpsi_nc,
     &i_ox_conversion_nc, !integer 
     &nrayelt_nc,         !integer
     &freqcy,             !scalar
c-----emission data-------------------------- 
     &wsn_nc,  
     &wcnpar_em_nc,
     &wcnper_em_nc,
     &wz_em_nc,
     &wr_em_nc,
     &wphi_em_nc,  
     &wal_emis_nc,
     &wj_emis_nc,            
     &wnray_nc,     
     &win_sn_nc,
     &win_0_nc,
     &w_specific_intensity_nc,
     &wi_0_nc,
     &wi_0sn_nc,
     &wtemp_em_nc,
     &wtemp_rad_em_nc,
     &wtemp_rad_fr_wall_nc,
     &waveraged_temp_rad_fr_nc,
     &waveraged_temp_rad_fr_wall_nc,
     &wtemp_rad_fr_nc,
     &wtemp_pl_fr_nc,
     &wr0_em_nc,
     &wz0_em_nc,
     &wrho0_em_nc,
     &wtaun_em_nc,
     &wfreq_nc,
     &wtau_em_nc,
     &wi_0t_nc,
     &wr_2nd_harm_nc,
     &wtemp_2nd_harm_nc,
     &freqncy0_nc, !scalar
     &nrayelt_emis_nc        ,
     &nrayelt_emis_initial_mesh_nc,
     &wr_emis_initial_mesh_nc,
     &wp_perpmax_dmc_nc,
     &wp_parmin_dmc_nc,
     &wp_parmax_dmc_nc,
     &wj_emis_x_nc,
     &win_sn_x_nc, 
     &win_0_x_nc,
     &wi_0sn_x_nc,
     &wi_0_x_nc,
     &w_specific_intensity_x_nc,
     &wdye_nc,
c------------------plasma profiles VS r at z=0  ------------------    
     &w_dens_vs_r_nc,
     &w_temp_vs_r_nc,
     &w_zeff_vs_r_nc,
     &w_r_densprof_nc,
c--------------------------------------------------
     &w_eff_nc,
     &w_theta_pol_nc,
c-----------------------------------------------------   
     &w_tot_pow_absorb_at_refl_nc,  !scalar
c-----------------------------------------------------
     &nharm,nrayl,irayl,
c------------------------------------------------------
cSAP100514
     &iray_status_nc  

c the data for mnemonic.nc file (output data in netcdf format)
c nrelta is the max number of the ray elements along any ray

c     nrayelt_nc(iray) is the number of 'iray' ray elements   

c     the data for the top of mnemonic.nc file:
c     freqcy,nharm,nrayl 
c 
c     vgr_z_nc(nrelta,nraya),vgr_r_nc(nrelta,nraya),vgr_phi_nc(nrelta,nraya)
c     are the group-velocity normalized to c (light speed)       
c
c     wn_z_nc(nrelta,nraya),wn_r_nc(nrelta,nraya),wn_phi_nc(nrelta,nraya)
c     are refructive index coordinates
c
c     The data for OX transmission.They are used for i_ox=2 case
c     transm_ox_nc(nraya),cn_par_optimal_nc(nraya),cnpar_ox_nc(nraya),
c      cn_b_gradpsi(nraya),i_ox_conversion_nc(nraya)
c

c     These data are to plot the resonance ellipse boundary
c     for several harmonics   
c     nrayelt_emis_initial_mesh_nc(nfreqa)
c     wr_emis_initial_mesh_em_nc(nrelta,nfreqa),
c     wp_perpmax_dmc_nc(nrelta,n_relt_harm1a:n_relt_harm2a,nfreqa),
c     wp_parmin_dmc_nc(nrelta,n_relt_harm1a:n_relt_harm2a,nfreqa),
c     wp_parmax_dmc_nc(nrelta,n_relt_harm1a:n_relt_harm2a,nfreqa)
    
c     waveraged_temp_rad_fr_nc(nfreqa), temperature averaged over rays
c     waveraged_temp_rad_fr_wall_nc(nfreqa),
cSAP080422
c     wj_emis_x_nc(nrelta,jx_kin_a,nfreqa),
ccccc wj_emis_x_nc(nrelta,jxa,nfreqa)  energy sdpectrum of j_emis at ray points
c                                      at each frequency for cenral ray
c     w_specific_intensity_nc(nrelta,jxa,nfreqa) energy spectrum

c     wdye_nc(nrelta,nfreqa) omega/omega_ce along the ray
c------------------plasma profiles vs r at z=0  ------------------    
c     +w_dens_vs_r_nc(nra,nbulka),w_temp_vs_r_nc(nra,nbulka),
c     +w_zeff_vs_r_nc(nra),w_r_densprof_nc(nra)
c--------------------------------------------------
c     +w_eff_nc(nrelta,nraya)  CD efficiency along all rays    
c     +w_theta_pol_nc((nrelta,nraya), poloidal angle [degree]
c------------------------------------------------------------------------
c w_tot_pow_absorb_at_refl_nc total power [MWatt] absorbed at all reflections
c                        at all rays
c------------------------------------------------------
cSAP100514
c     iray_status_nc(nray,nfreq)  array giving a status code for each type 
c                              of stopping of each ray
c------------------------------------------------------------------
