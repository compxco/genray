


c arrays for the ray data for the FP (CQL3D) code
c as a netcdf file
c     nrelta is max value of nrelt
c     nraya  is max value of nray
      integer nrelta1
      parameter(nrelta1=1)

      real*8
     +ws_1_nc,seikon_1_nc,
     +spsi_1_nc,
     +wr_1_nc,wphi_1_nc,
     +wz_1_nc,wnpar_1_nc,
     +wnper_1_nc,delpwr_1_nc,
     +sdpwr_1_nc,wdnpar_1_nc,
     +fluxn_1_nc,sbtot_1_nc,
     +sb_z_1_nc,sb_r_1_nc,
     +sb_phi_1_nc,
     +sene_1_nc,ste_1_nc,
     +salphac_1_nc,
     +salphal_1_nc,
     +salphas_1_nc,
     +vgr_z_1_nc,vgr_r_1_nc,
     +vgr_phi_1_nc,
     +flux_z_1_nc,flux_r_1_nc,
     +flux_phi_1_nc,
     +wdye_nc_1_nc

      double complex cwexde_1_nc,cweyde_1_nc,cwezde_1_nc
c     &cweps11_1_nc,cweps12_1_nc,cweps13_1_nc,
c     &cweps21_1_nc,cweps22_1_nc,cweps23_1_nc,
c     &cweps31_1_nc,cweps32_1_nc,cweps33_1_nc


      real*8
     +wn_r_1_nc,wn_z_1_nc,
     +wn_phi_1_nc
c     +transm_ox_nc,cn_par_optimal_nc,cnpar_ox_nc,
c     +cn_b_gradpsi_nc,
c     +freqcy

      integer
c     +i_ox_conversion_nc,
     +nrayelt_1_nc
c     +nharm,nrayl,irayl


      real*8
c     +wsn_nc,
c     +wcnpar_em_nc,
c     +wcnper_em_nc,
c     +wz_em_nc,
c     +wr_em_nc,
c     +wphi_em_nc,  
c     +wal_emis_nc,
c     +wj_emis_nc,            
c     +wnray_nc,     
c     +win_sn_nc,     
c     +win_0_nc,
c     +w_specific_intensity_nc,
c     +wi_0_nc,
c     +wi_0sn_nc,
c     +wtemp_em_nc,
c     +wtemp_rad_em_nc,
c     +wtemp_rad_fr_wall_nc,
c     +wtemp_rad_fr_nc,
c     +waveraged_temp_rad_fr_nc,
c     +waveraged_temp_rad_fr_wall_nc,
c     +wtemp_pl_fr_nc,
c     +wr0_em_nc,
c     +wz0_em_nc,
c     +wrho0_em_nc,
c     +wtaun_em_nc,
c     +wfreq_nc,
c     +wtau_em_nc,
c     +wi_0t_nc,
c     +wr_2nd_harm_nc,
c     +wtemp_2nd_harm_nc,
c     +freqncy0_nc,    
c     +wr_emis_initial_mesh_nc,
c     +wp_perpmax_dmc_nc,
c     +wp_parmin_dmc_nc,wp_parmax_dmc_nc,
c     +wj_emis_x_nc,win_sn_x_nc,win_0_x_nc,
c     +wi_0sn_x_nc,wi_0_x_nc,
c     +w_specific_intensity_x_nc,
c     +w_dens_vs_r_nc,w_temp_vs_r_nc,w_zeff_vs_r_nc,w_r_densprof_nc,
     +w_eff_1_nc,
     +w_theta_pol_1_nc
c-----data for power absorprion at reflections 
c     &w_tot_pow_absorb_at_refl_nc

c      integer nrayelt_emis_nc,nrayelt_emis_initial_mesh_nc

      common/writencdf/
     +cwexde_1_nc(nrelta1,nraya),cweyde_1_nc(nrelta1,nraya),
     +cwezde_1_nc(nrelta1,nraya),
     +ws_1_nc(nrelta1,nraya),seikon_1_nc(nrelta1,nraya),
     +spsi_1_nc(nrelta1,nraya),
     +wr_1_nc(nrelta1,nraya),wphi_1_nc(nrelta1,nraya),
     +wz_1_nc(nrelta1,nraya),wnpar_1_nc(nrelta1,nraya),
     +wnper_1_nc(nrelta1,nraya),delpwr_1_nc(nrelta1,nraya),
     +sdpwr_1_nc(nrelta1,nraya),wdnpar_1_nc(nrelta1,nraya),
     +fluxn_1_nc(nrelta1,nraya),sbtot_1_nc(nrelta1,nraya),
     +sb_z_1_nc(nrelta1,nraya),sb_r_1_nc(nrelta1,nraya),
     +sb_phi_1_nc(nrelta1,nraya),
     +sene_1_nc(nrelta1,nraya),ste_1_nc(nrelta1,nraya),
     +salphac_1_nc(nrelta1,nraya),
     +salphal_1_nc(nrelta1,nraya),
     +salphas_1_nc(nrelta1,nraya,nbulka),
     +vgr_z_1_nc(nrelta1,nraya),vgr_r_1_nc(nrelta1,nraya),
     +vgr_phi_1_nc(nrelta1,nraya),
     +flux_z_1_nc(nrelta1,nraya),flux_r_1_nc(nrelta1,nraya),
     +flux_phi_1_nc(nrelta1,nraya),
c     +cweps11_nc(nrelta1,nraya),
c     +cweps12_nc(nrelta1,nraya),
c     +cweps13_nc(nrelta1,nraya),
c     +cweps21_nc(nrelta1,nraya),
c     +cweps22_nc(nrelta1,nraya),
c     +cweps23_nc(nrelta1,nraya),
c     +cweps31_nc(nrelta1,nraya),
c     +cweps32_nc(nrelta1,nraya),
c     +cweps33_nc(nrelta1,nraya),
     +wn_r_1_nc(nrelta1,nraya),wn_z_1_nc(nrelta1,nraya),
     +wn_phi_1_nc(nrelta1,nraya),
c     +transm_ox_nc(nraya),cn_par_optimal_nc(nraya),cnpar_ox_nc(nraya),
c     +cn_b_gradpsi_nc(nraya),
c     +i_ox_conversion_nc(nraya),
     +nrayelt_1_nc(nraya),
c     +freqcy,
c     +nharm,nrayl,irayl,
c-----emission data-------------------------- 
c     +wsn_nc(nrelta1,nfreqa),  
c     +wcnpar_em_nc(nrelta1,nfreqa),
c     +wcnper_em_nc(nrelta1,nfreqa),
c     +wz_em_nc(nrelta1,nfreqa),
c     +wr_em_nc(nrelta1,nfreqa),
c     +wphi_em_nc(nrelta1,nfreqa),  
c     +wal_emis_nc(nrelta1,nfreqa),
c     +wj_emis_nc(nrelta1,nfreqa),            
c     +wnray_nc(nrelta1,nfreqa),     
c     +win_sn_nc(nrelta1,nfreqa),
c     +win_0_nc(nrelta1,nfreqa),
c     +w_specific_intensity_nc(nrelta1,nfreqa),
c     +wi_0_nc(nraya,nfreqa),
c     +wi_0sn_nc(nrelta1,nfreqa),
c     +wtemp_em_nc(nrelta1,nfreqa),
c     +wtemp_rad_em_nc(nrelta1,nfreqa),
c     +wtemp_rad_fr_wall_nc(nraya,nfreqa),
c     +waveraged_temp_rad_fr_nc(nfreqa),
c     +waveraged_temp_rad_fr_wall_nc(nfreqa),
c     +wtemp_rad_fr_nc(nraya,nfreqa),
c     +wtemp_pl_fr_nc(nraya,nfreqa),
c     +wr0_em_nc(nraya,nfreqa),
c     +wz0_em_nc(nfreqa),
c     +wrho0_em_nc(nraya,nfreqa),
c     +wtaun_em_nc(nrelta1,nfreqa),
c     +wfreq_nc(nfreqa),
c     +wtau_em_nc(nraya,nfreqa),
c     +wi_0t_nc(nraya,nfreqa),
c     +wr_2nd_harm_nc(nraya,nfreqa),
c     +wtemp_2nd_harm_nc(nraya,nfreqa),
c     +freqncy0_nc,
c     +nrayelt_emis_nc(nfreqa),
c     +nrayelt_emis_initial_mesh_nc(nfreqa),
c     +wr_emis_initial_mesh_nc(nrelta1,nfreqa),
c     +wp_perpmax_dmc_nc(nrelta1,n_relt_harm1a:n_relt_harm2a,nfreqa),
c     +wp_parmin_dmc_nc(nrelta1,n_relt_harm1a:n_relt_harm2a,nfreqa),
c     +wp_parmax_dmc_nc(nrelta1,n_relt_harm1a:n_relt_harm2a,nfreqa),
ccSAP080422
cc     +wj_emis_x_nc(nrelta1,jxa,nfreqa), 
cc     +win_sn_x_nc(nrelta1,jxa,nfreqa), 
cc     +win_0_x_nc(nrelta1,jxa,nfreqa),
cc     +wi_0sn_x_nc(nrelta1,jxa,nfreqa),
cc     +wi_0_x_nc(jxa,nfreqa),
cc     +w_specific_intensity_x_nc(nrelta1,jxa,nfreqa),
cc    +wj_emis_x_nc(nrelta1,jx_kin_a,nfreqa),
c     +win_sn_x_nc(nrelta1,jx_kin_a,nfreqa), 
c     +win_0_x_nc(nrelta1,jx_kin_a,nfreqa),
c     +wi_0sn_x_nc(nrelta1,jx_kin_a,nfreqa),
c     +wi_0_x_nc(jx_kin_a,nfreqa),
c     +w_specific_intensity_x_nc(nrelta1,jx_kin_a,nfreqa),
c     +wdye_nc(nrelta1,nfreqa),
c------------------plasma profiles VS r at z=0  ------------------    
c     +w_dens_vs_r_nc(NRA,nbulka),w_temp_vs_r_nc(NRA,nbulka),
c     +w_zeff_vs_r_nc(NRA),w_r_densprof_nc(NRA),
c--------------------------------------------------
     +w_eff_1_nc(nrelta1,nraya),
     +w_theta_pol_1_nc(nrelta1,nraya)
c-----------------------------------------------------   
c     &w_tot_pow_absorb_at_refl_nc,
c-----------------------------------------------------
c     &nharm,nrayl,irayl
c------------------------------------------------------
     

c the data for mnemonic.nc file (output data in netcdf format)
c nrelta is the max number of the ray elements along any ray

c     nrayelt_nc(iray) is the number of 'iray' ray elements   

c     the data for the top of mnemonic.nc file:
c     freqcy,nharm,nrayl 
c 
c     vgr_z_nc(nrelta1,nraya),vgr_r_nc(nrelta1,nraya),vgr_phi_nc(nrelta1,nraya)
c     are the group-velocity normalized to c (light speed)       
c
c     wn_z_nc(nrelta1,nraya),wn_r_nc(nrelta1,nraya),wn_phi_nc(nrelta1,nraya)
c     are refructive index coordinates
c
c     The data for OX transmission.They are used for i_ox=2 case
c     transm_ox_nc(nraya),cn_par_optimal_nc(nraya),cnpar_ox_nc(nraya),
c      cn_b_gradpsi(nraya),i_ox_conversion_nc(nraya)
c

c     These data are to plot the resonance ellipse boundary
c     for several harmonics   
c     nrayelt_emis_initial_mesh_nc(nfreqa)
c     wr_emis_initial_mesh_em_nc(nrelta1,nfreqa),
c     wp_perpmax_dmc_nc(nrelta1,n_relt_harm1a:n_relt_harm2a,nfreqa),
c     wp_parmin_dmc_nc(nrelta1,n_relt_harm1a:n_relt_harm2a,nfreqa),
c     wp_parmax_dmc_nc(nrelta1,n_relt_harm1a:n_relt_harm2a,nfreqa)
    
c     waveraged_temp_rad_fr_nc(nfreqa), temperature averaged over rays
c     waveraged_temp_rad_fr_wall_nc(nfreqa),
cSAP080422
c     wj_emis_x_nc(nrelta1,jx_kin_a,nfreqa),
ccccc wj_emis_x_nc(nrelta1,jxa,nfreqa)  energy sdpectrum of j_emis at ray points
c                                      at each frequency for cenral ray
c     w_specific_intensity_nc(nrelta1,jxa,nfreqa) energy spectrum

c     wdye_nc(nrelta1,nfreqa) omega/omega_ce along the ray
c------------------plasma profiles vs r at z=0  ------------------    
c     +w_dens_vs_r_nc(nra,nbulka),w_temp_vs_r_nc(nra,nbulka),
c     +w_zeff_vs_r_nc(nra),w_r_densprof_nc(nra)
c--------------------------------------------------
c     +w_eff_nc(nrelta1,nraya)  CD efficiency along all rays    
c     +w_theta_pol_nc((nrelta1,nraya), poloidal angle [degree]
c------------------------------------------------------------------------
c w_tot_pow_absorb_at_refl_nc total power [MWatt] absorbed at all reflections
c                        at all rays
c------------------------------------------------------
