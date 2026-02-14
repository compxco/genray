c      all input namelists

      include 'name_genr.i'               !/namelist /genr/
     
      include 'name_tokamak.i'            !/namelist/ tokamak/     

      namelist /wave/ frqncy,ioxm,ireflm,jwave,istart,delpwrmn,ibw,
     *i_vgr_ini,poldist_mx,ioxm_n_npar,
     &cnperp_plot_min,cnperp_plot_max,n_nperp_plot,
     &cN_perp_root_max,n_points_root,
     &i_look_roots,k_hot_root,
     &i_rho_find_hot_nperp_roots,
     &rho_step_find_hot_nperp_roots,rho_min_find_hot_nperp_roots, 
     &shift_rho_geom,
     &no_reflection,
     & istep_in_lcfs !YuP[2020-09-03] added

      namelist /scatnper/ iscat,scatd,rhoscat,
cSAP120511
     &iscat_lh_nicola
c----------------------

      namelist /dispers/ ib,id,iherm,iabsorp,iswitch,del_y,jy_d,
     *idswitch,iabswitch,n_relt_harm,n_relt_intgr,iflux,
     &i_im_nperp,i_geom_optic,ray_direction,errabs0,errrel0,navg,
     &diff_err,relres,iabsorp_collisional,coll_mult,refl_loss,
     &n_relt_harm1,i_salphal,ion_absorption,
     &iabsorp_ql

      namelist /numercl/ irkmeth,ndim1,isolv,idif,nrelt,
     * prmt1,prmt2,prmt3,prmt4,prmt6,prmt9,icorrect,iout3d,
     * maxsteps_rk,i_output,
     & nharm_refined_step, ![2025-12-16]
     & i_uh_switch,uh_switch,prmt6_uh_switch,
     &toll_hamilt,  
     + dL_step, dN_step, der_r, der_n, der_f,
     &i_power_switch_resonance, 
     &prmt6_power_switch_resonance,
     &n_power_switch_resonance,
     &y_power_switch_resonance,
     &del_y_power_switch_resonance,
     &i_resonance_curve_integration_method,epsi,  
     &eps_delta_pow

      namelist /output/ iwcntr,iwopen,iwj,itools,i_plot_b,i_plot_d,
     &n_plot_disp,r_plot_disp,id_plot_disp,
     &z_plot_disp,n_parallel_plot_disp,    
     &max_r_nperp_plot_disp,
     &min_r_nperp_plot_disp,
     &max_i_nperp_plot_disp,
     &min_i_nperp_plot_disp,
     &s_poloid_plot_disp,point_plot_disp,
     &i_plot_disp_cold,
     &n_plot_disp_cold,s_poloid_plot_disp_cold,r_plot_disp_cold,
     &point_plot_disp_cold,     
     &i_plot_wave_normal_cold,
     &number_map_points_real_nperp,number_map_points_image_nperp,
     &ratio_min_r_nperp,ratio_max_r_nperp,
     &ratio_min_i_nperp,ratio_max_i_nperp, 
     &n_contour_plot_disp,
     &r_freq,z_freq,alpha_freq,beta_freq,dist_freq,max_plot_freq,
     &nsteps_freq,n_ec_harmonics_freq,npar_freq

      namelist /plasma/ nbulk,izeff,idens,model_rho_dens,
     +  temp_scale,den_scale,
!     +  elx0,ely0,elz0, elax,elay,elaz,  ! YuP added !
!     +  dens0rr,dens0es,dens0ub,eltheta,sintt,costt, ! YuP added !
!     +  Rm0rr,Rm0es,rtau,r_ub_edge,      ! YuP added !
     +  ndens,
     &  nonuniform_profile_mesh,
     &  dendsk, tempdsk, tpopdsk,  ! YuP[2024-08-14] added !
     &  dens_read, temp_read, tpop_read !YuP[2024-08-14] added

      namelist /species/ charge,dmas

      namelist /varden/ var0,denn,denm,an,sigman

      namelist /denprof/ dense0,denseb,rn1de,rn2de   

      namelist /tpopprof/ tp0,tpb,rn1tp,rn2tp

      namelist /vflprof/vfl0,vflb,rn1vfl,rn2vfl

      namelist /tprof/ ate0,ateb,rn1te,rn2te

      namelist /zprof/ zeff0,zeffb,rn1zeff,rn2zeff

c----------------------------------------------------------
c     namelists for plasma profiles at uniform radial mesh:
c----------------------------------------------------------
      include 'name_uniform_mesh_profiles.i'     
c      namelist /dentab/    prof  
c      namelist /temtab/    prof  
c      namelist /tpoptab/   prof  
c      namelist /vflowtab/  prof  
c      namelist /zeftab/ zeff1
c      namelist /EdcTSCtab/ prof
c      namelist /JparTSCtab/prof
c---------------------------------------------------------
c     namelists for plasma profiles at non-uniform radial mesh
c     written by lines
c------------------------------------------------------------
      include 'name_non_uniform_mesh_profiles_line.i'  
c      namelist /dentab_nonuniform_line/   nj_tab,prof_2d,radii_2d
c      namelist /temtab_nonuniform_line/   nj_tab,prof_2d,radii_2d
c      namelist /tpoptab_nonuniform_line/  nj_tab,prof_2d,radii_2d
c      namelist /vflowtab_nonuniform_line / nj_tab,prof_2d,radii_2d
c      namelist /zeftab_nonuniform_line /   nj_tab,prof_2d,radii_2d
c-------------------------------------------------------------------

      namelist /read_diskf/ i_diskf,
     . netcdfnm,
     . rtem0,
     . rtail,r1t,r2t,ttail,rhot,r1h,r2h,thotpar,thotper,
     . hotexp,hotmnpar,hotmxpar,hotmnper,hotmxper,
     . rbeam,r1b,r2b,tbeam,ebeam,thbeam,
     . jx,iym,lrz,ngen,
CENM 31Aug05 Added (optional) parameters at the end if dispers namelist
C    to be used in the relativistic dispersion relation in abhay_disp.f
     & rvtovte1,rvtovte2,rtemp1,rtemp2,rtemp3
      namelist /emission/i_emission,tol_emis,nharm1,nharm2,nfreq,
     . freq00,freq01,wallr,i_rrind,i_r_2nd_harm,
     & i_emission_spectrum,jx_kin,max_kin_energy_kev

      include 'name_grill.i'              !namelist /grill/
     
      namelist /ox/ i_ox,
     &theta_bot,theta_top,i_ox_poloidal_max,eps_antenna,
     &eps_xe


      include 'name_eccone.i'              !namelist /eccone/
     
      include 'name_adj.i'                 !namelist /adj_nml/     
cSAP090203
      include 'name_edge_prof_nml.i'       !namelist /edge_prof_nml/

cSAP091221
      include 'name_lsc_approach_nml.i'    !namelist /lsc_approach_nml/
