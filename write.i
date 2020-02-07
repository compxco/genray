
c the arrays for the text and netCDF ray data for FP (CQL3D) code
c      nrelta is maximum value of nrelt
c      in all array nfreqa -> nfreq
       complex*16 cwexde,cweyde,cwezde !YuP[2020-01] double complex to complex*16
       complex*16 w_ceps  !YuP[2020-01] double complex to complex*16
cBH130508:  Variable nfreq and nray made it difficult to set up
cBH130508:  buffers for MPI send and receive.  Morevover, only
cBH130508:  a small amount of memory is saved here on this account.
cBH130508:  Therefore, reverting back to parameters nfreqa and nraya, 
cBH130508:  for the limited amount of nray and nfreq dimensioned ray 
cBH130508:  data below.

       real*8
     1              ws,seikon,spsi,
     2              wr,wphi,wz,wnpar,
     3              wnper,delpwr,sdpwr,
     4              wdnpar,fluxn,sbtot,
     4		    sb_z,sb_r,sb_phi,
     5              sene,ste,szeff,salphac,
     5              salphal,
     5              wye,wyi,wyi2,
     5              wxi,wxi2,
     6              xarr,yarr,
     7              rez,
     9              eff,wmtor,
     &              wn_z,wn_r,wn_phi,
     &              wvgr_z,wvgr_r,
     &              wvgr_phi, 
cSAP080906
     &              wtheta_pol,
c    tokman flux
     &              wflux_tokman_z, !z,r,phi components   
     &              wflux_tokman_r,wflux_tokman_phi,
     &              salphas, 
     &              wt,
c------------------------------------------------------
c    the data for the emission calculations           
     +              wcnz_em,wcnr_em,wcm_em,
     +              wz_em,wr_em,wphi_em,  
     +              wal_emis,wj_emis,            
     +              wnray,wsn,     
     +              win_sn,win_0,
c     +              wi_0,
     &wi_0sn,
     +              wtemp_em,wtemp_rad_em,
c     +              wtemp_rad_fr_wall,
c     +              wtemp_rad_fr,
c     +              wtemp_pl_fr,
c     +              wr0_em,wz0_em,wrho0_em,
     +              wtaun_em,
c     &wfreq,
c     +              wtau_em,
c     +              wi_0t,
c     +              wr_2nd_harm,
c     +              wtemp_2nd_harm,
     +              freqncy0,
c     +              wn_perp_ioxm_p,wn_perp_ioxm_m,
c     +              wn_perp_ioxm_n_npar_p,wn_perp_ioxm_n_npar_m,
c     +              wye_0,wxe_0,
c-------------------------------------------------------
c    for emission at O-X ebw jump case
c-------------------for O mode ray part 
     +              wsn_o,wz_o,wr_o,wphi_o,
     +              wcnz_em_o,wcnr_em_o,wcm_em_o,
c-------------------for X_EBW mode ray part 
     +              wsn_x,wz_x,wr_x,wphi_x,
     +              wcnz_em_x,wcnr_em_x,wcm_em_x,
c-------------------------------------------------------
c   the data for the determination of the N_par boundaries
     9		    wnrho,gnpar,wxe,wp,
     9		    wnparplb,wnparmnb,
c-------------------------------------------------------
     1              phiold,
     8              zold,rold,rhoold,cld,cvac,
c--------------------------------------------------------
c     These data are for OX transmission coefficient       
     &              transm_ox,cn_par_optimal,cnpar_ox,cn_b_gradpsi,
c-------------------------------------------------------
c     These data are to plot the resonance ellipse boundary
c     for several harmonics    
     &              wp_perpmax_dmvt,wp_parmin_dmvt,wp_parmax_dmvt,
     &              wp_perpmax_dmc,wp_parmin_dmc,wp_parmax_dmc,
cSm060303
     &              wp_0_dmvt,wdel_p_dmvt,
c-------------------------------------------------------------------
c     data for power absorprion at reflections
     &w_tot_pow_absorb_at_refl,
cSAP080902
c-----data to plot delta power_e and delta current along the ray
     &             delpow_e_ar,delcur_par_ar 

      integer       w_ires  

      integer       nrayelt,nrayelt_emis,ifreq0,i_ox_conversion,
     +              ifreq_write,iray_status_one_ray

cBH130508:      real*8, pointer ::
cBH130508:c-----data for the emission calculations     
cBH130508:     & wi_0(:,:),                    !(nray,nfreq)
cBH130508:     & wtemp_rad_fr_wall(:),         !(nfreq)
cBH130508:     & wtemp_rad_fr(:),              !(nfreq)
cBH130508:     & wtemp_pl_fr(:),               !(nfreq)
cBH130508:     & wr0_em(:),                    !(nfreq)
cBH130508:     & wz0_em(:),                    !(nfreq)
cBH130508:     & wrho0_em(:),                  !(nfreq)
cBH130508:     & wfreq(:),                     !(nfreq)
cBH130508:     & wtau_em(:,:),                 !(nray,nfreq)
cBH130508:     & wi_0t(:,:),                   !(nray,nfreq)
cBH130508:     & wr_2nd_harm(:),               !(nfreq)
cBH130508:     & wtemp_2nd_harm(:),            !(nfreq)
cBH130508:cBH130508:c-------------------to plot cold N_perp in initial ray point for test
cBH130508:     & wn_perp_ioxm_p(:),            !(nfreq)
cBH130508:     &wn_perp_ioxm_m(:),             !(nfreq)
cBH130508:     &wn_perp_ioxm_n_npar_p(:),      !(nfreq)
cBH130508:     &wn_perp_ioxm_n_npar_m(:),      !(nfreq)
cBH130508:     &wye_0(:),                      !(nfreq)
cBH130508:     &wxe_0(:)                       !(nfreq)

      real*8
c-----data for the emission calculations     
     & wi_0,                    !(nray,nfreq) !BH130508:Changed to (nraya,nfreqa)
     & wtemp_rad_fr_wall,         !(nfreq)      !etc.
     & wtemp_rad_fr,              !(nfreq)
     & wtemp_pl_fr,               !(nfreq)
     & wr0_em,                    !(nfreq)
     & wz0_em,                    !(nfreq)
     & wrho0_em,                  !(nfreq)
     & wfreq,                     !(nfreq)
     & wtau_em,                 !(nray,nfreq)
     & wi_0t,                   !(nray,nfreq)
     & wr_2nd_harm,               !(nfreq)
     & wtemp_2nd_harm,            !(nfreq)
cBH130508:c-------------------to plot cold N_perp in initial ray point for test
     & wn_perp_ioxm_p,            !(nfreq)
     &wn_perp_ioxm_m,             !(nfreq)
     &wn_perp_ioxm_n_npar_p,      !(nfreq)
     &wn_perp_ioxm_n_npar_m,      !(nfreq)
     &wye_0,                      !(nfreq)
     &wxe_0                       !(nfreq)

cBH130508: The common block /write/ was slightly rearranged to make
cBH130508: counting easier for the MPI subroutines sendrec_alloc,
cBH130508: senddata, recvdata.
cBH130508: 
cBH130508: Much more reorganzation of this data structure is warrented.

c***********************************************************************
c    ====>>>>     ***WARNING***     <<<<====
c    If you change this common block dimensions, then these changes need
c    to be accounted for in ./mpi/mpi.ins in order to maintain MPI
c     functionality. 
c***********************************************************************

      common/write_i02/             !-------------------> integers
     +             nrayelt,
     +             nrayelt_emis,
     +             ifreq0,    
     &             i_ox_conversion,    
     &             ifreq_write,
     &             iray_status_one_ray,     
     &             w_ires(nrelta, n_relt_harm1a:n_relt_harm2a)
                        ! nrelta*(n_relt_harm2a-n_relt_harm1a+1)
                           
                           
      common/write_c13/             !-------------------> complex*16
     &             cwexde(nrelta),
     &             cweyde(nrelta),
     &             cwezde(nrelta),            ! 3*nrelta complex*16
     &             w_ceps(3,3,nrelta)         ! 9*nrelta complex*16



      common/write/             ! <<<<==== Begin of common block /write/
csend-recvdata: begin r1()
     &ws(nrelta),
     &wt(nrelta),                   !RK time
     &seikon(nrelta),
     &spsi(nrelta),
     &wr(nrelta),
     &wphi(nrelta),wz(nrelta),wnpar(nrelta),
     3              wnper(nrelta),delpwr(nrelta),sdpwr(nrelta),
     4              wdnpar(nrelta),fluxn(nrelta),sbtot(nrelta),
     4		    sb_z(nrelta),sb_r(nrelta),sb_phi(nrelta),
     5              sene(nrelta),ste(nrelta),szeff(nrelta),
     5              salphac(nrelta),salphal(nrelta),
     5              wye(nrelta),wyi(nrelta),wyi2(nrelta),
     5              wxi(nrelta),wxi2(nrelta),
     6              xarr(nrelta),yarr(nrelta),
     7              rez(nrelta),
     9              eff(nrelta),wmtor(nrelta),
     &              wn_z(nrelta),wn_r(nrelta),wn_phi(nrelta),
     &              wvgr_z(nrelta),wvgr_r(nrelta),
     &              wvgr_phi(nrelta),
     &              wtheta_pol(nrelta),    ! 39 real*8 dimension nrelta
                                           ! variables.
csend-recvdata: begin r2()
     &              salphas(nrelta,nbulka),! 1 real*8 dim nrelta*nbulka 
     +              freqncy0,              ! 1 real*8
c    tokman flux
     &              wflux_tokman_z(nrelta), !z,r,phi components   
     &              wflux_tokman_r(nrelta),wflux_tokman_phi(nrelta),

c------------------------------------------------------
c    the data for the emission calculations           
     +              wcnz_em(nrelta),wcnr_em(nrelta),wcm_em(nrelta),
     +              wz_em(nrelta),wr_em(nrelta),wphi_em(nrelta),  
     +              wal_emis(nrelta),wj_emis(nrelta),            
     +              wnray(nrelta),wsn(nrelta),     
     +              win_sn(nrelta),win_0(nrelta),
     +              wi_0sn(nrelta),
     +              wtemp_em(nrelta),wtemp_rad_em(nrelta),
     +              wtaun_em(nrelta),      ! 19 real*8 dimension nrelta
c    The following pointered variables are dimensioned in dinit.f.
     &wtau_em(nraya,nfreqa),           !pointer (nray,nfreq),
     &wi_0t(nraya,nfreqa),             !pointer (nray,nfreq),
     &wi_0(nraya,nfreqa),              !pointer (nray,nfreq),  !3 nray*nfreq real*8
     &wtemp_rad_fr_wall(nfreqa), !pointer (nfreq),
     &wtemp_rad_fr(nfreqa),      !pointer (nfreq),
     &wtemp_pl_fr(nfreqa),       !pointer (nfreq),
     &wr0_em(nfreqa),            !pointer (nfreq),
     &wz0_em(nfreqa),            !pointer (nfreq),
     &wrho0_em(nfreqa),          !pointer (nfreq),
     &wfreq(nfreqa),             !pointer (nfreq),
     &wr_2nd_harm(nfreqa),       !pointer (nfreq),
     &wtemp_2nd_harm(nfreqa),    !pointer (nfreq),
c-------------------to plot cold N_perp in initial ray point for test
     &wn_perp_ioxm_p(nfreqa),        !pointer (nfreq),
     &wn_perp_ioxm_m(nfreqa),        !pointer (nfreq),
     &wn_perp_ioxm_n_npar_p(nfreqa), !pointer (nfreq),
     &wn_perp_ioxm_n_npar_m(nfreqa), !pointer (nfreq),
     &wye_0(nfreqa),                 !pointer (nfreq),
     &wxe_0(nfreqa),                 !pointer (nfreq),       !15 nfreq real*8
c-------------------------------------------------------
c    for emission at O-X ebw jump case
c-------------------for O mode ray part 
     +              wsn_o(nrelta),
     +              wz_o(nrelta),wr_o(nrelta),wphi_o(nrelta),
     +              wcnz_em_o(nrelta),wcnr_em_o(nrelta),
     +              wcm_em_o(nrelta),
c-------------------for X_EBW mode ray part 
     +              wsn_x(nrelta),
     +              wz_x(nrelta),wr_x(nrelta),wphi_x(nrelta),
     +              wcnz_em_x(nrelta),wcnr_em_x(nrelta),
     +              wcm_em_x(nrelta),
c-------------------------------------------------------
c   the data for the determination of the N_par boundaries
     +		    wnrho(nrelta),gnpar(nrelta),wxe(nrelta),wp(nrelta),
     +		    wnparplb(nrelta),wnparmnb(nrelta),  !20 nrelta real*8
c-------------------------------------------------------
     +              phiold,
     +              zold,rold,rhoold,cld,cvac,          !6 real*8
c--------------------------------------------------------
c     These data are for OX transmission coefficient       
csend-recvdata: begin r3()
     &              transm_ox,cn_par_optimal,cnpar_ox,cn_b_gradpsi, !4 r8
c-------------------------------------------------------
c     These data are to plot the resonance ellipse boundary
c     for the several harmonics    
csend-recvdata: begin r4()
     &              wp_perpmax_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a),
     &              wp_parmin_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a),
     &              wp_parmax_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a),
     &              wp_perpmax_dmc(nrelta,n_relt_harm1a:n_relt_harm2a),
     &              wp_parmin_dmc(nrelta,n_relt_harm1a:n_relt_harm2a),
     &              wp_parmax_dmc(nrelta,n_relt_harm1a:n_relt_harm2a),
cBH130508     &              w_ires(nrelta,n_relt_harm1a:n_relt_harma),
cSm060303
     &              wp_0_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a),
     &              wdel_p_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a),
                           !8 r8 nrelta*(n_relt_harm2a-n_relt_harm1a+1)


c-------------------------------------------------------------------
c     data for power absorption at reflections 
csend-recvdata: begin r5()
     &              w_tot_pow_absorb_at_refl,         !1 r8
c-----data to plot delta power_e and delta current along the ray
     &             delpow_e_ar(nrelta),delcur_par_ar(nrelta)  !2*nrelta r8
c-------------------to plot data in initial ray point 



c the data for mnemonic.txt file
c nrelta is the max number of the ray elements along any ray
c---------------------------------------------------------
c     ws poloidal ray distance (cm)
c----------------------------------------------------------
c   wal_emis is the absorption(1/cm)           along the ray
c   wj_emis  is the emmisivity (egr*sec/cm**3) along the ray
c   wcnz_em,wcnr_em,wcm_em   the refractive index coordinates along the ray
c   wnray_em           the ray refractive index N_ray along the ray  
c   win_sn_em          emission I_n at the detector side of nth bin at s=s_n
c   win_0_em           emission at the plasma boundary s=s_1 from one nth bin s_n
c   wi_0_em            emission I_0 at the plasma boundary from the ray at s=s(1)=0
c   wsn_em total ray distance (cm)
c   wtemp_em           temperature along the ray
c   wtemp_rad_em       radiation temperature along the ray from the ratio
c                      wtepm_rad_em=(2.d0*pi)**3*(clight)**2/(omega*cnray)**2
c                     *j_emis/al_emis/1.6022d-9
c
c            
c   nrayelt_emis       the number of points along the ray with the additional
c                      emission points
c   wi_0sn(n)          sum{k=1,n}[in_0(k)]! emission at the plasma boundary
c                      from the part 0<s<sn(n) of the ray
c  wtaun_em(n)         tau_n from n-th bin (for ploting only)
c  wfreq(nfreqa)       the array for the frequencies.It is used for the
c                      emission calculations
c  ifreq0              is the number of the frequency in wfreq(ifreq) more
c                      closed to the second elecron gyro-frequency
c                      at plasma center
c                      wfreq(ifreq0) gives min{dabs(wfreq(ifreq)-2*freqncy0)}
c
c  wtemp_rad_fr()       radiation temperature from multi-frequence case
c                      =2.d0*pi*clight**2/1.6022d-9*wi_0_em()/omega**2       
c
c  wtau_em(nraya,nfreqa) is the tau(optical length) from one ray pass
c
c  wi_0t               is the flux wi_0 divided by the coefficient with
c                      wallr wall reflection coefficient
c  wtemp_rad_fr_wall   radiation temperature from malti freq. case
c                      with wallr coefficient
c  wtemp_pl_fr(nfreqa) the plasma temperature at the point along the ray
c                      where In_0/ds has the maximal value 
c
c  wr0_em(nfreqa)      the space coordinates asnd the small radius 
c  wz0_em(nfreqa)      of the ray point where In_0/ds has the maximal value 
c  wrho0_em(nfreqa)
c
c  wr_2nd_harm(nfreqa) the major radius of EC 2nd harmonic resonace points
c                       at Z=0  
c  wtemp_2nd_harm(nfreqa) the bulk plasma temperature at EC 2nd harmonic
c                         points T(z=0,wr_2nd_harm,phi=0)
c
c  w_ceps(3,3,nrelta)   complex dielectric tensor along the ray, for checking 
c
c  wn_z(nrelta),wn_r(nrelta),wn_phi(nrelta) are the refructive index
c                                           coordinates N_r,N_z,N_phi
c
c  wvgr_z(nrelta),wvgr_r(nrelta),wvgr_phi(nrelta) are the group velocity
c  components normalized to c

c  Tokman wave energy flux vector normalized to c 
c  wflux_tokman_z(nrelta)
c  wflux_tokman_r(nrelta)
c  wflux_tokman_phi(nnrelta)

c  For iabsorp=3, individual ion absorption damping coeffs.
c                 [Could be extended to some additional iabsorp]:
c                  \salphas(nrelta,nbulka)
c-----------------------------------------------------------
c          transm_ox OX transmission coefficient
c          i_ox_conversion=1 was the jump in the radial direction
c                         =0 was not OX conversion
c----------------------------------------------------------------
c   cn_par_optimal=dsqrt(Y_abs/(Y_abs+1)) is the optimal N parallel
c                  for OX conversion. It is used for mnemonic.nc file
c-----------------------------------------------------------------
c   cnpar_ox is N parallel before OX conversion procedure.
c            It is used to prepare the data for mnemonic.nc file
c-----------------------------------------------------------------
c   cn_b_gradpsi=N^*[b^*gradpsi^]/|[b^*gradpsi^]| is the refractive index
c                component before OX conversion.
c                It is perpendicular to the magnetic field b^ and
c                to the gradiend(psi)
c----------------------------------------------------------------
c     These data are to plot the resonance ellipse boundary
c     for the several harmonics    
c
c wp_perpmax_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a), maximal value of
c                                                     p_perp/ (m_e*V_thermal)
c                                                     at resonace allopse
c
c wp_parmin_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a)  minimal and the maximal 
c wp_parmax_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a) p_parallel/(m_e*V_thermal)
c                                        at resonance ellipse      
c
c w_ires(nrelta,n_relt_harm1a:n_relt_harm2a) =0 no resonance
c                              1 ellipse
c                              2 parabola
c                              3 hyperbole
c ellipse center for ires=1 case
c wp_0_dmvt(nrelata,n_relt_harm1a:n_relt_harm2a),
c shift from the ellipse center for  ires=1 case
c wdel_p_dmvt(is,n)nrelata,n_relt_harm2a:n_relt_harm2a)
c
c-------------------to plot cold N_perp in initial ray point for test
c                   wn_perp_ioxm_p(nfreqa) is nperp at ioxm=+1
c                   wn_perp_ioxm_m(nfreqa) is nperp at ioxm=-1
c                   wn_perp_ioxm_n_npar_p(freqa) is nperp at ioxm_n_npar=+1
c                   wn_perp_ioxm_n_npar_m(nfreqa) is nperp at ioxm_n_npar=-1
c                   wye_0(nfreqa) ye at initial poin
c                   wxe_0(nfreqa) xe at initial point
c                   ifreq_write
c------------------------------------------------------------------------
c w_tot_pow_absorb_at_refl total power [MWatt] absorbed at all reflections
c                        at all rays
c------------------------------------------------------------------------
c-----data to plot delta power_e and delta current along the ray
c                  delpow_e_ar(nrelta),delcur_par_ar(nrelta)
c---------------------------------------------------------------------------
c     wtheta_pol(nrelta) poloidal angle [degree] -180< wtheta_pol<+180
c---------------------------------------------------------------------------
c  iray_status_one_ray   giving a status code for each type of stopping
c                        of one ray
c                     =1  t.gt.poldist_mx
c                     =2  in prep3d delpwr(is).lt.delpwrmn*delpwr(1)
c                     =3  irefl.ge.ireflm
c
c                     =4  ray is at toroidal limiter boundary
c                         for max_limiters.ge.1
c                     =5  nstep_rk.gt.maxsteps_rk
c                     =6  nmode.gt.cnmax
c                         cnmax was set in output.f cnmax=10000.d0
c                     =7 
c                         for ((id.eq.1,2).and.(uh.gt.1.d0).and.
c                         ((uh-1.d0).lt.del_uh)).and.(cnper.gt.cnper_max_ebw))
c                         The ray is close to upper hybrid resonance
c                         It can be Xmode to close to the UH resonance'
c                     =8  if ((vgrmods.gt.1.1).and. 
c                         ((id.ne.10).and.(id.ne.12).and.(id.ne.13)))
c                         and 
c                         if (id.eq.14)
c                     =9  D/(N|gradD|) > toll_hamilt
c                     =10 in prep3d.f argexp>0
c                     =11 in prep3d fluxn(is).lt.0.0
c                     =12 nrayelt=nrelt
c                     =13 in RK drkgs2 small time step h.lt.1.d-11 
c                         can not get the given accuracy prmt4 
c------------------------------------------------------------------------
