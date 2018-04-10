
c the arrays for the text ray data for FP (CQL3D) code
c      nrelta is maximum value of nrelt
       double complex cwexde,cweyde,cwezde
       double precision wal_emis,wj_emis
       double complex w_ceps
       common/write/cwexde(nrelta),cweyde(nrelta),cwezde(nrelta),
     1              ws(nrelta),seikon(nrelta),spsi(nrelta),
     2              wr(nrelta),wphi(nrelta),wz(nrelta),wnpar(nrelta),
     3              wnper(nrelta),delpwr(nrelta),sdpwr(nrelta),
     4              wdnpar(nrelta),fluxn(nrelta),sbtot(nrelta),
     4		    sb_z(nrelta),sb_r(nrelta),sb_phi(nrelta),
     5              sene(nrelta),ste(nrelta),salphac(nrelta),
     5              salphal(nrelta),
     5              wye(nrelta),wyi(nrelta),wyi2(nrelta),
     5              wxi(nrelta),wxi2(nrelta),
     6              xarr(nrelta),yarr(nrelta),
     7              rez(nrelta),
     9              eff(nrelta),wmtor(nrelta),
     9              effpol(nrelta),efftor(nrelta),
     &              wn_z(nrelta),wn_r(nrelta),wn_phi(nrelta),
     &              wvgr_z(nrelta),wvgr_r(nrelta),
     &              wvgr_phi(nrelta),
c   for eps checking
     &              w_ceps(3,3,nrelta), 
c    tokman flux
     &              wflux_tokman_z(nrelta), !z,r,phi components   
     &              wflux_tokman_r(nrelta),wflux_tokman_phi(nrelta),
     &              salphas(nrelta,nbulka), 

c------------------------------------------------------
c    the data for the emission calculations           
     +              wcnz_em(nrelta),wcnr_em(nrelta),wcm_em(nrelta),
     +              wz_em(nrelta),wr_em(nrelta),wphi_em(nrelta),  
     +              wal_emis(nrelta),wj_emis(nrelta),            
     +              wnray(nrelta),wsn(nrelta),     
     +              win_sn(nrelta),win_0(nrelta),
     +              wi_0(nraya,nfreqa),wi_0sn(nrelta),
     +              wtemp_em(nrelta),wtemp_rad_em(nrelta),
     +              wtemp_rad_fr_wall(nfreqa),
     +              wtemp_rad_fr(nfreqa),
     +              wtemp_pl_fr(nfreqa),
     +              wr0_em(nfreqa),wz0_em(nfreqa),wrho0_em(nfreqa),
     +              wtaun_em(nrelta),wfreq(nfreqa),
     +              wtau_em(nraya,nfreqa),
     +              wi_0t(nraya,nfreqa),
     +              wr_2nd_harm(nfreqa),
     +              wtemp_2nd_harm(nfreqa),
     +              freqncy0,
c-------------------------------------------------------
c   the data for the determination of the N_par boundaries
     9		    wnrho(nrelta),gnpar(nrelta),wxe(nrelta),wp(nrelta),
     9		    wnparplb(nrelta),wnparmnb(nrelta),
c-------------------------------------------------------
     1              phiold,
     8              zold,rold,rhoold,cld,cvac,nrayelt,
     +              nrayelt_emis,ifreq0,
c--------------------------------------------------------
c     These data are for OX trasmition coefficient       
     &              transm_ox,cn_par_optimal,cnpar_ox,cn_b_gradpsi,
     &              i_ox_conversion


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
c  wflux_tokman_phi(nrelta)

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
