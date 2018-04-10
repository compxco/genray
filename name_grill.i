c-------------------------------------------------------
c     namelist grill
c-------------------------------------------------------

      namelist /grill/ i_n_poloidal,n_theta_pol,ksi_nperp,
     *i_rho_cutoff,rho_step_find_LHFW_cutoff,
     *rho_initial_find_LHFW_cutoff, 
     *ngrilld, !the old name for ngrill for old genray.in files
     *ngrill,igrillpw,igrilltw,rhopsi0,thgrill,
     *phigrill,height,nthin,anmin,anmax,nnkpar,powers,
     &antormin,antormax,nnktor,anpolmin,anpolmax,nnkpol,
     &ilaunch,r0launch,z0launch,phi0launch,i_grill_pol_mesh,
     &i_grill_npar_ntor_npol_mesh
     
