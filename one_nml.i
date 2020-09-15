
c This include file contains the variables which are used 
c in most of the subroutines and functions:  
c magnetic field; small radius; wave frequency;
c sin and cos of the angle gam between the refractive index and
c the magnetic field;
c the constants w and v  for the calculations Y and X;
c the characteristics of the profiles of the density, temperature, tpop,
c vflow, zeff;
c the characteristics for the toroidal density variations; 
c the values of the different switches
c
c It has only namelist data.

      real*8  
c-----from namelist /genr/ 
     &            r0x,b0,     
c-----from namelist /tokamak/           
     &            deltripl,psifactr,
     &            r_wall,z_wall,r_limiter,z_limiter,phi_limiter, 
     &            h_add_wall, 
c-----from namelist /wave/
     &            frqncy,delpwrmn,poldist_mx,cnperp_plot_min,
     &            cnperp_plot_max, 
     &            cN_perp_root_max,
     &            rho_step_find_hot_nperp_roots,
     &            rho_min_find_hot_nperp_roots,
     &            shift_rho_geom,
cSAP090304 move to edge_prof_nml
c     &            sigmedgn,sigmedgt,
c-----from namelist /dispers/
     &            del_y,ray_direction,errabs0,errrel0,diff_err,
     &            coll_mult,refl_loss,
c-----from namelist /numercl/
     &            prmt1,prmt2,prmt3,prmt4,prmt6,prmt9,
     &            uh_switch,prmt6_uh_switch,toll_hamilt,
     +            dL_step, dN_step, der_r, der_n, der_f,
     &            y_power_switch_resonance,
     &            prmt6_power_switch_resonance,
     &            del_y_power_switch_resonance,
     &            epsi,
     &            eps_delta_pow,
c-----namelist /output/
c-----namelist /plasma/
     &            temp_scale,den_scale,
c-----namelist /species/
c-----from namelist /varden/
     &            var0,denn,denm,an,sigman,
c-----from namelist /denprof/
     &            dense0,denseb,rn1de,
     &            rn2de,
c-----from namelist /tpopprof/ 
     &            tp0,tpb,rn1tp,rn2tp,
c-----from namelist /vflprof/
     &            vfl0,vflb,rn1vfl,
     &            rn2vfl,
c-----from namelist /tprof/ 
     &            ate0,ateb,rn1te,rn2te,
c                 namelist did not work at PC 
c                 with the names te0 and teb
c                 dimension ate0,ateb
c-----from namelist /zprof/
     &            zeff0,zeffb,rn1zeff,rn2zeff,
c-----from namelist /read_diskf/
     &            rtem0,
     &            rtail,r1t,r2t,ttail,rhot,r1h,r2h,thotpar,thotper,
     &            hotexp,hotmnpar,hotmxpar,hotmnper,hotmxper,
     &            rbeam,r1b,r2b,tbeam,ebeam,thbeam,
     &            rvtovte1,rvtovte2,rtemp1,rtemp2,rtemp3,
c-----from namelist /emission/
     &            tol_emis,freq00,freq01,wallr,    
c-----from namelist /grill/ 
     &            rho_step_find_LHFW_cutoff,
     &            rho_initial_find_LHFW_cutoff,
c-----from namelist /ox/
     &            theta_bot,theta_top,
     &            eps_antenna,
     &            eps_xe 
c---------------------------------------------------------
  
      integer
c-----from namelist /tokamak/           
     &            indexrho,ipsi,ionetwo,ieffic,nloop,i_ripple,
     &            n_wall,max_limiters,n_limiter,ieffic_mom_cons,  
c-----from namelist /wave/
     &            ioxm,ireflm,jwave,istart,ibw,i_vgr_ini,ioxm_n_npar,  
     &            n_nperp_plot,                                       
     &            n_points_root,             !number of  N_perp mesh points
     &            i_look_roots,k_hot_root,
     &            i_rho_find_hot_nperp_roots, no_reflection,
     &            istep_in_lcfs,
     
c-----from namelist /dispers/ 
     &            ib,id,iherm,iabsorp, iswitch,jy_d,idswitch,iabswitch,
     &            n_relt_harm,n_relt_intgr,iflux,i_im_nperp,
     &            i_geom_optic,navg,relres,iabsorp_collisional,
     &            n_relt_harm1,i_salphal,
     &            iabsorp_ql, 
c-----from namelist /numercl/
     &            ndim1,
     &            irkmeth,isolv,idif,nrelt,icorrect,i_output,
     &            maxsteps_rk,i_uh_switch,
     &            i_power_switch_resonance,
     &            n_power_switch_resonance,
     &            i_resonance_curve_integration_method,
c-----from namelist /output/
     &            iwcntr,iwopen,iwj,itools,i_plot_b,i_plot_d,
     &            i_plot_disp_cold,i_plot_wave_normal_cold, 
c-----from namelist /plasma/
     &            nbulk,izeff,idens,ndens,
c-----namelist /species/
c-----namelist /varden/
c-----namelist /denprof/
c-----namelist /tpopprof/
c-----namelist /vflprof/ 
c-----namelist /tprof/ 
c-----namelist /zprof/
c-----from namelist /read_diskf/
     &            i_diskf,
     &            jx,iym,lrz,ngen,
c-----from namelist /emission/ 
     &            i_emission,nfreq,nharm1,nharm2,
     &            i_rrind,i_r_2nd_harm,
     &            i_emission_spectrum,     
c-----namelist /grill/ 
c-----from namelist /ox/
     &            i_ox,     
     &            i_ox_poloidal_max
c-----------------------------------------
      character
c-----from namelist /numercl/
     &iout3d*20,           
c-----from namelist /genr/  
     &outdat*20,stat*3,     
     &outnetcdf*8, outprint*8, outxdraw*8,
     &mnemonic*128,                 
     &rayop*8,                        
     &dielectric_op*8,              
     &partner*32,                   
c-----from namelist /tokamak/  
     &eqdskin*256,          
c-----from namelist /dispers/
     &ion_absorption*8, 
c-----from namelist /plasma/
     &nonuniform_profile_mesh*8,    
c-----from namelist /read_diskf/
     &netcdfnm*128


      common /one_nml/
c-----real*8  
c-----from namelist /genr/ 
     &            r0x,b0,     
c-----from namelist /tokamak/           
     &            deltripl,psifactr, 
     &            r_wall(n_wall_a),z_wall(n_wall_a),
     &            r_limiter(n_limiter_a,max_limiters_a),
     &            z_limiter(n_limiter_a,max_limiters_a),
     &            phi_limiter(2,max_limiters_a), 
     &            h_add_wall, 
c-----from namelist /wave/
     &            frqncy,delpwrmn,poldist_mx,cnperp_plot_min,
     &            cnperp_plot_max, 
     &            cN_perp_root_max,
     &            rho_step_find_hot_nperp_roots,
     &            rho_min_find_hot_nperp_roots,
     &            shift_rho_geom,
cSAP090304 move to edge_prof_nml
c     &            sigmedgn,sigmedgt,
c-----from namelist /dispers/
     &            del_y,ray_direction,errabs0,errrel0,diff_err,
     &            coll_mult,refl_loss,
c-----from namelist /numercl/
     &            prmt1,prmt2,prmt3,prmt4,prmt6,prmt9,
     &            uh_switch,prmt6_uh_switch,toll_hamilt,
     +            dL_step, dN_step, der_r, der_n, der_f,
     &            y_power_switch_resonance(n_power_switch_resonance_a),
     &            prmt6_power_switch_resonance,
     &            del_y_power_switch_resonance,
     &            epsi,
     &            eps_delta_pow,
c-----namelist /output/
c-----namelist /plasma/
     &            temp_scale(nbulka),den_scale(nbulka),
c-----namelist /species/
c-----from namelist /varden/
     &            var0,denn,denm,an,sigman,
c-----from namelist /denprof/
     &            dense0(nbulka),denseb(nbulka),rn1de(nbulka),
     &            rn2de(nbulka),
c-----from namelist /tpopprof/ 
     &            tp0(nbulka),tpb(nbulka),rn1tp(nbulka),rn2tp(nbulka),
c-----from namelist /vflprof/
     &            vfl0(nbulka),vflb(nbulka),rn1vfl(nbulka),
     &            rn2vfl(nbulka),
c-----from namelist /tprof/ 
     &            ate0(nbulka),ateb(nbulka),rn1te(nbulka),rn2te(nbulka),
c                 namelist did not work at PC with the names te0 and teb
c                 dimension ate0(nbulka),ateb(nbulka) 
c-----from namelist /zprof/
     &            zeff0,zeffb,rn1zeff,rn2zeff,
c-----from namelist /read_diskf/
     &            rtem0,
     &            rtail,r1t,r2t,ttail,rhot,r1h,r2h,thotpar,thotper,
     &            hotexp,hotmnpar,hotmxpar,hotmnper,hotmxper,
     &            rbeam,r1b,r2b,tbeam,ebeam,thbeam,
     &            rvtovte1,rvtovte2,rtemp1,rtemp2,rtemp3,
c-----from namelist /emission/
     &            tol_emis,freq00,freq01,wallr,    
c-----from namelist /grill/ 
     &            rho_step_find_LHFW_cutoff,
     &            rho_initial_find_LHFW_cutoff,
c-----from namelist /ox/
     &            theta_bot(nconea),theta_top(nconea),
     &            eps_antenna,
     &            eps_xe, 
c---------------------------------------------------------  
c      integer
c-----from namelist /tokamak/           
     &            indexrho,ipsi,ionetwo,ieffic,nloop,i_ripple, 
     &            n_wall,max_limiters,n_limiter(max_limiters_a),
     &            ieffic_mom_cons,
c-----from namelist /wave/
     &            ioxm,ireflm,jwave,istart,ibw,i_vgr_ini,ioxm_n_npar,  
     &            n_nperp_plot,                                       
     &            n_points_root,             !number of  N_perp mesh points
     &            i_look_roots,k_hot_root,
     &            i_rho_find_hot_nperp_roots,
     &            no_reflection,
     &            istep_in_lcfs,

c-----from namelist /dispers/ 
     &            ib,id,iherm,iabsorp, iswitch,jy_d,idswitch,iabswitch,
     &            n_relt_harm,n_relt_intgr,iflux,i_im_nperp,
     &            i_geom_optic,navg,relres,iabsorp_collisional,
     &            n_relt_harm1,i_salphal(nbulka),
     &            iabsorp_ql, 
c-----from namelist /numercl/
     &            ndim1,
     &            irkmeth,isolv,idif,nrelt,icorrect,i_output,
     &            maxsteps_rk,i_uh_switch,
     &            i_power_switch_resonance,
     &            n_power_switch_resonance,
     &            i_resonance_curve_integration_method,
c-----from namelist /output/
c-----from namelist /output/
     &            iwcntr,iwopen,iwj,itools,i_plot_b,i_plot_d,
     &            i_plot_disp_cold,i_plot_wave_normal_cold, 
c-----from namelist /plasma/
     &            nbulk,izeff,idens,ndens,
c-----namelist /species/
c-----namelist /varden/
c-----namelist /denprof/
c-----namelist /tpopprof/
c-----namelist /vflprof/ 
c-----namelist /tprof/ 
c-----namelist /zprof/
c-----from namelist /read_diskf/
     &            i_diskf,
     &            jx,iym,lrz,ngen,
c-----from namelist /emission/ 
     &            i_emission,nfreq,nharm1,nharm2,
     &            i_rrind,i_r_2nd_harm,
     &            i_emission_spectrum,     
c-----namelist /grill/ 
c-----from namelist /ox/
     &            i_ox,     
     &            i_ox_poloidal_max,
c-----------------------------------------
c     character
c----------------------------------------
c-----from namelist /numercl/
     &iout3d,           
c-----from namelist /genr/  
     &outdat,stat, outnetcdf, outprint, outxdraw,    
     &mnemonic,                 
     &rayop,                        
     &dielectric_op,              
     &partner,                   
c-----from namelist /tokamak/  
     &eqdskin,                 
c-----from namelist /dispers/
     &ion_absorption,
c-----from namelist /plasma/
     &nonuniform_profile_mesh,    
c-----from namelist /read_diskf/
     &netcdfnm

c--------------------------------------------------------------------
c     a           - small radius(normalized)
c     b0          - characteristic magnetic field (tl)
c     q0,pq	  -
c     rho	  -normalized small radius
c     rho2	  -rho2=rho**2
c     psifactr    =<1, is used in gr2new to find the last closed flux surface
c                  psi(r,z)=psimag+(psilim-psimag)*psifactr
c     ax	  -
c     r0x         - characteristic radius(m)
c     frqncy      - wave friquency (GHz)
c     pi	  -	  4*datan(1.d0)
c     gam         - angle between magnetic field and refractive index(rad)
c     bz,br,bphi  - components of the magnetic field
c     bmod        - module of the magnetic field
c     psid        - poloidal flux
c     dpdzd,dpdrd - derivatives from psi by! i_look_roots=0    !do not plot D(N_perp) and do not calculate all 
!                   !hot roots
!             =1    !plot D(N_perp) and do calculate all 
!                   !hot roots and do not calculate ray   
!             =2    !calculate hot root, use the root with number k_hot_root
!                   !as the initial ray condition and calculate ray z and r
c     dbzdz,dbzdr,dbzdph -derivatives from z-component of magnetic field
c     dbrdz,dbrdr,dbrdph -derivatives from r-component of magnetic field
c     dbphdz,dbphdr,dbpdph- derivatives from toroidal magnetic field
      !YuP: careful! There is {dbphdz,dbphdr,dbpdph}. There is no dbphdph !
c     dqdrho
c     dbmdz,dbmdr,dbmdph -derivatives from the module of magnetic field
c     dc2dz,dc2dr,dc2dph -derivatives from cos(gam)**2 by z,r and phi
c     dc2dnz,dc2dnr,dc2dm -derivatives from cos(gam)**2 coordinates of the
c                         refractive index n_z,n_r and m
c     w(4) - (omega_p(j)/omega)**2 on the magnetic axis
c     w(4) - (omega_B(j)/omega)    on the magnetic axis
c     dense0,denseb - central and bound electron densities (in 10**13/cm**3)
c     den(rho)=denseb+(dense0-denseb)*(1-rho*rn1de)**rn2de
c     te0,teb - central and bound electron temperatures (keV)
c     tempe(rho)=teb+(te0-teb)*(1-rho**rn1te)**rn2te !average temperatura
c     tpop(rho)=tpb+(tp0-tpb)*(1-rho**rn1tp)**rn2tp  !T_perp/T_parallel
c     vfow(rho)=vflb+(vfl0-vflb)*(1-rho**rn1vfl)**rn2vfl  !drift velocity || B
c     zeff0,zeffb - central and bound effective charge
c     zeff0(rho)=zeffb+(zeff0-zeffb)*(1-rho**rn1zeff)**rn2zeff
c     powtott  - total power summed over antennas
c     powini   - initial power in ray_iray chanal
c     delpwrmn - Minimum power in each ray, as a fraction of
c                starting power in the ray, after which ray is
c                stopped.  
c     nbulk -number of the plasma components
c     ioxm - sign before square root in dispersion relation
c           N=N(theta)
c           root=(-b+ioxm*sqrt(b**2-4a*c))/2a
c     ioxm_n_npar - sign before square root in dispersion relation
c           N=N(N_parallel)
c           root=(-g+ioxm_n_npar*sqrt(g**2-4f*w))/2g
c     icount -
c     ib - index of plasma component for which delib=1-yib may be
c                equal zero inside the plasma,
c                dispersion relation is multiplied by delib
c     id - index of the form of the dispersion relation
c     isolv -index of the solution method of ray trasing equations
c     idif  -index of analytical(=1) or numerical (=2) differentiation
c                  of Hamiltonian
c     irkmeth- index of the modification of Runge-Kutta  procedure
c     ireflm - maximum number of ray reflections
c     irefl  - number of ray reflections   
c
c     no_reflection !=1 switch off the artificial reflection from 
c               !the last closed flux surface
c               !=0 (default) switch on the artificial reflection from 
c               !the last closed flux surface
c
c     sigmedgn =0.02 (default) normalized exponential density fall off dist
c               outside LCFS starting at rho=1 density, for no_reflection=1
c     sigmedgt =0.02 (default) normalized exponential temperature fall off dist
c               outside LCFS starting at rho=1 temperature, for no_reflection=1
c
c     istart - index to choose wave start is outside the plasma (EC)
c                    or inside plasma (LH,FW)
c     idens  - index of spline or quasiparabolic density,temperature
c                     and z_eff  approximation
c     indexrho -index to choose radial coordinate
c     ionetwo - index to calculate power and current drive profiles
c     jwave   -	number of the wave harmonic
c     i_      - number of the file for wr! i_look_roots=0    !do not plot D(N_perp) and do not calculate all 
!                   !hot roots
!             =1    !plot D(N_perp) and do calculate all 
!                   !hot roots and do not calculate ray   
!             =2    !calculate hot root, use the root with number k_hot_root
!                   !as the initial ray condition and calculate rayiting 3D data
c     i1_     - number of the file for writing ray data
c     izeff   =0,then zeff will be calculated from the given plasma
c              components; =1 then zeff are given, but the ions
c              components(density) will be calculated
c		iboundb  this flag controls (in subroutine b.for)
c              is the point inside the plasma or not:
c              =-1 the ray point is inside th   
c	       =2 r<rmin+epsbnd or r>rmax-epsbnd  (outside the plasma)
c	       =1 z<zlimiter_min(r) or z>zlimiter_plus(r)(outside the plasma)
c		drhodzb  is a derivative ! i_look_roots=0    !do not plot D(N_perp) and do not calculate all 
c                   !hot roots
c             =1    !plot D(N_perp) and do calculate all 
c                   !hot roots and do not calculate ray   
c             =2    !calculate hot root, use the root with number k_hot_root
c                   !as the initial ray condition and calculate ray of the small radius from z.
c              (For the case when the point is outside the plasma).
c		It was calculated in b.for.
c     iabsorp -choice of Imag(N_perp) 
c             (N_perp is the perpendicular refractive index)
c              iabsorp=1 for EC waves from Mazzucato solver
c               =2 for LH waves, =3 for FW waves 
c--------------------------------------------------------------------
c   iabsorp_ql
c   iabsorp_ql =0  do not use QL flux for absorption calculations
c   iabsorp_ql =1  to use QL flux for electron absorption calculations
c               In this case QL flux will be calculated for harmonics
c               numbers nharm in the following interval:
c               n_harm_adj_min =< nharm   =< n_harm_adj_max 
c               n_harm_adj_min   number of minimal and maximal harmonics
c               n_harm_adj_max   for power and CD calculations  
c               Electric field polarization will be calculated according to
c               the value of the index: iabsorp 
c               Energy flux "fluxn" will be calculated according to 
c               the value of the index: iflux 
c------------------------------------------------------------------
c For use with Abhay Ram\'s dispersion relation (id=14) parameters to
c control the integration routipne for the Trubnikov integral.
c
c relres offers 4 choices for the resolution to use for the relativistic
c dispersion relation using the Trubnikov integral:
c relres=1 low resolution, errabs0=1.d-4, errrel0=1.d-4, navg=3, diff_err=0.1
c       =2 medium res., errabs0=1.d-5, errrel0=1.d-5, navg=12, diff_err=1.d-3
c       =3 high res., errabs0=1.d-6, errrel0=1.d-6, navg=25, diff_err=1.d-6
c       =4 user-defined res., set the following parameters manually:
c The Trubnikov one-dimensional (complex) integral is performed by splitting
c up a region from 0 to 1.d8 into 10^6 pieces, and each piece is integrated
c using the SLATEC adaptive quadrature routine dqag. errabs0 and errrel0 are
c the absolute and relative error tolerances passed directly to dqaq.
c Then the adjacent pieces are compared (it is an oscillatory integrand)
c and using navg number of pieces, when the average difference between them
c are less then diff_err, the integration is presumed finished (Thus it may
c finish long before the upper limit of 1.d8).
c
c errabs0 - absolute error for dqag integration routine
c errrel0 - relative error for dqag integration routine
c navg - number of adjacent integration intervals to use in comparison
c diff_err - error tolerance using navg pieces, when the average difference
c         is less than diff_err, then the integration is done.
c--------------------------------------------------------
c       flux=B~.B+E~.d(omega*eps_herm)/(domega).E
c    iflux=1 the flux will be calculculated using the the group velocity from
c           the choosed disperson relation (with given id) and the electric
c           field calculated for choosed iabsorp                  
c    iflux=2 the flux will be calculated using V_gr for the electron cold
c            plasma dispersion and polarization (using subroutine  grpde2)     
c-------------------------------------------------------
c     ipsi    -=1 calculation of conturs psi(z,r)=const,
c              =0 reading these contours from psi.bin file 
c                 (see genray.in file)
c     dpsimax -if in eqdsk (psimag.gt.psilim) then dpsimax=-1. else dpsimax=1.
c               It was determined in equilib.for: subroutine input
c     nrelt     Maximum number of ray elements per ray.
c               Must be .le. nrelta (a parameter)
c     iwcntr   =1 genray will call mk_grapc(iray) for the  
c                 calculations of contours
c               wb_c=i,
c              =0 call mk_grapc will be commented
c iwcntr =1 genray.f will calculate the contours wb_c=n
c        =0 genray.f will not do it
c iwopen =1 mk_grapc will calculte open contours wb_c=n	(using contrb1)
c         2 mk_grapc will calculte close contours wb_c=n (using contrb2)
c iwj    =mk_grapc will calculate contours wb_cj=n, here j is a kind of plasma
c         component must be.le.nbulk, j=1 for the electron hirofrequency
c         j=2,  for the ion (j kind) hirofrequency
c ieffic  chooses the formula for the current drive efficiency
c        =1 asymptotic symple formula (homogenius, nonrelativictic)
c        =2 asymptotic symple formula (East-Karney )
c        =3 asymptotic symple formula (curba subroutine)
c ibw    =0 the launch is not for Bernstein waves
c        =1 the launch of electron Bernstein wave from dhot tensor
c           The last case works only for istart=2 and grill_lh conditions 
c------------------------------------------------------------------------
c The change of the dispersion relation and absorption
c near the gyro-frequency points
c-------------------------------------------------
c iswitch=1   To use the change of the dispersion relation and
c             absorption
c        =0   Do not use the change of the dispersion relation and
c             absorption 
c     del_y   If the difference |1-nY(jy)|<del_y 
c             (jy=1-nbulk ,n=...-2,-1,0, 1,2,3,...)
c             then switch on the new 
c             given type of the disperion and absorption.
c   jy_d      is the type of plasma species 1<=jy_d<=nbulk
c   idswitch  is the type of the dispersion function near the 
c             gyro-frequency points
c             It can be equal 1,2,3,4,5,6
c   iabswitch is the type of the absorption near the gyro-frequency point  
c-------------------------------------------------------------------
c   The parameters for the density variation     
c   var0,denn,denm,an,sigman (see x.f: vardens)
c   var0 is an amplitude of the density variation (del_n_0) (see 3.37...)
c   denm is the number of the poloidal mode in the density variation(l_theta)
c   denn is the number of the toroidal mode in the density variation(l_phi)
c   an   is the radial localization of the variation (rho_0)
c   sigman is the parameter that characterizes the radial thickness
c          of the density fluctuation    
c----------------------------------------------------------------
c   the switch of the correction in outpt
c   icorrect=0 switch off the correction
c           =1 switch on the correction
c----------------------------------------------------------------
c  itools =0 do not use mkgrtool in genray.f
c         =1 use mkgrtool in genray.f
c---------------------------------------------------------------
c  nstep_rk    is the number of the Runge-Kutta time step 
c
c  maxsteps_r  is the maximal number of the time steps of Runge-Kutta
c              solver
c----------------------------------------------------------------
c  i_diskf=0 no usage of the non-maxwellian electron distribution 
c         =1 reading the file diskf
c         =2 readinf the file netcdfnm.nc 
c         =3 analitical calculation of the non-maxwellian disribution
c         =4 analytic calculation of continuous non-Maxwellian distribution
c            with three temperatures in  three energy ranges
c-------------------------------------------------------------- 
c deltripl is an amplitude of the ripple magnetic field at the
c          last flux surface (at rho=1)
c nloop    is the number of the magnetic coils
c i_ripple is the index to chose the ripple model
c-------------------------------------------------------------
c iout3d    ='enable'  to write output 3d.dat wile   [DEFUNCT]
c           ='disable' not write 3d.dat file         [DEFUNCT]
c--------------------------------------------------
c i_vgr_ini =+1 the wave is directed into the plasma (in the initial point)
c           =-1 the wave is directed out the plasma (in the initial point) 
c--------------------------------------------------
c outdat is the name of the output file that contains the data along the ray
c        r,z,phi,n_r,n_z,m (usually it outdat='zrn.dat')
c stat   is the status of outdat file (usually stat='new')        
c netcdfnm !name of the input *.nc file with 3d non-maxwellian distribution
c outnetcdf='enabled' (by default), or 'disabled' (to suppress writing data to *.nc)
c outprint='enabled' (by default), or 'disabled' (to suppress printout to screen)
c outxdraw='enabled' (by default), or 'disabled' (to suppress saving files for xdraw)
c------------------------------------------------------------------------
c The data for the creation of the analytical
c non-Maxwellian electron distribution:
c     rtem0,
c     rtail,r1t,r2t,ttail,rhot,r1h,r2h,thotpar,thotper,
c     hotexp,hotmnpar,hotmxpar,hotmnper,hotmxper,
c     rbeam,r1b,r2b,tbeam,ebeam,thbeam
c     jx,iym,lrz,ngen   
c The data fornalytic calculation of continuous non-Maxwellian distribution
c with  three  temperaures in three energy  ranges
c     rvtovte1,rvtovte2,rtemp1,rtemp2,rtemp3
c--------------------------------------------------
c  n_relt_harm1 is the lowest, i.e., minimum harmonic used in the
c               anti-hermitian dielectric tensor calculations.
c               It man be positive of negative.
c               Default value is +9999, in which case this input
c               is ignored.
c   n_relt_harm (.ge.1) gives the number of EC harmonics used 
c               in anti-hermitian dielectric tensor calculations
c               If n_relt_harm1=9999, then harmonics from 
c                 -n_relt_harm to +n_relt_harm are used.
c               If n_relt_harm1.ne.9999, then harmonics from
c                  n_relt_harm1 to n_relt_harm1+n_relt_harm are used.
c     It is necessary that the harmonics used in this calculation
c     be within the range of parameters [n_relt_harm1a,n_relt_harm2a]
c-----------------------------------------------------------
c     the data for the emission calculations
c     i_emission=0 no emission
c               =1 emission calculations
c---------------------------------------------------------
c   0<tol_emis=<1 tolerance parameter to add the new mesh point s_n
c                if in_0(n)>tol_emis*i_0 
c--------------------------------------------------------
c   nfreq=<nfreqa (see param.i) is the number of frequences
c         nfreq=1 gives detailed plots of emission from a single ray
c         ngraq.gt.1 gives spectra covering the specified frequency
c          range
c   freq00=<freq01 are ratios of the minimal and maximal emission
c         frequencies to the central electron gyro-frequency
c         f_ce(rho=0) 
c------------------------------------------------------------- 
c   nharm1=<nharm=<nharm2 give the number of used EC emission haremonics 
c
c   wallr is the wall reflection coefficient
c
c     i_rrind chooses the subroutine to calculate N_ray 
c     i_rrind=0  N_ray=N
c     i_rrind =1 from cold electron plasma using rrind   
c     i-rrind =2 from hot non-relativistic dispersion relation 
c---------------------------------------------------------------
c     i_r_2nd_harm=1 to calculate the major radius of the EC 2nd harmonic
c                 =0 do not calculate (default =0) 
c---------------------------------------------------------------
c poldist_mx is the maximal polidal distance (cm) along the ray
c            default=1.d+5
c---------------------------------------------------
c   temp_scale(ncomp),den_scale(ncomp) are the parameters to multiply
c   the given temperature and density profiles
c--------------------------------------------------------------------
c i_plot_b =1 create figures for the magnetic field,density and temperature 
c             profiles in plot.ps file using subroutine map_b based on PGplot
c i_plot_b =0 do not write the b,n,T figures to plot.ps file 
c
c i_plot_d =1 create the dispersion function contours d(ReNperp,ImN_perp)
c             in plot.ps file using PGplot
c i_plot_d =0 do not create the dispersion function contours
c-----------------------------------------------------------
c i_im_nperp chois of the metode to find Im_N_perp for hot plasma (iabsorp=4)
c-----------------------------------------------------------------------
c i_im_nperp=1 Im_N_perp=abs(ImD_full/(dD_hermitian/dReN_perp)) 
c i_im_nperp=2 (Re_N_perp,Im_N_perp) is the complex root 
c              (of the complex dispersion relation)
c              calculated by Newton iterations with the numerical
c              derivatives (the chord method)
c---------------------------------------------------------------
c  ndens is the number of point for the radial profiles 
c        density and temperature
c---------------------------------------------------------------
c k_root,       the number of the dispersion tensor
c               eigenvalue, used for id=10, determined in cninit.f
c----------------------------------------------------------- 
! i_geom_optic sets  the form of the ray equations
!              =1  integration in time (default):
!                  ray-tracing equations right hand side=
!                  dr^/dt=-(dD/dN^)/(dD/domega)
!                  dN^/dt=+(dD/dr^)/(dD/domega)
!                  In this case rside1 gives v_group.
!              =2  integration is space,
!                  ray-tracing equations right hand side=
!                  dr^/dl=- ray_direction * (dD/dN^)p
!                  dN^/dl=  ray_direction * (dD/dr^)p
c                  p=1.d0/dsqrt(deru(1)**2+deru(2)**2+(r*deru(3))**2)=
c                  deru(1)=dD/dN_z,deru(2)=dD/dN_r,deru(3)=dD/dCM,
c                  N_phi=cm/r
c----------------------------------------------------------------------
c ray_direction =+1 as default it
c               -1
c It is a multiplier in right hand side of ray-tracing equations
c It is used for i_geom_optic=2 case
c----------------------------------------------------------------------
c i_ox =0 /default/ do not use the vertex coordinates calculations 
c      =1 use EC cone vertex cordinates calculations
c      =-2 creates ray jump at OX conversion  area
c-----------------------------------------------------------------------------
c theta_bot(icone)<theta_top(icone) 
c are the poloidal angle boundaries (degree) at
c Xe=1 surface.They are used in genray.f to find the optimal ray
c direction at the given EC cone vertex icone=1,...,ncone
c------------------------------------------------
c eps_antenna       is the accuracy that the ray launched from the
c                   flux surface Xe(rho)=1 at the given poloidal angle
c                   reached the given cone vertex.
c                   This accuracy is calculated using the difference of 
c                   the vertex radial coordinate:
c                   dabs(r_st_ox-r_st)< eps_antenna
c------------------------------------------------
c i_ox_poloidal_max is the maximal number of the used poloidal angles
c                   in the bisection method. This method calculates the
c                   the poloidal angle theta_pol of the point M at the
c                   flux surface Xe(rho_=1 : M(poloidal_angle,rho).  
c                   The ray lauched from M to the plasma edge will
c                   go to the EC cone vertex.
c---------------------------------------------------------
c eps_xe   Is the parameter which sets the vicinitity
c          of the O-mode cutoff surface
c          If xe < (1-eps_xe) then the subroutine ox_conversion will 
c          creates the ray jump in small radius direction
c          and finds X mode.
c-------------------------------------------------------------------            
c partner  = 'disabled' to use input profiles from genray.dat/genray.in
c          = 'onetwo'   to use input plasma profile data from the file: 
c                       genray_profs_in created by transport code
c                       and to write the output file of power absorption
c                       and current drive genray_profs_out
c                       for the transport code.
c-------------------------------------------------------------------
c i_output =1 output at equal poloidal distance
c          =2 output at equal total distance
c------------------------------------------------------------------
c prmt6 time period or poloidal distance for output
c------------------------------------------------------------------
c i_uh_switch=1    if uh=dsrt(xe+ye**2) < uh_switch then change 
c                  the output step prmt(6) for prmt6_uh_switch 
c            =0    do not change the output step prmt(6)
c prmt6_uh_switch  is the output step for i_uh_switch=1 case
c uh_switch        if uh<uh_switch then change the output step for    
c                  i_uh_switch=1 case
c---------------------------------------------------------------------
c      iabsorp_collisional =0 no additional collisional absorption
c                          =1 collisional absorption  using formula
c                             Im(N)=dabs(nu_ei/(gr_perp))*clight/omega)
c      coll_mult =1.d0(default), multiplies above coll absorp expression
c---------------------------------------------------------------------------
c i_salphal (nbulka)  sets which damping will be in salphal_nc
c	     Default:i_salphal(1)=1,i_salphal(2:nbulk)=0 electron damping only  
c            A particular species contribution to salphal_nc is added if
c            i_salphal(species_number)=1. That is, damping coefficients
c            for all species with i_salphal(k).ne.0 are summed into saplhal_nc 
c--------------------------------------------------------------------------- 
c i_plot_disp_cold    ! =1 to plot D_cold(N_perp) to plot.ps file 
c                     ! =0 do not plot  
c                     ! It is used only  in grill_lh to create
c                     ! the plot in the initial point
c-----------------------------------------------------------------------
c i_plot_wave_normal_cold! for plotting cold wavenormal to plot.ps file
c                        !  =1 to plot wavenormal
c                        !  =0 do not plot
c------------------------------------------------------------------------
c nonuniform_profile_mesh= 'enabled' use nonuniform small radius mesh for input
c                                spline profiles (works for idens=1 only)
c                = 'disabled'    do not use nonuniform mesh
c--------------------------------------------------------------------------------  
c  Calculation of the small radius value near the plasma edge
c  where LH or FW have cutoff:  
c  i_rho_cutoff=0 (default) no calculations
c              =1 use these calculations
c--------------------------------------------------------------------
c  rho_step_find_LHFW_cutoff  is the non dimensional small radius step
c                            used in subroutine  rho_ini_LHFW                            
c                            It is used at i_rho_cutoff=1
c---------------------------------------------------------------------
c  rho_initial_find_LHFW_cutoff  is the initial small radius point. 
c                            As default rho=1- rho_step_find_LHFW_cutoff
c                            It is used at i_rho_cutoff=1
c--------------------------------------------------------------------------------
c           i_emission_spectrum      !to calculate emisssion spectrum   
c-----------------------------------------------------------------------------
c           cnperp_plot_min,cnperp_plot_max !max and min Nperp to plt D(Nperp)      
c           n_nperp_plot,                   !number of Nperp points to plot D(Nperp)
c-----------------------------------------------------------------------------------
c            cN_perp_root_max               !max value of n_perp to
c                                           !seek hot roots 
c            n_points_root                  !number of  N_perp mesh points
c                                           !to find hot plasma roots
c-------------------------------------------------------------------------------------
c i_look_roots=0    !do not plot D(N_perp) and do not calculate all 
c                   !hot roots
c             =1    !plot D(N_perp) and do calculate all 
c                   !hot roots and do not calculate ray   
c             =2    !calculate hot root, use the root with number k_hot_root
c                   !as the initial ray condition and calculate ray
c k_hot_root        is the number of the hot plasma root
c                   N_perp_root_ar(k_hot_root)
c                   which will be used for ray initial condition
c                   It works for i_look_roots=2 case only
c-------------------------------------------------------------------------
c For subroutine  rho_ini_hot_nperp_roots to find hot plasma roots
c i_rho_find_hot_nperp_roots,
c rho_step_find_hot_nperp_roots,
c rho_min_find_hot_nperp_roots,
c--------------------------------------------------------------
c  measure error in the dispersion realation
c  If D/(N|gradD|) > toll_hamilt stop ray calculation
c
c  toll_hamilt 
c--------------------------------------------------------------------------
c
c i_power_switch_resonance   =1  to use  prmt6_power_switch_resonance
c                            =0  do not change the output step prmt(6)
c
c prmt6_power_switch_resonance   is the output step for 
c                                i_power_switch_resonance=1 case
c                                 in resonace area
c
c n_power_switch_resonance   is the number of different used EC resonance 
c
c y_power_switch_resonance(n_resonance_a) are used ratios omega_ce/omega
c
c del_y_power_switch_resonance determines the resonance area
c The condition for resonance area 
c abs(Ye-y_power_switch_resonance(k))< del_y_power_switch_resonance
c k=1,...,n_power_switch_resonance 
c n_power_switch_resonance_a is max of  n_power_switch_resonance
c                            It is set in param.i
c-------------------------------------------------------------------
c       i_resonance_curve_integration_method=1
c       i_resonance_curve_integration_method=2
c       i_resonance_curve_integration_method=3 !trapezoidal 
c       i_resonance_curve_integration_method=4 !adaptive Simpson 
c       It used in subroutine intgr_rl to choose the numerical integration
c       method  o calculate anti-hermitian relativistic dielectric tensor
c-------------------------------------------------------------------
c       epsi is the accuracy used in integration by adaptive Simpson
c-------------------------------------------------------------------
c       eps_delta_pow   ! If the power variation along a ray at one output step
c                       ! bigger than eps_delta_pow:
c                       ! (delpwr(is-1)-delpwr(is))/delpwr(is-1) > eps_delta_pow
c                       ! and if the ray power is bigger than initial ray power
c                       ! multiplied by eps_delta_pow:
c                       ! delpwr(is-1)> delpwr(1)*eps_delta_pow 
c                       ! then the code will reduce the output step prmt(6).
c                       !
c                       ! Then the code will return the output step prmt(6)
c                       ! to the initial value prmt6 given in the input file.
c--------------------------------------------------------------------------
c      ion_absorption ="enabled"  add ion absorption
c                     ="disabled" no ion absorption
c------------------------------------------------------------------
c      to read wall and limiter positions
c
c      n_wall number of wall points =< n_wall_a
c             if n_wall= 0 then
c                no wall to be used, no reflection from the wall
c
c             if n_wall>0 then 
c                wall coordinates"  r_wall, z_wall,
c                and reflect rays from straight line segments between
c                the given points.  Count the number of coord pairs.
c                The coords must begin and end at the same physical
c                point.
c
c      max_limiters number of limiters. It should be =< max_limiters_a
c      n_limiter(1:max_limiters) numbers of limiter points at each limiter
c             if n_limiter=0 the
c                 no limiter to be used, no reflection from the limiter
c             if n_limiter>0 then 
c                Reflect from the
c                chamber wall consisting of the wall coordinates and
c                those limiter coordinates to the plasma side of the
c                wall.
c
c      r_wall(n_wall_a),z_wall(n_wall_a) wall coordines [m]
c      r_limiter(n_limiter_a,max_limiters_a),
c      z_limiter(n_limiter_a,max_limiters_a) limiters coordines [m]
c
c      phi_limiter(2,max_limiters_a)       limiter toroidal angles [degree]
c
c      h_add_wall    the distance between wall points of the wall mesh
c                    with additional points
c------------------------------------------------------------------------------

