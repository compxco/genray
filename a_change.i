
c     a_change.i
c
c
c***********************************************************************
c
c     This file records changes in the code 
c
c
c***********************************************************************

!------------------------------------
[233] version="genray_v12.1_260114"
!------------------------------------

[233] Corrected several bugs. Search [2026-01-14].
 Most important - in the call of subr. 
 calc_emission_spectrum() - 
 the argument wi_0 is replaced with wi_0(iray_loc,ifreq),
 It should be a scalar. 

!------------------------------------
[232] version="genray_v12.0_251217"
!------------------------------------

[232]  Added an option for refining the step of integration
  around cyclotron resonances.
  Presently it is done for irkmeth=2 (subr. drkgs2) 
  and irkmeth=3 (subr. drkgs_auto).
  The step is reduced around omega = n*omega_c resonances
  where the harmonic number is scanned in the range
  n= [1; nharm_refined_step], 
  where nharm_refined_step is the new namelist variable.
  Each species in the 1:nbulk list is checked.
  If condition (1.0 - n*|omega_c/omega|) .le. 0.1
  is met, then the step is reduced by 10x (hardwired, for now).
  Away from the resonance, the step is restored.
  The advantage is that a user can start 
  with fairly large prmt4 and prmt6 
  (for irkmeth=2, i.e. subr.drkgs2),
  but around resonance they are reduced.
  Similarly - for irkmeth=3 (subr. drkgs_auto)
  where dL_step is refined.
  As a result, the number of total steps
  (ray elements) is reduced.
  This is especially useful for runs in cold plasma
  where the effective width of power absorption
  around resonance is very narrow.
  See [2025-12-15], [2025-12-16]

!------------------------------------
[231] version="genray_v11.2_240814"
!------------------------------------

[231] Added new namelist variables to enable capability for using
  3D profiles for density, temperature and Tpar/Tperp values
  over general (Z,R,phi) grid.
  See description for model_rho_dens=7
  in read_write_genray_input.f.
  Track changes related to this option by searching "[2024-08-14]".

!------------------------------------
[230] version="genray_v11.1_240725"
!------------------------------------

[230] Fixed subr. absorpfd_091016,
  to avoid NaN. Maybe need more fixing...
  Search for "YuP[2022-08-17]", see details in comments. 

[229] Corrected functioning of ion_absorption.eq.'disabled' option
  (used for test purposes). For disabled ion absorption
  need to set cnprim_i=0.d0 in coding. See "YuP[2022-04-22]" for details.

[228] Impose a lower limit of 1.0 for Zeff after 
   call read_nonuniform_line_profile() which makes splines for profiles.
   Reason: spline-related "oscillations" in Zeff(rho)
   when tabulated values have a sudden change from one radial
   point to next point. [2022-03-22]

[227] Fixed a bug in subr.grill_lh, for the option igrillpw=1.
   This option is supposed to give a flat power spectrum over Npar0.
   The problem was occuring when the number of points in Npar0
   spectrum, nnkpar, was set to an odd number. In this case,
   the power for the central point in Npar0 spectrum
   was set to a larger value comparing to other points.
   See YuP[2022-03] in grill_lh.f
 

[226] In subr. efKarney (called for ieffic=2 setting), 
   added islofa==jwave=0 condition 
   for Landau damping (was: jwave.ne.-1 for Landau damping);
   For any other jwave, except -1 and 0, print warning message
   before starting ray tracing. Subr.efKarney is not supposed
   to be called with jwave>0, but if it is called with jwave>0, 
   it will be using same settings as for jwave=0 
   (although results are not expected to be physical in this case).
   YuP[2021-02-11] 

[225] Added (cn2-cnpar2.le.0.01) condition 
   into subr.find_maximal_OX_transmission, as in 
      if( (xe.gt.(x0-eps_xe) .or. (cn2-cnpar2.le.0.01))
         .and. (xe.le.x0)                ) then
   So, the search of the O-X conversion could be triggered
   by value of Nperp becoming too low ~0.1,
   which helps to detect the proximity of evanescent layer.
   YuP[2021-02-04]
   
[224] Improved graphics in PGplot-generated *.ps file, 
   see subroutine plot_fcefuh. YuP[2021-01]

[223] Adjusted subroutine find_rho_x, to avoid infinite loop
   that happens at certain conditions in OXB runs,
   particularly when condition Xe>1 is not met in plasma.
   YuP[2021-01-27]

[222] version="genray_v10.15_201206"

[222] Added printout of |Eplus|/E, |Eminus|/E, |Ez|/E polarizations
   at ray starting point to screen 
   (or log file). BH,YuP[2020-11-30]

[221] Adjusted logic in searching for the proper root in function hotnp 
    (used for id=6), which sometimes was finding a wrong root of the
    dispersion in near vacuum (when two cold roots are very close).
    We appreciate documentation on this issue by Edmund Simpson,
    King\'s College, London,  BH,YuP[2020-11-25] 

[220] pgconst.i was added to the distribution, to make variables used
[220] by the pgplot plot library subroutines explicitly REAL*4.
[220] EG: PGSVP(.2,.8,.3,.95) is  converted to PGSVP(R4P2,R4P8,R4P5,R4P95).
[220] This is prepatory to converting all REAL to explicit REAL*4, so
[220] that -fdefault-real-8 can be applied during compilation, to 
[220] increase accuracy and consistency of the output.  BH201124

[219] Several adjustments are made in the code to avoid over-flow errors
   AYP, YuP [2020-11-25]. 
   In grill_lh, the renormalization powers(i)=powers(i)*1.d13 
   is removed because in case of i_ox=1 this subroutine
   can be called repeatedly. Instead of this renormalization,
   the factor 1.d13 (conversion to erg/sec) 
   is used where needed [2020-11-25].
   
[218] version="genray_v10.14_200831"
[218] The file genray_help has been added to the distribution,
[218] and will subsequently be maintained as the primary documentation
[218] of namelist variable definitions to be used in genray.in and
[218] and genray.dat namelist input files. [BH200914]

[218] New namelist variable:  istep_in_lcfs
   For istart=1 (EC cone) starting method, choose
   how the EC ray is launched with respect to the LCFS:
   istep_in_lcfs=1 means Step inside LCFS (slightly).
   This is the default option, for reverse consistency with original coding.
   istep_in_lcfs=0 means Start ray-tracing directly from {rst,zst} even
   if it is far outside  of the LCFS. That is, no stepping inside LCFS.
   This is more physical, because there is
   some plasma outside of LCFS (exponential drop-off of density),
   and this option is required in cases when O-X conversion
   happens just outside of LCFS.   !YuP[2020-09-03]

[217] More adjustments in subroutines in forest.f,
   related to nperp-->0 limit, which is allowed now.
   Need to use asymptotic expansion in expressions like In(x)/x
   when x-->0 (where In(x) is the mod-ed Bessel function).
   YuP[2020-09-02]

[216] Imposed a limit for Vgroup/c components.
   It gives a better stability for integration.
   Note: although each component of Vgroup/c is limited now,
   the magnitude of |Vgroup/c| can still be larger than 1.0
   (an unphysical condition which can occur with EBW ray tracing).
   YuP[2020-08-30]

[215] Revised subroutine solvnperp for better convergence
   and a more clear logic flow; added many comments.
   Also, now the value of accurcy0 is set to a fixed value,
   independent from prmt4. See comments in that subroutine
   around accurcy0=1.d-4 line.  YuP[2020-08-26]


[214] Added new option for ray integration: irkmeth=3.
   It is faster and more stable than irkmeth=2.
   See read_write_genray_input.f for description.
   YuP[2020-08-26]

[213] Merged improvements from GENRAY-C version into this version.
   Mostly in subroutines in forest.f. Now n_par=0 is allowed
   in many subroutines related to hot plasma.
   YuP[2020-08-24]
   
[212] Modified subroutine efield1, similarly to that in GENRAY-C.
   The changes are related to issue of choosing a proper
   row (one of three) in matrix 
   ((-D_tensor)).(E_vector)=0, [see Eqn(1.6) in Genray manual].
   Usually a row with largest element is selected [see Eq.(8.3-8.6)],
   but there can be cases when there are two rows with nearly same 
   elements (by magnitude). See comments in subroutine efield1.
   YuP[2020-08-21]
   
    
[211] Added checks of grad(psi)=0 (may happen outside of LCFS).
  YuP[2020-08-20]
  

[210] Corrected error in procedure for raypatt='diskdisk'.
   The rays were guided from bin centers of the launching disk
   to bin boundaries of the focusing disk.
   Now they are guided from bin centers to bin centers. 
   !YuP[2020-08-07]

[209] Corrected bug in subr.disk_beam_rays_initial_launching_data.
   there was a typo in 1st argument. !YuP[2020-08-06]

[208] In subr. effcurb(), added check of prho>0 in 
      if(prho.gt.0.d0)then   
        ctheta=(r-xma)/prho
        stheta=(z-yma)/prho
      else
        ctheta=1.d0
        stheta=0.d0
      endif
   !YuP[2020-08-04]

[207] In subr. effcurb() and eff_Lin_Liu(),
   adjusted the definition of lh.
   It was lh=int(sqrt(1d0-cnpar**2)/abs(ye))+1
   When cnpar>1, it results in 
   large negative values for lh, which gives INF 
   in other subroutines.  
   Adjusted to sqrt(max(0.d0, 1d0-cnpar**2)), i.e.
    lh= int( sqrt(max(0.d0,1d0-cnpar**2))/abs(ye) ) +1
   !YuP[2020-08-04] 

[206] In subr.subroutine transmit_coef_ox:
   Added check of grad_x>0. In some cases (outside of LCFS)
   it can become 0, resulting in L_n=x(z,r,phi,1)/grad_x = INF.
   If it happens, set L_n to a large value, 
   L_n= 100.d0 ![cm] any large number - 
   !  it will yield small efficiency in 
   !  transmission_ox_p=dexp(-pi*omega*L_n/clight*...)
   !YuP[2020-08-04] 

[205] In subr. npernpar:
   adjustment for a very small neg.discriminant.
   !YuP[2020-08-04] 


[204] In subr. owconvr, adjusted 
      rhoright=1.d0 !YuP[2020-07-29] This may fail: If edge density is high
                    !Xe=1 layer could be just outside of rho=1.0
   to
      rhoright=1.5d0 !YuP[2020-07-29] Give extra room to find Xe=1.

[203] Noticed a problem in units, in subr.gr_OX_optimal_direction, in
   " if(dabs(f_left*rho_st*1.d+2).le.eps_antenna) then ...".
   (There are three lines like this.)
   Both rho_st and eps_antenna are in meters,
   so why are we using 1.d2 factor? Looks like units-issue.
   Not corrected yet - it may cause a slightly different result
   for old tests (like test4).  Not really critical.
   !YuP[2020-07-29] 

[202] Small "safety" modification in subr.antenna_vertex,
   in definition of 
   gamma= dacos((r_st**2+delta2-r_0**2)/(2.d0*r_st*delta)),
   to avoid case of delta=0.
   !YuP[2020-07-29]
   
[201] Significant change for the i_ox=1 option 
   (a search of optimal launching angle).
   Changes are mostly done in subr. bound, boundc.
   First, for i_ox=1, skip sections related to no_reflection=1 flag.
   It was giving a false reflection/stopping of ray during i_ox=1
   procedure. Now the procedure does not depend on no_reflection. 
   Second, during search, check rho at ray element: 
   is it larger than rho_plasm_vac ?
   rho_plasm_vac is rho at plasma-vacuum border;
   normally, it would be at rho=1.0, but if edge density is high,
   Xe=1 conversion layer could be just outside of rho=1.0
   Then, we define rho_plasm_vac as 
        if(rhoconv.lt. 1.d0)then
          rho_plasm_vac=1.d0
        else ! rhoconv=1.d0 or larger
          rho_plasm_vac= rhoconv+0.1d0
        endif
   (or could be slightly different, see sub.bound)
   where rhoconv is rho at Xe=1.0 conversion layer,
   calculated in sub.ox_conversion_grill_in_poloidal_point.
   !YuP[2020-07-29]

[200] Small modification in subr.cone_ec. 
   Compare average tor. aiming angle (alpha_avg) 
   over rays on the cone with tor. aiming angle alphaj(1) 
   of central ray. If they are opposite, redefine tor. aiming 
   angle of the central ray:  alphaj(1)= alphaj(1) -2.d0*pi
   It is not important for computations,
   but it is good to do this for saving data --> plots.
   !YuP[2020-07-24] 

[199] Made modifications in subr.dinit_1ray. Skip using subr.plasmray 
   and edgecor, which move the ray starting pooint onto LCFS.
   Just start the ray in ~vacuum or wherever it is set in namelist.
   (Could be made into namelist, as an option?) 
   This modification is important for modern tokamaks with high
   edge density, when the O-X conversion may happen outside of LCFS,
   so we should not start rays at exactly the LCFS position.
   A related modification is made 
   in subr.ox_conversion_grill_in_poloidal_point;
   set wide limits for searching O-X conversion point:
   rmn=1.d-5 
   rmx=min(xeqmax,abs(xeqmin)) !basically, equilibrium grid limits.
   !YuP[2020-07-23]

[198] Added more checks in subr. bound and boundc, edgcor. 
   Check if the ray is out of (R,Z)-equilibrium grid. YuP[2020-07-23]

c[197] Corrected call of ccast. Second argument should be (0.d0,0.d0), 
c[197] it was 'zero'. Because of that, some of complex arrays were not 
c[197] properly initialzed, and appeared as very small values
c[197] or QNAN in genray.nc file. 
c[197] This has no effect on physical results of calculations,
c[197] only on values of arrays in *.nc file that suppose to be 0.  
c[197] YuP[2020-03-10]

c[196] Added message+stopping if no eqdsk file was found.

c[196] In param.i, tried to increase nlimit to 501 (was 101) - max number
c[196] of points for LCFS arrays. Also nves - to 201 (was 62).
c[196] However, changing nlimit modified the results (e.g., in test7),
c[196] just a little,
c[196] probably because of accuracy for tracing surfaces in subr.gr2new,
c[196] which uses nteta=nlimit-1.
c[196] Need to keep it in mind when performing auto-tests 
c[196] (verification against data in genray.nc). Changed nlimit back to 101.

c[196] Adjusted logic in sub.equilib for reading 
c[196] (rlimit(i),zlimit(i),i=1,nnlim) and (rves(i),zves(i),i=1,nnves)
c[196] arrays.

c[196] Adjusted logic in subr.calc_r_z_psi_theta_binary  for calculation
c[196] of flux surfaces. If too many iterations, the subr. will stop the run
c[196] (instead of going into INF loop, before this adjustment) 
c[196] with a message 'Probably errors in eqdsk.  STOPPING'.
c[196] See YuP[2020-03-09]

========================================================================
c[195] version="genray_v10.13_200117"
c[195] Converted some of "implicit integer ..." to "implicit none" with 
c[195] explicit declaration of variables. Still many to go...
c[195] See YuP[2020-01], YuP[2020-01-14]

c[194] Adjusted fourb.i, oxb.i, transport_prof.i, writencdf.i --
c[194] Collected same-type variables into separate blocks,
c[194] for better alignment in common blocks. 
c[194] Also, in some of *.i files, changed "double complex" to "complex*16"
c[194] (There are still many "double complex" in *.f files, but not in *.i).
c[194] See YuP[2020-01]

c[193] Adjusted mpi.ins: For wk_pwr_cur() array, changed real to real*4.
c[193] The default size, for a declaration such as "real val", 
c[193] can be altered by compiling with options -dbl or -r8, 
c[193] and we better have it always real*4, for 
c[193] MPI_SEND(wk_pwr_cur,***, MPI_REAL, ***)
c[193] See YuP[2020-01]

c[193.1] Added coll_mult multiplier to Bonoli LH collisional damping.
c[193.1] Paul Parks reexamined coll damping and says it should be about
c[193.1] half to eighth of what Paul Bonoli obtained. 
c[193.1] Needs further checking out.  [BH191017]

c[192] version="genray_v10.12_180912"
c[192] Fixed bugs related to MPI implementation 
c[192] of ADJ-related subroutines (i_adj=1).
c[192] The bugs resulted in zero current. 
c[192] With this fix, the file 'adjout' (usually a large file)
c[192] is only recordered by one core (myrank=0), and 
c[192] read also by one core, and then MPI-broadcasted to other cores.
c[192] See "write (iout5,...)" and "read (iout5,...)".
c[192] YuP[2018-09-12].

c[191] Fixed an argument list dimensioning in call to
c[191] lh_bonoli_disp.f, which made no change in results for test7.1
c[191] of this dispersion relation.  [BH180620]

c[190] version="genray_v10.12_180529"
c[190] Reversed changes made in [187]. 
In the manual, the two cold plasma roots N**2 
(for dispersion id=2) are described as
N2p= (-B +sqrt(B^2-4AC))/(2A)      (4.12a) 
N2m= (-B -sqrt(B^2-4AC))/(2A)      (4.12b) 
In the code, these equations are modified by
multiplying each of A,B,C coeffs by a resonance 
delta=(1-Y) factor (where Y=omega_c/omega)
which cancels out with corresponding 1/(1-Y), so that
a=A*delta, b=B*delta, c=C*delta do not contain 
the diverging 1/(1-Y) factor (either at ECR or ICR).
Then the equations (4.12) should cast into
N2p= (-b +sign_del*sqrt(b^2-4ac))/(2a)      
N2m= (-b -sign_del*sqrt(b^2-4ac))/(2a) 
or, in general form,  
N^2= (-b +ioxm*sign_del*sqrt(b^2-4ac))/(2a) 
where sign_del= sign of (1-Y), and ioxm=+1 or -1.
The sign of (1-Y) changes across the corresponding resonance
(across ECR when ib=1, or across ICR 
when ib=2 or some other value>1).
In the code, the sign_del factor is omitted, 
which looks like an error.
However, a detailed analysis shows that 
when sign_del is included,
the solution (root) jumps from one branch to another
after crossing the corresponding Y=1 resonance.
So, it is better to keep it in the original form:
N^2= (-b +ioxm*sqrt(b^2-4ac))/(2a)
At the same time, it is important to keep in mind that
the selection between two branches for N^2 depends on both 
the value of ioxm and the value of ib.
In general, for ECR frequency range, we should set ib=1.
Then, ioxm=-1 selects the X-mode, and ioxm=+1 selects the O-mode.
[But notice that if you keep ib=1 and go to the ICR range,
there is a jump (switching) between ioxm=+1 
and -1 branches across the ICR layer.]
For the ICR region, we should set ib=2 
(or corresponding ion species in case of nbulk>2).
Then, ioxm=-1 selects the Fast wave, 
and ioxm=+1 selects the Slow wave.
[But notice that if you keep ib=2 and go to ECR range,
there is a jump (switching) between ioxm=+1 
and -1 branches across the ECR layer.]
The most difficult is the case of omega_ci<omega<omega_ce.
Both ib=1 and ib=2 could be selected, 
if none of EC or IC resonances are accessible.
However, depending on ib value, the ioxm value 
corresponds to different branches.
So, in applications like LH, Helicon or whistler waves
it is advised to try different ib values,
and try both ioxm=-1 and +1 cases,
then compare the initial Nperp^2 values.
Example of wave launch in LH frequency range,
at omega_ce/omega ~ 20-30,
omega_ci/omega ~ 0.005-0.008,
and (omega_pe/omega)^2 ~ 100-700:
The LH (Slow wave with large nperp~20-40) can be launched with
  {ib=2, ioxm_n_npar=+1} 
[ioxm value is not important when ioxm_n_npar is not 0] 
or 
  {ib=2, ioxm=-1}  [ioxm_n_npar should be set to 0].
The Fast wave (smaller nperp~5-10) can be launched with
  {ib=2, ioxm_n_npar=-1} [and ioxm value is not important] or 
  {ib=1, ioxm=+1}  [ioxm_n_npar should be set to 0].

c[189] Changed the STOP command in subr. GRPDE2, 
c[189] when 'nper2 is too negative' (related to calc. 
c[189] of group velocity) to setting vgrpdc=0.d0.
c[189] It will stop the ray, but would not halt the run 
c[189] as with the STOP command. This is particularly important 
c[189] in MPI runs, when a STOP generated at one CPU core
c[189] effectively hangs the run - other cores 
c[189] (having no STOP for their rays) are waiting 
c[189] until the wall time runs out.
c[189] Also added ibmx=min(ib,nbulk) - 
c[189] a safety check, not to exceed nbulk.
c[189] YuP[2018-05-23]


c[188] version="genray_v10.10_180518"

c[188] Added a namelist variables for control/suppressing of printout,
c[188] saving data and saving plots:
c[188]  outnetcdf (to suppress writing data to *.nc)
c[188]  outprint (to suppress printing to the screen)
c[188]  outxdraw (to suppress saving files for xdraw).
c[188] Many lines with printout are permanently commented out,
c[188] which were used for debugging printout. For example, 
c[188] see lines that start with 'cyup' in dinit.f. 
c[188] Also added more detailed printout of cpu time spent
c[188] in different parts of the code, like cpu_time(time_drkgs2_2).
c[188] YuP[2018-01-17] 


c[187] A problem noticed and fixed by YuP[07-2017] (later adjusted).
c[187] In equation for the two roots [Eqn. numbers refer to Genray manual],    
c[187]	    N2p=(-B +sqrt(B^2-4AC))/(2A)      (4.12a) (O)
c[187]	    N2m=(-B -sqrt(B^2-4AC))/(2A)      (4.12b) (X)
c[187] the sign in front of sqrt() [defined as ioxm] determines the mode: 
c[187] '+' is for O-mode, '-' is for X-mode .
c[187] However, after A,B,C are mult-ed by delib=(1-Yib) factor,
c[187] to eliminate the resonance denominator in eps1 and eps2,
c[187] the result depends on the sign of delib=(1-Yib).
c[187] If (1-Yib) is positive, nothing changes 
c[187] in correspondence of ioxm=+/-1 and the two modes. 
c[187] But for the negative (1-Yib), e.g. (1-Yib)=-1, we get
c[187]    (-b +ioxm*sqrt(b^2-4ac))/(2a) = [use a=(1-Yib)*A= -A, etc] 
c[187]  = (+B +ioxm*sqrt(B^2-4AC))/(-2A)= 
c[187]  = (-B -ioxm*sqrt(B^2-4AC))/(+2A)
c[187] Thus, for ioxm=+1, we are getting the branch (4.12b), 
c[187] which is the X mode. THE MEANING OF MODES IS REVERSED!
c[187] To correct this problem, we simply need to further adjust 
c[187] the a,b,c cofficients:
c[187]    sign_delib=sign(1.d0,delib) 
c[187]    a=a*sign_delib 
c[187]    b=b*sign_delib 
c[187]    c=c*sign_delib 
c[187] So, effectively, we use a=|1-Yib|*A, etc.
c[187] Then, the meaning (mode type) defined through a,b,c
c[187] will remain the same as that defined through the original A,B,C.
c[187] This is done in many subroutines of cninit.f file,
c[187] also in subr. abc, and also in dddrz1, rside1 and others where 
c[187] the analytical derivatives of the dispersion equation (id=2)
c[187] are setup. Now the LHCD type of runs DO NOT depend on value of
c[187] ib variable (which corresponds to the species number   
c[187] that can have ECR or ICR in plasma, 
c[187] and which determines the (1-Yib) factor).
c[187] Besides, the search of a launching wave type is more stable now.
c[187] For the Fast wave / Slow wave selection, when the initial 
c[187] value of Npar is given, it is recommended 
c[187] to use ioxm_n_npar option instead of ioxm.
c[187] Use ioxm_n_npar=-1 for FW, and +1 for SW.
c[187] Search "YuP[07-2017]" for all corresponding changes. 



c[186] version="genray_v10.10_170301"
c[186] Small bugs fixed.  Netcdf output of w_dens_vs_r_nc units fixed.
c[186] [BH and YuP, 170720].

c[185] version="genray_v10.10_151110.2"
c[185] Bug fix in prep3d.f, enabling i_ox=1 (optimized OX mode
c[185] conversion, test4 to work). [BH160227].

c[184] version="genray_v10.10_151110.1"
c[184] Another bug fix for i_grill_npar_ntor_npol_mesh.eq.1 option
c[184] for launching of rays from inside the plasma on a grill of
c[184] ntor and npol values [See cShiraiwa151113].

c[183] version="genray_v10.10_151110"
c[183] Bug fix in grill_lh.f for i_grill_npar_ntor_npol_mesh.eq.1 option
c[183] for launching of rays from inside the plasma on a grill of
c[183] ntor and npol values [Shiraiwa151107, BH151110].

c[182] version="genray_v10.10_151015"
c[182] The EC current drive harmonic number may now be specified
c[182] automatically through formula, origially given in TORBEAM.
c[182] by specifying jwave=0 when ieffic=3 or 4.. See notes in prep3d.f.
c[182] Addition proposed by Nicola Bertelli.  
c[182] A bug in reading temperature profiles outside the LCFS when
c[182] i_edge_dens_rz_mesh=2 was fixed, as proposed by Syunichi
c[182] Shiraiwa  [BH151018].

c[181] version="genray_v10.9_150615"
c[181] YuP [06,2015] In files prep3d.f, lsc_approach.f, lh_ql_flux.f
c[181] Added a check for rho_loc=0.d0 condition (can happen if a ray
c[181] goes through the magnetic axis).
c[181] If rho_loc=0.d0 is detected, the values of cos and sin 
c[181] of pol.angle are set to 
c[181] cos_theta_pol=1.d0 and sin_theta_pol=0.d0

c[180] YuP[06,2015] In subroutine read_all_namelists :
c[180] (In files read_write_genray_input.f  
c[180]     and   read_write_genray_input_prep.f )
c[180] Setting i_adj to 0 for all ieffic that do not need it.
c[180] Otherwise some subroutines are called that are related to ADJ,
c[180] and it may result in undefined values.
c[180] So, now i_adj=1 is allowed only for (ieffic=5).or.(ieffic=6).

c[179] YuP [June 2015] In subroutine CD_adj_LH_efficiency_1 :
c[179] In some cases, b_ratio=0 (print-out: at rho_small=0)
c[179] Skip calc. of current drive in such case.


c[178] version="genray_v10.8_141024"
c[178] MPI version: modified mpi.ins file (in /mpi/ folder).  Changed
c[178]  real*4 wk_pwr_cur(5+NR*(6+nbulk))  to
c[178]  real wk_pwr_cur(5+NR*(6+nbulk))
c[178] and changed  MPI_FLOAT to MPI_REAL .
c[178] Before that, an error was triggered when using PathScale (PPPL)
c[178] (MPI_ERR_TYPE: invalid datatype).
c[178] Probably, in PathScale,
c[178] real*4 and float (MPI_FLOAT) are not the same. [YuP141016].

c[177] MPI version is updated - changes in genray.f, within do 21 
c[177] loop over rays. The bug was related to logics/control of 
c[177] SEND/RECV communication.  For "short" rays (bad initial 
c[177] condition or near zero initial power) the SEND and RECV
c[177] commands, and procedure of confirmation of receiving the data
c[177] were logically disconnected.   [YuP141007]. 

c[176] version="genray_v10.8_140909"
c[176] Added immediate stopping of a ray, if initial power is zero
c[176] [.lt.1.e-100], to save compute time. iray_status_nc=14. 
c[176] [BH140909].

c[175] version="genray_v10.7_140801"
c[175] Added momentum conservation in the Lin-Liu subroutine
c[175] (lin-liu_curnt.f), controlled by new namelist variable 
c[175] ieffic_mom_cons (which works with ieffic=4).
c[175] The implementation is the same as the implementation done
c[175] in the TORBEAM code. More specifically, the model is described
c[175] in Marushchenko et al, PoP 18, 032501 (2011) and Marushchenko
c[175] et al, NF 48, 054002 (2008). [Nicola Bertelli, 140716]

c[174] Fixed write of lsc-approach namelist, per Andre and F. Poli.
c[174] Made adjustments to template file, for clarity.
c[174] Added three test cases (test8.1-3) comparing various methods for
c[174] calculation of O-mode ITER ECCD. New  ieffic=4,ieffic_mom_cons=1
c[174] compares well with ieffic=6/i_adj ADJ calcs  and CQL3D calcs.
c[174] Adjusted makefiles regarding new ieffic-mon-cons .f90 files.
c[174] Zeroed out remainder of delpwr_lsc, removing previous iteration 
c[174] data in lsc_approach.f:subroutine LSC_power_iterations.
c[174] [BH140729]

c[173] version="genray_v10.6_140421"
c[173] Several small bugs which affected execution/use at PPPL were
c[173] fixed according to Rob Andre presecriptions [BH140421].

c[172] version="genray_v10.5_131024"
c[172] Modified character length of 512 character variables down
c[172] to 256, in response to a compiler limitation [BH131024].

c[171] Fixed a bug in the outpt ray testing subroutine which enabled
c[171] ray lengths to greatly exceed the nominal max, nrayelt
c[171] [BH, 131024].

c[170] Added delim='apostrophe', in several open() of namelist files, 
c[170] to ensure that any namelist writes include apostrophe delimiters
c[170] in the character variables.  Default for PGI and Intel compilers
c[170] is delim='none'  [BH 130911].

c[169] version="genray_v10.4_130726",  YuP+BH re-implemented MPI.
c[169] The MPI feature had been broken during introduction of 
c[169] dynamic dimensioning, item [162] and ongoing re-dimensioning.
c[169] [YuP and BH, to July 27, 2013.]

c[168] version= "genray_v10.2_120629"
c[169] Added Bonoli-type scattering of LH ray (lh_scattering.f),
c[169]  controlled by new namelist variable iscat_lh_nicola
c[169]  [Nicola Bertelli+SAP+BH, 120529].

c[167] version="genray_v10.1_120424"
c[167] Added collisional power absorption to the mnemonic.nc[BH].
c[167] Fixed several bugs:cBH120421 in forest.f; 
c[167]       cYuP+BH120424 in oxb.f and dinit.f.

c[166] version="genray_v10.1_120308"
c[166] Added plotting of dispersion D(ReN_perp,ImN_perp) contours
c[166] At given RZ N_parallel points before and code stops 
c[166] before ray calculations [SAP, 120308].

c[165] Add LH dispersion relation including Bonoli/Englade warm plasma
c[165] term, id=16 [SAP, 120114].

c[164] Added collisional power absorption to the mnemonic.nc output
c[164] file.  To be used with SWIM IPS.  [BH, 111212].

c[163] Code version="genray_v10.0_111101".
c[163] Added ilaunch=2 grill ray launch option, to launch rays at a
c[163] list of R,Z points in the poloidal plane.  This will enable
c[163] specific modeling of a physical tokamak antenna setup. 
c[163] [SAP and BH, 110331].

c[162] Code version="genray_v9.1_w_lsc-option_111011".
c[162] Adjusted sytem for input of density/temp (per c[160])
c[162] to reading fixed format files.  Made fine RZ density/temp
c[162] grid dimensioning dynamic.     See template file
c[162] genray.in_template_MKSA_111011 re i_edge_dens_rz_mesh.   
c[162] Also, moved coding from read_write_genray_input.f related to
c[162] to i_edge_dens_rz_mesh to equilib.f, fixing it so
c[162] it can be readily used with the SWIM IPS software system.
c[162] [SAP with BH, 111011 and 111025].

c[161] Code version="genray_v9.0_w_lsc-option_110801"
c[161] Fixed bug which added a radial oscillation to power and current 
c[161] drive radial profiles under some conditions of ray step size.
c[161] [SAP, 110503].

c[160] Density and temperature outside the LCFS can now be entered 
c[160] via an external file.  This data is bi-cubic-splined in genray
c[160] and used in the ray tracing.  A file can be output of density
c[160] and temperature throughout the RZ plane, this file can be
c[160] modified for use outside the LCFS, and read back in for further
c[160] See namelist variable i_edge_dens_rz_mesh [SAP, ~110503 and 110314].

c[159] Substantial changes in edge density outside LCFS, in dense.f.
c[159] Need to discuss with SAP.  Code version 110503, compared to
c[159] 110314.

c[158] Fixed bugs in cone_ec.f.  Ray pattern width function had
c[158] a square left out, widening the pattern somewhat. [TJenkins, 110419].
c[158] real*8 n_r,n_phi,n_z, na1r, na2r added, with discernible effects.
c[158] [YuP and BH, 110419].

c[157] Fixed bug in dimensioning for new wall reflection in chamber_wall.f.
c[157] Had effect when nyeqd_add .ne. nxeqd_add [SAP, 110403].

c[156] Code version="genray_v9.0_w_lsc-option_110314"
c[156] Corrected extra calling arguments of several subroutines,
c[156] with no effect of regression test suite (10 at this time)
c[156] [SAP110323].

c[155] Major addition to code: LSC 1D-in-vpar FP determination
c[155] of LHCD, following the approach in D.W. Ignat, E.J. Valeo,
c[155] and S.C. Jardin, Nucl. Fusion 34, 837 (1994).

c[154] Added subroutine r4bcast(a,val,n) for real*4 [SAP]

c[153] Added array iray_status giving status of ray ending
c[153] iray_status wil be print out at the screen
c[153] and to the output netcdf file SAP100517

c[152] Corrected times tau_e,tau_n,tau_ei in lh_ql_flux.f file [SAP].

c[151] version="genray_v8_3_20100507"
c[151] Added defaults for all namelist, so can run with empty namelist
c[151] input file.  Added screen printout of all parameters.

c[150] version="genray_v8_2_20100105"
c[150] Calculate (r,z) points at LCFS and 
c[150] arrays rpsi(j,i),zpsi(j,i) at given (psi,theta)
c[150] using new subroutine based on binary method.
c[150] With these subroutines the code can work at psifact
c[150] =0.999 and =1.00 [SAP091218].

c[149] Added new input variable shift_rho_geom
c[149] for itart=1 case
c[149] It shifts the initial point of the ray from the plasma
c[149] bounday (LCFS) inside the plasma.
c[149] SAP091127

c[148] Calculate arrays rpsi(j,i),zpsi(j,i) at given (psi,theta)
c[148] using field lines. It works for all nnlim values
c[148] for eqdsk with and without limiter points.
c[148] SAP091116  Now this part of the code is switched off.
c[148] rpsi(j,i),zpsi(j,i) are calculated using binary solver

c[147] Added new subroutines (in the file limit.f) which work 
c147[] for eqdsk without limiter points, at  nnlim=0.
c[147] Subroutines in limiter.f file determine the limiter points
c[147] from the solution of the field line equation
c[147] Using these subroutines the code works at psifactr values 
c[147] close to one [SAP091104].

c[146] Made corrections in equilib.f file in subroutine limitr()
c[146] for calculation of the monotonic plasma boundaries:
c[146] zzp() from rrp() and zzp() from rrm() [SAP091104].

c[145] The current and absorbed power profile are calculated if 
c[145] the input parameter onetwo=1. If ionetwo=0,then some 
c[145] pointers will not be  allocated (which has zero values for
c[145] these profiles), and call dnonetwo will give segmentation
c[145] fault.  It was added to genray.f 
c[145] if(ionetwo.eq.1) then call dnonetwo [SAP091031].

c[144] If partner is used than the current and power profiles      
c[144] are used.   For partner calculations it should be ionetwo=1.
c[144] It was added to read_write_genray_input.f:
c[144] if (partner.ne.'disabled') ionetwo=1
c[144] SAP091031

c[143] Made corrections in grill_lh.f inside subroutine rho_ini_LHFW
c[143] In the case when the ray was launched outside the plasma
c[143] (outside LCFS at rho_loc.gt.1.d0), the code calculated 
c[143] space coordinates(r,z) from the given (psi,theta_poloidal)
c[143] using subroutine  zr_psith(psi,theta,z,r) based on spline. 
c[143] Spline coefficients in zr_psith(psi,theta,z,r) were created 
c[143] at the uniform mesh psimag < arpsi(1:npsi) < psilim.
c[143] At the point outside LCFS (psi > psilim) spline works like 
c[143] the extrapolation procedure.
c[143] Now if the point is outside LCFS subroutine subroutine
c[143] rho_ini_LHFW uses the geometrical small radius 
c[143] rho_loc_geom=sqrt((z-yma**2)+(r-xma)**2)) 
c[143] to find (z,r) coordinatres 
c[143] z=yma+rho_loc_geom*dsin(theta)
c[143] r=xma+rho_loc_geom*dcos(theta).
c[143] Same coordinates are used for density and temperature profiles
c143[] outside LCFS.
c[143] SAP091027

c[142] Made correction in dinit.f using new local name of argument
c[142] n_toroidal=cnphi in call rho_ini_LHFW
c[142] The problem was that subroutine 
c[142] call rho_ini_LHFW(..,cnphi,..)
c[142] can change the input argument cnphi.
c[142] Now the code uses call rho_ini_LHFW(..,n_toroidal,..)
c[142] SAP091027


c[141] Made some corrections for work with the grill conditions
c[141] at i_n_poloidal=3 case
c[141] The corrections are in cninit.f.
c[141] SAP091027

c[140] version="genray_v8_1_091021"
c[140] Fixed bug in write of power density profiles, powden_s, for
c[140] indivual ions to netcdf output file [SAP, 091021].

c[139] Added ability for code to continue with zero densities in
c[139] spline input, e.g., near boundary [SAP, 091020].

c[138] Fixed Zeff in nonuniform-mesh input system.  There was bug
c[138] in transfer of input data to internal regular mesh which
c[138] the code uses [SAP, 091016].

c[137] Modified procedure to calc modified Bessel function from separate
c[137] In(x) times exp(-x) calculation to subroutine calculating the
c[137] product directly, removing NaN in absorption.    Increased number 
c[137] of gyroharmonics from 100 to 200, removing a problem that stopped 
c[137] the code, but probably causes no effective change in results
c[137] [SAP, 091016].

c[136] Added pointers for 
c[136] adj_no_mnl.i,
c[136] dskin.i
c[136] emissa_no_nml.i
c[136] onetwo_no_nml.i,
c[136] six_no_nml.i,
c[136] spline_distrib.i
c[136] write.i
c[136] writencdf.i
c[136] [SAP090831].

c[135] Modified ion edge density such that namelist dens_min_edge gives
c[135] that ion density ratios at the LCFS are maintained throughout the
c[135] SOL edge region, resulting in zeff constant in the SOL[SAP090527].

c[134] mk_graph.f has new subroutine which plots to plot.ps file
c[134] plots of the the cold plasma dispersion relation  roots at the 
c[134] initial ray point (z,r,phi) with given N_parallel.
c[134] The dispersion relation has the following form:
c[134] eps*N_perp**4-a2*N_perp**2+a0=0
c[134] The roots are
c[134] N_perp_p**2 (frequency) = (a2+cqrt(a2**2-eps*a0))/(2*eps)
c[134] N_perp_m**2 (frequency = (a2 -cqrt(a2**2-eps*a0))/(2*eps)
c[134] [SAP090523]

c[133] Fix to density outside LCFS [SAP090516].


c[132] Corrected subroutines in cninit.f for i_n_poloidal.eq.3 case.
c[132] Corrected commented goto in dddrz1.f causing problem for
c[132] an id=4 case [SAP090505].

c[131] Corrected condition in subroutine drkgs2 in rk_new.f  
c[131] for small time step
c[131] if (dabs(h).lt.1.d-11) then [SAP090502].

c[130] Added subroutines to create exponential density fall
c[130] near near the vertical limiter boundaries 
c[130] perpendicular to the toroidal direction  [SAP090425].

c[129] Added stopping the ray at the vertical limiter boundary,
c[129] which is perpendicular to the toroidal direction e_phi
c[129] Added xdraw plots of the toroidal limiter boundaries
c[129] at XY plane  [SAP090414]. 

c[128] Added subroutines to create exponential density fall
c[128] near near limiters
c[128] SAP090413

c[127] Added subroutines to create exponential density fall
c[127] near the chamber wall.
c[127] SAP090407

c[126] Version="genray_v7-16_090406". Updating svn and tagged version.
c[126] Fixed bug in screen output accounting for collisional power
c[126] [SAP090406].   

c[125] Fixed bug in calc of drho_dz and drho_dr functions in SOL.
c[125] SAP090331.

c[124] Add refinement of the wall coords, inserting additonal
c[124] wall coords to aid test of rays passing outside the chamber
c[124] wall, using new namelist h_add_wall [SAP090325].

c[123] moved sigmedgn/sigmedgt from &wave to &edge_prof_nml.  Add
c[123] namelist variable temp_min_edge.  SAP090304

c[122] Code version="genray_v7-15_090228"
c[122] Corrected rk_new.f subroutine drkgs2 Runge Kutta integrator.
c[122] If step size  h<1.d-11 then stop trajectory calculations.
c[122] In files dxdz.f dxdy.d x.f files,  corrected the
c[122] condition: if the ray point is outside the LCFS.
c[122] SAP090229

c[121] Added namelist edge_prof_nml to specify plasma density outside
c[121] the LCFS. Added calculations of plasma density outside LCFS
c[121] The dependence of the plasma density on the poloidal angle can
c[121] be set by the table (for spline approximation) or by the given
c[121] analytic formula.  
c[121] Files added: edge_prof.i,edge_prof_nml.i,
c[121]              edge_prof_no_nml.i
c[121] SAP090214

c[120] Added namelist i_chi_interpolation to chose interpolation method
c[120] for calculation of derivatives from the chi function.
c[120] Splines have given some problems near the trapped-passing bndry.
c[120] Default is splines [SAP090130].
c
c[119] Added calculations of power absorption using QL flux.
c[119] Added new input variable iabsorp_ql.  This verifies the QL
c[119] calculation against alternate absorption calculations.
c[119] iabsorp_ql =0  do not use QL flux for absorption calculations
c[119] iabsorp_ql =1  to use QL flux for electron absorption calculat
c[119] SAP081203

c[118] Corrected prep3d.f file at ion_absorption='disabled'
c[118] SAP081111

c[117] Renamed psisep for psilim
c[117] SAP081110

c[116] The code was modified for i_rho_cutoff=1 
c[116] to work in the situation that a ray
c[116] LH or FW cutoff is not found [SAP,081028].

c[115] Added input data to namelist/tokamak/
c[115] to read chambar wall and limiter points
c[115] n_wall,max_limiters,n_limimer,r_wall,z_wall,r_limiter,z_limiter,
c[115] phi_limiter.
c[115] Added subroutine to write these input data and eqdsk poloidal flux
c[115] to genray.nc file
c[115] [SAP081024]

c[114] Added -i into reps(1,2)=-ci*xe/ye for iabsorp=3 (Chiu polarizations)
c[114] bringing iabsorp HHFW polarizations into approx agreement with
c[114] iabsorp=4 (Stix) for an NSTX case. [SAP, 080916]
c
c[113] Add poloidal angle of ray on flux surface to cd.bin and genray.nc
c[113] output files.
c[112] Code at version="genray_v7-10_080906".  [SAP, 080906] 
c
c[112] Changed procedure of the control in b.f file that the point is
c[112] inside the last closed flux surface.
c[112] Made new function rho_zr(z,r) which uses the given control procedure. 
c[112] Code at version="genray_v7-9_080902".  [SAP,080817]

c[111] Added namelist variable prmt9 giving accuracy of 
c[111] Hamiltonian in the Runge-Kutta subroutine for irkmeth=2.
c[111] This will reduce Runge-Kutta time step if dabs(ham).ge.eps_ham)
c[111]       =1 by default
c[111] [SAP080805]

c[110] Added namelist variables sigmedgn and sigmedgt, controlling
c[110] edge profiles for no_reflection=1 case [BH080731].
c[110] Code at version="genray_v7-8_080801". 
c
c[109] CD efficiency at rho.gt.1 set to zero.
c[109] Errors in calc of power profile for collisions fixed (prep3d.f).
c[109] Error in transformation from genray.in to genray.dat format.
c[109] [SAP, 080731]

c[108] Add input variablel:
c[108]  no_reflection !=1 switch off the artifitial reflection from 
c[108                 !the latest closed flux surface
c[108]                !=0 (default) switch on the artifitial reflection from 
c[108]                !the latest closed flux surface
c[108][SAP080729]

c[107] Modified read_write_genray_input.f for usability with Plasma
c[107] State, giving version="genray_v7-6_080714" [SAP and BH, 080714].
c
c[106] Added partner='genray_profs_out.nc' namelist option, which
c[106] obtains all plasma profiles from the genray.dat or genray.in file
c[106] but writes genray output data to genray_profs_out.nc.
c[106] giving version="genray_v7-7_080707"[BH, 080707].
c
c[105] Designating this version="genray_v7-3_080609".
c[105] Add input variable:
c[105] i_calculate_or_read_adj_function  =1 (by default) calculate the 
c[105]                                       adjoint function chi
c[105]                                   =0 read chi function from
c[105]                                           the file: adjout
c[105] [SAP080602]
c
c[104] Changed names of some dimensions in the main netcdf output file
c[104] to removed redundancy and increase consistency:
c[104] New nbulkm=previous nispec, old nispecdim now nbulkm,
c[104] old nbulkdim now nbulk, old nrhodim now nrho,
c[104] old nrhodim1 now nrhom [BH, 080605].
c
c[103] Add input variable ion_absorption to switch on/off ion absorption
c[103] at iabsorp=3,9,91,92 cases.
c[102] [SAP080303]

c[102] Change subroutine efield1.f to remove jumps of
c[102] electric field polarization. 
c[102] [SAP080228]

c[101] Designating this version="genray_v7-1_080128".
c[101] Fixed bug in analytic density profiles.
c[101] Add namelist variable refl_loss, giving fractional power lost
c[101] at each ray reflection.  [SAP080128].

c[100] Designating this version="genray_v7-0_071009".
c[100] Write DC electric field conductivity (calculated using ADJ
c[100] function) at npsi0 radial mesh points to output netcdf file.
c[100] This is for a test of ADJ[SAP 071116].

c[95] Write CD efficiency along all rays to output netcdf file
c[95] [SAP 071113.]

c[99] Add for ieffic=6 calculate efficiency using ADJ function from
c[99] Karney code general case (arbitrary form of resonance curve
c[99] for arbitrary N_parallel and EC harmonics)
c[99] [SAP 071025]

c[98] Add for ieffic=5 calculate LH CD efficiency using ADJ
c[98] function from Karney code
c[98] [SAP 070925]

c[97] Add reading input data from the file genray.in in MKSA system
c[97] Add the template file genray.in_template_MKSA_070911
c[97] [SAP 070911]

c[96] Modify output.f to stop cold plasma rays for which
c[96] k_perp*rho_larmor_e becomes .gt. 1 , and in vicinity
c[96] of UHR. Prevents singularity in cold plasma rays
c[96] [SAP and BH, 070808].

c[95] Fixed bug which caused n_theta_pol != 0 to stop code.
c[95] Fixed bug in poloidal current calculation  [SAP, 070807].

c[94] Designating this version="genray_v6-2_070801".
c[94] Fixed bug in electron temperature profiles for nbulk.gt.1,
c[94] which crept into genray_v6-0.

c[92] Modified calc of error parameter in output.f for numerical
c[92] integration of ray eqns to only include refractive indices,
c[92] removing ray position.   This regularized the integration in
c[92] some way [SAP 070728].

c[91] Designating this version="genray_v6-1_070726".
c[91] Add new input variables in namelist /numercl/ (and functionality):
c[91] i_resonance_curve_integration_method to choose integration method
c[91] epsi  is the absolute accuracy used in adaptive Simpson.
c[91] This functionality applies to iabsorp=6,7 and emission calc.
c[91] [SAP 070725]

c[90] Designating this version="genray_v6-0_070712".
c[90] Improvements in printed output to plots of dispersion relation
c[90] through &output namelist parameters (in plot.ps)
c[90] [SAP 070712].

c[89] Add subirectory /read_write_namelists which contains
c[89] files to read and write all namelist data from and
c[89] to input genray.dat file. The subroutine prepare_genray_input from
c[89] read_write_genray_input.f file  reads and writes data
c[89] from and to genray.dat.
c[89] Example program 'test' with call prepare_genray_input is in test.f
c[89] [SAP 070427].

c[88] Add new subroutines into read_write_genray_input.f file: 
c[88] subroutine read_all_namelists reads all input namelists from genray.dat
c[88] subroutine write_all_namelist writes all namelists to genray.dat file
c[88] [SAP 070427]


c[87] Add facility to input plasma profiles at nonuniform small radius meshs.
c[87] These profiles and radial meshes are set by line form rather than 
c[87] column [SAP 070427].

c[86] Include files *.i were divided into two groups of files: 
c[86]     *_nml.i and *_no_nml.i
c[86] *_nml.i group contains declarations and common block *_nml 
c[86] for input namelist\'s data, which were previously used in 
c[86] common block *.i.   This rearrangement of namelist is for the
c[86] SciDAC SWIM IPS (Integrated Plasma Simulator).
c[86] [SAP 070427]

c[85] Add new input namelist to genray.dat file switching the output step 
c[85] size prmt6 near EC resonance points for accurate power absorption
c[85] calculation:  i_power_switch_resonance,n_power_switch_resonance,
c[85] y_power_switch_resonance,del_y_power_switch_resonance
c[85] prmt6_power_switch_resonance.
c[85] [SAP 070419]

c[84] Add adaptive Simpson integration along resonance curve 
c[84] in the velocity space for calculations of anti-hermitian relativistic
c[84] tensor which used for ray absorption. 
c[84] [SAP 070419]

c[83] Corrected small radius rho and derivatives 
c[83] d_rho/d_r,d_rho/d_z, d_density/d_r,  d_density/d_z
c[83] culculations outside the last closed flux surface 
c[83] [SAP 070301]


c[82] Add subroutine refractive_index_relative_error and 
c[82] input parameter to namelist /numercl/ toll_hamilt
c[82] to measure error in the dispersion realation
c[82] If D/(N|gradD|) > toll_hamilt stop ray calculation
c[82] [SAP 070223]

c[81] Add Lin-Liu subroutine to calculate current drive efficiency
c[81] at ieffic=4.Subroutine TorGA_curgap and used subroutines were
c[81] transformed from Fortran 90 to Fortran 77
c[81] [SAP 070221]

c[80] Add: write frequencies fpe,fce,j*fce,fuh,fpi,fci,j*fci,flh
c[80] to frequency1.bin and  frequency_ion.bin files
c[80] Add drawfreq1.in and drawfreqion.in files to plot frequencies
c[80] using xdraw
c[80] [SAP 070217]

c[79] Add output N_perp roots and the polarization
c[79] to the file: hot_roots_at_given_point.dat
c[79] in subroutine calculate_hot_nperp_roots
c[79] [SAP 070214]

c[78] Add calulate output N_perp roots and the polarization along
c[78] the small radius in points, where the hot dispersion function has
c[78] one two or three roots. Write these data to the file: find_hot_roots.dat
c[78] in subroutine calculate_hot_nperp_roots
c[78] New input data in /wave/ of genray.dat:
c[78] i_rho_find_hot_nperp_roots,
c[78] rho_step_find_hot_nperp_roots,rho_min_find_hot_nperp_roots
c[78] [SAP 070213]


c[77] Add new input parameters to genray.dat file
c[77] to plot D(Nperp) and to calculate and choose hot plasma roots:
c[77] i_look_roots,n_nperp_plot,cnperp_plot_min,cnperp_plot_max 
c[77] cN_perp_root_max,n_nperp_plot,k_hot_root 
c[77] [SAP 070207]

c[76] Add hot plasma dispersion function root solver
c[76] [SAP 070207]

c[75] Add PGLPLOT plots of D(Nperp) hor cold and hot plasma
c[75] [SAP 07012]


c[74] Add i_diskf=4 case
c[74] analytic calculation of continuous non-Maxwellian distribution
c[74] with three temperatures in three energy ranges
c[74] [SAP 061227]

c[73] Add calculation of emission spectrums.
c[73] Put emission spectrums to output genray.nc file.
c[73] New input variables in genray.dat file for spectrum calculations:
c[73] i_emission_spectrum,jx_kin,max_kin_energy_kev 
c[73] [SAP 061226]

c[72] Add raypatt='diskbeam' case in genray.dat file 
c[72] [SAP 061209]

c[71] Add the spline coeff1 and terp1 for plasma profiles
c[71] [SAP 061207] 

c[70] Add iabsorp=91 and 92 cases for absorption calculations
c[70] [SAP 061129]

c[69] Add ioxm_n_npar to .nc output file [SAP 061113]


c[68] Add FW absorption options iabsorp=91 [Pinsker and Porkolab
c[68] ions, generalized 05/02/11 email and e-absorption from Chiu 
c[68] et al, NF (1989)], and iabsorp=92 [e and i absorption according 
c[68] to Pinsker and Porlolab email] [SAP, 061127].

c[67] Add namelist input variable nonuniform_profile_mesh
c[67] to use input radial profiles at nonuniform radial grids
c[67] Add namelists to set profiles at nonuniform grids
c[67] dentab_nonuniform, temtab_nonuniform,tpoptab_nonuniform, 
c[67] vflowtab_nonuniform, zeftab_nonuniform [SAP 061104]

c[66] Add namelist input variable rho_step_find_LHFW_cutoff
c[66] to set step in search for initial radius of propagation
c[66] of LH or FW (i_rho_cutoff=1 case). [SAP,BH 060904].
c
c[65] Added input variable ioxm_n_npar to choose mode
c[65] N(N_parallel,ioxm_n_npar) from cold plasma dispersion
c[65] relation [SAP 060925].

c[64] Added: input variables in namelist &output to plot the cold plasma
c[64] dispersion function D_cold(N_perp) and cold plasma wave normal 
c[64] surfaces along the ray at the give arrays
c[64} of poloidal distances or major radii. [SAP 060918].

c[63] Added:  input variable i_salphal(nbulka) in genray.dat file
c[63] It sets which damping will be in salphal_nc in output netcdf file
c[63] [SAP 060915].

c[62] Submitted to CVS at transpgrid.pppl.gov, as version
c[62] genray_v5-0_060804 [BH, 060805].

c[61] Fixed bug in edge treatment of LH in output.f (commented out
c[61] goto120).  Changed parameter name ncompa to nbulka, for 
c[61] consistency [SAP, BH 060804].

c[61] Added namelist coll_mult, multiplier of collisional damping
c[61] [BH, 050715]

c[60] Added plot of contours of D(ReN_perp,Im_N_perp)
c[60] at points with given major radius specified by new namelist
c[60] variables :n_plot_disp,id_plot_disp,r_plot_disp,s_poloid_plot_disp,
c[60] point_plot_disp [SAP060630].

c[59] Added BobH files for partner [SAP060415].

c[58] Added Eric files for iabsorp=12 and id=15 cases
c[58] [SAP060414].

c[57] Added angle integration (like in horace) for integral calculations
c[57] along the resonance curves in emission calculations[SAP060412].
c[57] The code has variable i_resonance_curve_integration_method
c{57] inside emission.f and relat_tens.f
c[57] It chooses integration formulas: =1 for angle integration,
c[57] =2 for p_perp integration with rectangle formula
c[57] =3 for p_perp integration with trapezoidal formula
c[57] Now the code uses
c[57] i_resonance_curve_integration_method=1
c[57] for ellipse case abs(N_par)<1 and
c[57] i_resonance_curve_integration_method=3
c[57] for hyperbola and parabola cases abs(N_par)>=1.

c[56] Added the Bessel function table in emission calc.
c[56] This sped up the emmision calculation by ~4x [SAP060322].


c[55] Changed the boundaries in the integral for resonance curves,
c[55] in the emission calculation[SAP060311]

c[54] Added trapesoidal formula for interal calculations along the
c[54] resonance curves in emission calculations[SAP060308].

c[53] Add namelist control for using a range of harmonic numbers, 
c[53] n_relt_harm1,n_relt_harm in genray dat file [SAP060307 SAP060414]].

c[52] Corrected bug in emission.f: in 'additional points calculation'
c[52] for iabsorp=7 case the electric field calculations did not allow
c[52] for iabsorp =7 case. Fixed it in emission.f file. [SAP060224]

c[51] Corrected bug in prep3d.f calculations of zfacgeom using 
c[51] b_average(psi) instead of using bmin(psi) and aspect ratio,
c[51] for calculation of j_ONETWO [SAP 0602011].


c[50] Corrected call of subroutine write_transport_prof for partner case
c[50] changing the input argument for s_cur_den_onetwo
c[50][ SAP 0602010]

c[49] Add the emission calculations at multiple-rays
c[49] at each frequency. 
c[49] The code writes the data for this case into netcdf file
c[49] [SAP 060209] 

c[48] Add raypatt='diskdisk' case to launch EC rays from disk 
c[48] to disk [SAP 060209]

c[47] Add the simple collisional e-i absorption ~nu_ei/vgroup
c[47] for ray-tracing and for emission calculations [SAP 060119].

c[46] Modification (one line) in plasmray.f to enable ray to
c[46] start in a case with small plasma region outside the the
c[46] LCFS (HSX test case).  [SAP, 051019].

c[45] Addition of UHR to genray pgplot output [EricNM,BH 051018].

c[44] Modification (one line)in forest.f to enable ray to start in
c[44] the correct mode (O-mode) near plasma edge in an
c[44] OXB test case [SAP, 051018].

c[43] Facility to switch space/time step along the ray to a
c[43] shorter value in the vicinity of the upper hybrid layer
c[43] was added. This enables resolution of the UHL when the
c[43] it is very close to the plasma boundary, and has been
c[43] found useful to OXB. New namelist variables are
c[43] i_uh_switch,uh_switch,prmt6_uh_switch [SAP, 051011].
c
c[42] Removed old calculation of RF driven currents,
c[42] scurntor,scurntpol,curdentr,curdenpl which were output
c[42] to the netcdf file.  These have been replaced by
c[42] s_cur_den_xxx, where xxx=parallel,toroidal,poloidal,onetwo,
c[42] which have more accurate accounting of geometry [SAP, BH 050928]. 

c[41] Add the absorption calculations iabsorp=11 using
c[41] the A.Ram relativistic tensor, dispersion function 
c[41] D = determinant and the projection
c[41] method for Im(N_perp)=-=-Im(D)/(dD/dRe(N_perp))  [SAP, 050928]

c[40] Changed sign of output currents to conform with positive
c[40] current in positive toroidal direction, except for the
c[40] parallel current density which is referenced to the dirn
c[40] of the magnetic field.   Thus toroidal and "onetwo" current
c[40] are referenced to toroidal angle dirn.  [SAP, BH 050923].

c[39] Added emission related data to netcdf output file for
c[39] case of i_emission=1  [SAP, 050923].

c[38] Add namelist relres to set accuracy level in Abhay Ram
c[38] dispersion relation. Added additional namelist parameters
c[38] if user wants finer control of Ram dispersion (errabs0,
c[38] errrel0,navg,diff_err), and uncombined Nelson-Melby and
c[38] Ram dispersion tensors [ENM,050903].

c[37] Corrected the calculation of the length along the ray,
c[37] which is used in emission calculations. Affects non-perp
c[37] viewing [SAP, 050817].

c[36] Added xdraw plotting (drawemres.in) of the resonanse ellipse data
c[36] p_perp_max,p_parallel_min,p_parallel_max. 
c[36] The code uses data for all harmonics, but writes the data
c[36] in em_res.bin file (from mk_graph.f) only for the second harmonic.
c[36] Also, added a new parameter, n_relt_harma, which is the maximal  
c[36] number of harmonics used for relativistic anti-hermitian tensor
c[36] calculations[SAP, 050817].

c[35] Add the input parameter i_output (in nmlist) to choose the output
c[35] step, using the poloidal or the total length along the ray [SAP].

c[34] BXO emission calc for single frequency,  expanded to
c[34] multiple frequencies using i_emission=1,i_ox=2 [SAP,BH 050721].  

c[33] New option controlled by i_ox and related variables in 
c[33] namelist section /ox/ gives optimum launch angles for OXB
c[33] from a give antenna location (i_ox=1) in file ECcone_optimal.dat
c[33] and then with i_ox=2 can launch rays (using optimized angles or
c[33] others in /eccone/ ) obtaining OXB ray trajectories and
c[33] transmission coefficients [SAP, BH 050623]. 

c[32] Add capability for variable number of radial bins NR (new nmlst)
c[32] up to NRA (parameter) for tabulation of power deposition and
c[32] current drive  [SAP, 050614].

c[31] Re-worked istart=3 case for OXB launch at the optimal n_par,
c[31] launching the X-mode at Xe .gt.1.  Need ioxm=-1, instead of +1.
c[31] [SAP, 050510]. 

c[30] The relativistic absorption using Stix formula and relativistic
c[30] dielectric (electron plasma) tensor and electric field polarization 
c[30] iabsorp=10  [SAP, 050408]

c[29] The hot plasma absorption using Stix formula and hot plasma
c[29] polarization iabsorp=9  [SAP, 050408]

c[28] If (partner.eq.transport) the code reads the file created by
c[28] transport code: genray_profs_in and writes the power and current
c[28] density profiles in the output file: genray_profs_out 
c[28] Works for electron only (nbulk=1) of multi-species [SAP,BH 050408]

c[27] The coordinates of the ray starting points will be written in
c[27] output *.nc file for grill launching [SAP, 050408]  

c[26] New option added for calculation of grill ray starting conditions:
c[26]
c[26] For i_grill_pol_mesh=2 GENRAY creates non-uniform poloidal mesh
c[26] for the power distribution according to the used poloidal power
c[26] distribution. All rays at this mesh have the same power.
c[26]
c[26] For i_grill_npar_ntor_npol_mesh=2 GENRAY creates non-uniform
c[26] meshes for the specified power spectrum:
c[26] a) N_parallel mesh (for i_n_poloidal=1,2,3) for power(N_parallel)
c[26] b) N_toroidal and N_poloidal meshes (for i_n_poloidal=4) for
c[26] power(N_toroidal)*power(N_poloidal) .
c[26]
c[26] All rays will have the same power at these meshes.  [SAP,050316].

c[25] Added further option for small radius coordinate:
c[25] rho=(R_max(psi)-R_min(psi))/(R_max(psi_lim)-R_min(psi_lim))
c[25] accessed through namelist value  indexrho=6  [SAP, 050310].

c[24] Added option for the fully relativistic combined dispersion 
c[24] function (id=14) constructed by Nelson-Melby, which uses
c[24] 1) Eric Nelson-Melby tensor based on Weiss method (npar.le.0.3):
c[24] Isaac Weiss, JCP 61, p. 403-416 (1985).
c[24] 2) Abhay Ram tensor (npar.gt.0.3). 
c[24] To use these dispersion functions some Slatec subroutines 
c[24] were added into genray.  [SAP, 050302].

c[23.0] Boundary condition for the spline of the radial density
c[23.0] and temperature profiles was adjusted to remove negative
c[23.0] values for sharp edge profiles near the edge.
c[23.0] [SAP, 050226]

c[23] A trajectory shift in small radius direction was added, accounting
c[23] for OX mode-conversion. The Mjolhus transmission coefficient
c[23] is used, and is added to the .nc file  [SAP, 050218]. 
  
c[22] The calculations of the OX optimal direction of EC cone central ray 
c[22] Namelist OX data were added into the input file genray.dat 
c[22] for these calculations.
c[22] Parameters of the optimal direction of the central O mode are
c[22] output to the ECcone_optimal.dat file.    [SAP, 050218].

c[21] For iabsorp=3 (Chiu et al FW absorp), added cnprim_s(i),i=2,nbulk 
c[21] calculation of individual ion component absorptions.
c[21] Added breakdown of individual ion absorption along rays to .nc. 
c[21] Output radial profiles of individual ion absorption and
c[21] powers, radial profiles to density, temp, and zeff, bin 
c[21] volume/area/pol-length, and some wave specification parameters 
c[21] to the .nc file. [BH041022].

c[20] Put rewind(unit=1) before each namelist section read, so
c[20] order of sections in the namelist input can be arbitrary
c[20] [BH,041006].

c[19] CVS tagged version genray_v3-4, and moving CVS repository to
c[19] to PPPL  [BH, 040917].

c[18] For Westerhof-Tokman ray equations, poloidal damping calculated 
c[18] Tokman fluxes.  Also, added magnetic field components to .nc 
c[18] output files.  Updated CVS, designated version 
c[18] genray_cvs_v3-3_040915    [SAP, 040908; BH, 040915]. 

c[17] Added id=13, Westerhof-Tokman disp with Stix hot-plasma tensor
c[17]              iherm =1 hermit dielectric tensor, 2-full
c[17] [SAP, 040907;  Updated 040907, BH] 

c[16] Added id=11, dispersion with Eric Nelson-Melby tensor
c[16]              iherm =1 hermit dielectric tensor, 2-full
c[16]         =12 Westerhof-Tokman dispersion with Nelson-Melby tensor
c[16]              iherm =1 hermit dielectric tensor, 2-full
c[16] [SAP, 040830;  Updated 040904, BH] 

c[15] Added iabsorp=7 option, providing absorption using cold
c[15] plasma polarizations together with the numerical anti-Hermitian
c[15] calculation.  This is relevant for 3rd and higher EC harmonic
c[15] absorption.  Updated from previous SAP GR020408 code, which
c[15] got left behind. Updated CVS, designated version 
c[15] genray_cvs_v3-1_040827 [BH,040827].

c[14] Modified name of variable to feqd=R*Bphi from f, to avoid
c[14] confusion it with other f variables. Changed output
c[14] of 3d.dat and 3d.nc files from mod(sbtot) to negative sbtot when
c[14] feqd is negative at the magnetic axis, in accord with toray.
c[14] Add new namelist variable, eqdskin, default="equilib.dat".
c[14] Added new namelist char variable mnemonic, changing name of 3d.dat
c[14] and 3d.nc files to mnemonic.txt and mnemonic.nc.  Output is
c[14] controlled by additional namelist variable rayop.
c[14] Added namelist variable dielectric_op controlling eps output.
c[14] Namelist input file can now be genray.in or genray.dat.
c[14] Updated CVS, designated version genray_cvs_v3-0_040820.
c[14] [BH, 040822].

c[13] Added Westerhof-Tokman ray tracing plus additional term in
c[13] energy flux vector, accounting for effect of finite 
c[13] anti-Hermitian dielectric tensor, id=10.  This uses the
c[13] Mazzucato-Fidone-Granata disp tensor. Added optional 
c[13] netcdf output of dielectric tensor elements [SAP,040804].

c[12] Added power density and calculated current density (from
c[12] TorGA_curba) to netcdf output file.
c[12] Obtained agreement between GENRAY and TORAY on ECCD, using 
c[12] Lin-Liu memo \"Driven Current Density and Toroidal Current
c[12] in TORAY-GA [SAP,040729]\".

c[11] Added TORAY system for Gaussian EC ray pattern, to facilitate
c[11] direct comparisons.  [BH, 040412].

c[10] Added capability to read nbulk.le.nbulka, ndens.le.ndensa,
c[10] where nbulk,ndens are namelist variables and nbulka,ndensa
c[10] maximumal parameters.
c[10] It is necessary to have nbulk and ndens defined in genray.in.
c[10]  ==> genray_cvs_v2-3_040122. [SAP,040122]

c[9] Smirnov has fixed interpolation of radial set of distributions,
c[9] for emission, and instituted max dimensions in v,theta,rho 
c[9] directions for input (nonthermal) distributions
c[9]  ==> genray_cvs_v2-3_040116. [SAP,040116]

c[8] Smirnov brings GENRAY up to date with changes made in 
c[8] San Diego, Jan\'03. See ./code/readme030123.
c[8] Changes included: 
c[8] The new b.f uses new splines: coeff1, coeff2, terp1 and terp2p
c[8] to calculate f(psi) [now, feqd]  and psi(r,z).  [SAP, 030124]

c[7] Reconciling version of genray brought from Moscow to SD by
c[7] Smirnov with genray_cvs version as of 021231.
c[7] Main changes related to: 
c[7] -Corrections, derivs in hot plasma corrected for eV
c[7] -Ono dispersion for FW added
c[7] -Plotting using PGPLOT of complex dispersion relation
c[7] -Plotting of magnetic field and plasma quantities using PGPLOT
c[7] [bobh,smirnov, 030122].

c[6] Substituted smirnov_021021 bug fixes in absorpfw.f (Z==>Z_0).
c[6] Still need to integrate new prep3d.f, to get Chiu polarization
c[6] output with iabsorp=3 damping of FW [bobh, 021103].

c[5] Adding "character*64 version" to param.i, printing it to std out,
c[5] and putting in netcdf output file.  Also, printing out time
c[5] of execution, and machine.  Cvs tagged and exported version
c[5] genray_v1-1.       [bobh, 020828]. 

c[4] Fixed coding inadaquacy in poloidal distn of power, grill_lh.f.
c[4] Added igrilltw=1 option for uniform pol distn of power, grill_lh.f.
c[4] Added igrillpw=3 option for guassian distn of power vs npar,
c[4] grill_lh.f.   [bobh, 020826].

c[3] Added file, a_change.i, to keep track of significant changes.
c[3] [bobh, 020824].

c[2] Added capability of eqdsk arrays for feqd=R*Bphi,ff',p,pp',q to
c[2] have smaller dimensionality than the radial array, nxeqd, for psi.
c[2] Fixed a bug in coeff1 spline routine, and in equilib.f which caused
c[2] code not to work for nxeqd.ne.nxmax. This feature still needs work.
c[2] [bobh, 020824].

c[1] GENRAY put under control of CVS at compxco.com.
c[1] [bobh, 020814]

c[0] GENRAY is described in Genray_manual_CompX-2001-1-V1.pdf, by 
c[0] Smirnov and Harvey. Additional functionality for electron cyclotron
c[0] emission has been added in Jan-Feb\'02 by Smirnov and Harvey, not
c[0] described in the manual.  [This is status as of July\'02].
