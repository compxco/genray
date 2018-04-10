c-----input data for LSC approach
c-----from namelist/lsc_approach_nml/

c     velocity mesh is normalized at v_e_thermal which is different
c     at different radial points (1:n_psi_TSC)
   
      real*8  v_par_lsc_dev_vte_max !maximal normalized velosity 
                                    !of the velocity mesh
      integer
     &i_lsc_approach,   !=0 do not use lrs_approach
                        !=1 use lrs_approach
     &nv_lsc_ql_lh,     !number of the velocity mesh points 
     &n_psi_TSC,        !number of radial points for TSC arrays 
     &icount_iter_max,  ! max number of iteratins to get QL diffusion
                        ! coefficients. 
                        ! It is used in subroutine LSC_power_iterations
     &n_power_steps,    !number of iterations in initial power
                        !initial power at each ray will be
                        !(i_power_step*1.d0/n_power_steps)*powinilh(iray) 
     &n_dim_radii_f_plots,   ! number of radii points from TSC mesh
                             ! rho_bin_center(1:n_psi_TSC)
                             ! at which obtained distribution function 
                             ! ln(f_e(v_parallel)) will be plotted 
                             ! into the plot.ps file
                             ! It should be  n_dim_radii_f_plots.le.n_psi_TSC 
                             ! 
                             ! If n_dim_radii_f_plots=0 then
                             ! the codfe will not plot f_e
     &n_radii_f_plots_ar     !(1:n_dim_radii_f_plots) numbers of radii points
                             ! from TSC mesh
                             ! rho_bin_center(1:n_psi_TSC)
                             ! at which obtained distribution
                             ! function ln(f_e(v_parallel))
                             ! will be plotted 
                             ! into the plot.ps file


      real*8
     &EdcTSC_1D,        ! E_DC along B field   [V/m]                  
     &JparTSC_1D,       ! current density parallel B fielf [A/m**2]
     &EdcTSC_min,       ! minimal value of  abs(E_DC) [V/m]   
     &eps_diff_d_ql     ! accuracy of iteratins to get QL diffusion
                        ! coefficients:
                        ! diff_d_q < eps_diff_d_ql.   
                        ! diff_d_q=||d_ql_lsc_lh_ar_old- d_ql_lsc_lh_ar||
                        ! It is used in subroutine LSC_power_iterations


      common/lsc_approach_nml/
     &EdcTSC_1D(1:ndensa),          !(1:n_psi_TSC)
     &JparTSC_1D(1:ndensa),         !(1:n_psi_TSC)
     &EdcTSC_min,
     &v_par_lsc_dev_vte_max,    
     &eps_diff_d_ql,
     &nv_lsc_ql_lh,
     &i_lsc_approach,
     &n_psi_TSC,
     &icount_iter_max,
     &n_power_steps,   
     &n_dim_radii_f_plots,  
     &n_radii_f_plots_ar(1:ndensa)  !(1:n_psi_TSC) 
