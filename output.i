
c     The data for  subroutine outpt
    
      include 'output_nml.i'
      include 'output_no_nml.i'


c----------------------------------------------------------------------------
c    id_initial !to change Ono dispersion
c
c    first      !to change Ono dispersion
c
c    was_not_ox_conversion !for OX jump
c    n_plot_disp is the number of major radius points to plot D 
c    r_plot_disp(n_plot_dispa) are major radius where the dispersion
c    function contouts plot D(ImN_perp,ReN_perp) will be ploted
c    id_plot_disp is the dispresion function type for contours plots
c
c    s_poloid_plot_disp() are poloidal distences [m] where the dispersion
c    function contouts plot D(ImN_perp,ReN_perp) will be created
c
c    id_plot_disp is the dispresion function type for contours plots
c
c    point_plot_disp ='poloidl_dist' to create D contours at given
c                      s_poloid_plot_disp()
c                    ='major_radius'  to create D contours at given
c                      r_plot_disp()
c--------------------------------------------------------------------------
c    n_plot_disp_cold    is the number of major radius or poloidal
c                        distance points to plot D_cold(Re_nperp)
c 
c    s_poloid_plot_disp_cold() are poloidal distences [m] where the dispersion
c                              function function plot D(ReN_perp) will be 
c                              created
c
c    r_disp_cold() are major radius [m] where the dispersion
c                   function function plot D(ReN_perp) will be created
c
c    point_plot_disp_cold ='poloidl_dist' to create D(Re_N_perp) plots at given
c                      s_poloid    integer n_contour_plot_disp !is the number of contour_plot_disp_cold()
c                    ='major_radius'  to create D(ReN_per) plots at given
c                      r_plot_disp_cold()
c---------------------------------------------------------------------------
c     for plotting dispersion function contours D(ImN_perp,ReN_perp)
c---------------------------------------------------------------------------
c number_map_points_real_nperp  is the number of map points
c                               in Real(N_perp) direction
c
c number_map_points_image_nperp  is the number of map points
c                                in Image(N_perp) direction
c
c ratio_min_r_nperp,ratio_max_r_nperp  set the ratio of
c          minimal and maximal map boundaries in Real N_perp direction
c          (min_r_nperp < Real(N_perp) < max_r_nperp) to the value of
c          Real(N_perp_ray_=cnper along the ray:
c          min_r_nperp= Real(N_perp_ray)*ratio_min_r_nperp
c          max_r_nperp= Real(N_perp_ray)*ratio_max_r_nperp
c          These parameters should be: 
c           0 =< ratio_min_r_nperp < 1
c           1 <  ratio_max_r_nperp 
c
c ratio_min_i_nperp,ratio_max_i_nperp  set the ratio of
c          minimal and maximal map boundaries in Image N_perp direction
c          (min_i_nperp < Image(N_perp) < max_i_nperp) to the value of
c          Image(N_perp_ray)=cnprim along the ray:
c          min_i_nperp= Image(N_perp_ray)*ratio_min_i_nperp
c          max_i_nperp= Image(N_perp_ray)*ratio_max_i_nperp
c          These parameters should be: 
c           0 =< ratio_min_i_nperp < 1
c           1 <  ratio_max_i_nperp 
c
c          If Image(N_perp_ray) < 1 then the code will set
c          following map boundaries: min_i_nperp=0 and man_i_nperp=1. 
c
c n_contour_plot_disp !is the number of contours for D(ReN_perp,ImN_perp)
c----------------------------------------------------------------------------
c    frequencies plot along the straight line
c                   r_freq,z_freq, the line edge point [m]
c                   alpha_freq,beta_freq are angles[degree]
c                   dist_freq is the line length   [m]
c                   nsteps_freq, the number of points for plot
c                   n_ec_harmonics_freq  is the number of plotted 
c                                        ec harmonics
c                   max_plot_freq 
c                   npar_freq   N_parallel to plot X mode cutoof
c------------------------------------------------------
