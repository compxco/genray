      integer
     & i_adj, ! =1 to use adj calculations
              ! =0 do not use adj_calculastions 
     & npsi0, ! number of radial points  
     & nthp0, ! number of poloidal points for the integration
              ! along b filed line
     & nmax_chi,  ! the number of grid points in u_0 direction (for adj mesh)
     & imax_chi,  ! the number of grid points in pitch angle 
     & lmax_chi,  ! is the maximal number of Legandre harmonic
              ! lmax=1 corresponds to collision with a Maxwellian.
              !      Usually lmax=21
     & tmax,  ! is the maximal number of time steps
cSm071009
     & n_harm_adj_min,  ! minimal and maximal harmonics for power and CD
     & n_harm_adj_max,  ! calculations using ADJ function
     & n_relt_intgr_adj, ! number of points for Power and CD integration
                         ! in p_perp0 direction along resonce curve
     & i_resonance_curve_integration_method_adj,!chooses integration method
                                                !for power and CD calculation
     & i_calculate_or_read_adj_function,        ! =1 calculate chi function
                                                ! =0 read chi function from
                                                !    the file: adjout
     & i_chi_interpolation            !to chose interpolation method
                                      !for derivatives from chi function
                                      ! =1 using spline
                                      ! =2 linear interpolation
                                      !    of numerical derivatives

 
      real*8
     &t,      ! the temperature in units of m*u_n**2
              ! If you want velocities normalized to sqrt(T/m), set t=1.
     &ze,     ! ze is the electron charge state. It is usually =1, 
              ! but for Lorentz limit, ion charge Zi ==> infinity,
              ! set  ze=0 and zi=1
     &dt,     ! time step in adj solver, usually =1
     &umax,   ! the maximal value of u_0. 
              ! The grid is spaced between 0 and umax.
     &alpha,  ! is a weighting factor in the time stepping algorithm.
              !       =0 for explicit algorithm
     &rho_,   !governs the treatment  of the boundary at u=umax
              !      for rho_=1 the second derivative from the adjoint function chi
              !      d^2(chi)/dx^**2(at u=umax)=0
              !      for rho_=0 the first derivative =0
     &aerrmx, ! some error for subadj 
     &rerrmx,  ! some error for subadj 
     &epsi_adj !the accuracy used in integration by adaptive Simpson
               !Works at i_resonance_curve_integration_method_adj=4 case only

      common /idj_nml/
     & t,ze,umax,dt,  
     & aerrmx,rerrmx,epsi_adj,
     & npsi0,nthp0,nmax_chi,imax_chi,lmax_chi,tmax,alpha,rho_,
     & i_adj,
cSm071009
     & n_harm_adj_min,n_harm_adj_max,n_relt_intgr_adj, 
     & i_resonance_curve_integration_method_adj,
     & i_calculate_or_read_adj_function,
     & i_chi_interpolation          
!-----------------------------------------------------------------------
! The following parameters are for adj function by Charles Karney.
! The velocity is normalized to arbitrary u_n (which is usually either
! sqrt(T/m) or c), times to tau_n(G^e/e/u_n**3)**-1
!-----------------------------------------------------------------------
!-------------------------------------------------------
! t is the temperature in units of m*u_n**2
! If you want velocities normalized to sqrt(T/m), set t=1.
!-----------------------------------------------------------------------
! c2 is the square of speed of light (normalized to u_n**2).
! If t=1 is chosen, then c2=512/Te_keV.
!-----------------------------------------------------------------------
! ze is the electron charge state. It is usually =1, 
! but for Lorentz limit, ion charge Zi ==> infinity,
! set  ze=0 and zi=1
!-----------------------------------------------------------------------
! zi is ion charge state
!-----------------------------------------------------------------------
!! eps is the inverse aspect ratio r/R. This is for specific 
!    cases
!-----------------------------------------------------------------------
! nmax_chi is the number of grid points in u_0 direction (for adj mesh)
!-----------------------------------------------------------------------
! umax is the maximal value of u_0. 
! The grid is spaced between 0 and umax.
!-----------------------------------------------------------------------
! imax_chi is the number of grid points in pitch angle 
!-----------------------------------------------------------------------
! lmax_chi is the maximal number of Legandre harmonic
!      imax=1 corresponds to collision with a Maxwellian.
!      Usually lmax=21
!-----------------------------------------------------------------------
! kmax  is the number of integration points for field line integration.
!       ~1000  
!-----------------------------------------------------------------------
! dt is the time step
!-----------------------------------------------------------------------
! tmax is the maximal number of time steps
!-----------------------------------------------------------------------
! alpha is a weighting factor in the time stepping algorithm.
!       =0 for explicit algorithm
!-----------------------------------------------------------------------
! rho_ governs the treatment  of the boundary at u=umax
!      for rho_=1 the second derivative from the adjoint function chi
!      d^2(chi)/dx^**2(at u=umax)=0
!      for rho_=0 the first derivative =0 
!-----------------------------------------------------------------------
!  i_resonance_curve_integration_method_adj=1 !rectangle integration
!                                              !over angle,
!                                              !for ellipse case only
!       i_resonance_curve_integration_method_adj=2 !rectangle formula,
!                                              !p_perp0 integration
!                                              !elliptic or hyperbolic case
!       i_resonance_curve_integration_method_adj=3 !trapezoidal formula,
!                                              !p_perp0 integration
!                                              !elliptic or hyperbolic case
!       i_resonance_curve_integration_method_adj=4 !adaptive Simpson integration 
!                                              !p_perp0 integration
!                                              !elliptic or hyperbolic case
!
!       i_resonance_curve_integration_method_adj is used 
!       to choose the numerical integration method for 
!       power and CD calculation.
!-----------------------------------------------------------------------
!  i_calculate_or_read_adj_function          ! =1 calculate chi function
!                                            ! =0 read chi function from
!                                            !    the file: adjout
!--------------------------------------------------------------------------
! i_chi_interpolation            !to chose interpolation method
!                                !for derivatives from chi function
!                                ! =1 using spline
!                                ! =2 linear interpolation
!                                !    of numerical derivatives
!---------------------------------------------------------------------------
