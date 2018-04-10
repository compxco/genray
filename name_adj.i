c------------------------------------------------------
c     namelist /adj_nml/  
c-------------------------------------------------------

      namelist /adj_nml/ 
     & i_adj,
     & t,ze,umax,dt,
     & aerrmx,rerrmx,epsi_adj,
     & npsi0,nthp0,nmax_chi,imax_chi,lmax_chi,tmax,alpha,rho_,
     & n_harm_adj_min,n_harm_adj_max,
     & n_relt_intgr_adj,i_resonance_curve_integration_method_adj,
     & i_calculate_or_read_adj_function,    
     & i_chi_interpolation  
