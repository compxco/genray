c-----------------------------------------------------------
c     The spline coefficients for the emission calculations
c     z,r,phi,N_z,N_r,M are the coordinates of the points along the ray
c                       in 6D phase space (3R,3N)
c      Here are no data from namelists.
c-----------------------------------------------------------     
c-----------------------------------------------------------
c     The spline coefficients for the emission calculations
c     z,r,phi,N_z,N_r,M are the coordinates of the points along the ray
c                       in 6D phase space (3R,3N)
c      Here no data from namelists.
c-----------------------------------------------------------     


      real*8, pointer ::
     &cx_z(:),     !(nrelta4),
     &tr_z(:),     !(nrelta4),
     &cx_r(:),     !(nrelta4),
     &tr_r(:),     !(nrelta4),
     &cx_phi(:),     !(nrelta4),
     &tr_phi(:),     !(nrelta4),
     &cx_cnz(:),     !(nrelta4),
     &tr_cnz(:),     !(nrelta4),
     &cx_cnr(:),     !(nrelta4),
     &tr_cnr(:),     !(nrelta4),
     &cx_cm(:),      !(nrelta4),
     &tr_cm(:),      !(nrelta4),
c-----for ebw-x mode ray part 
     &cx_z_x(:),     !(nrelta4),
     &tr_z_x(:),     !(nrelta4),
     &cx_r_x(:),     !(nrelta4),
     &tr_r_x(:),     !(nrelta4),
     &cx_phi_x(:),     !(nrelta4),
     &tr_phi_x(:),     !(nrelta4),
     &cx_cnz_x(:),     !(nrelta4),
     &tr_cnz_x(:),     !(nrelta4),
     &cx_cnr_x(:),     !(nrelta4),
     &tr_cnr_x(:),     !(nrelta4),
     &cx_cm_x(:),      !(nrelta4),
     &tr_cm_x(:),      !(nrelta4),    
c-----for o mode ray part 
     &cx_z_o(:),       !(nrelta4),
     &tr_z_o(:),       !(nrelta4),
     &cx_r_o(:),       !(nrelta4),
     &tr_r_o(:),       !(nrelta4),
     &cx_phi_o(:),     !(nrelta4),
     &tr_phi_o(:),     !(nrelta4),
     &cx_cnz_o(:),     !(nrelta4),
     &tr_cnz_o(:),     !(nrelta4),
     &cx_cnr_o(:),     !(nrelta4),
     &tr_cnr_o(:),     !(nrelta4),
     &cx_cm_o(:),      !(nrelta4),
     &tr_cm_o(:),      !(nrelta4),    
c-------------------------------      
     &kinetic_energy_kev (:)   !(jx_kin_a)
  

      common/emissa/
     &cx_z,
     &tr_z,
     &cx_r,
     &tr_r,
     &cx_phi,
     &tr_phi,
     &cx_cnz,
     &tr_cnz,
     &cx_cnr,
     &tr_cnr,
     &cx_cm,
     &tr_cm,
c-----for ebw-x mode ray par     
     &cx_z_x,
     &tr_z_x,
     &cx_r_x,
     &tr_r_x,
     &cx_phi_x,
     &tr_phi_x,
     &cx_cnz_x,
     &tr_cnz_x,
     &cx_cnr_x,
     &tr_cnr_x,
     &cx_cm_x,
     &tr_cm_x, 
c-----for o mode ray part      
     &cx_z_o,
     &tr_z_o,
     &cx_r_o,
     &tr_r_o,
     &cx_phi_o,
     &tr_phi_o,
     &cx_cnz_o,
     &tr_cnz_o,
     &cx_cnr_o,
     &tr_cnr_o,
     &cx_cm_o,
     &tr_cm_o,     
c-----for emission spectrum
     &kinetic_energy_kev
     
