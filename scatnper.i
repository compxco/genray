
c     data for the scattering of the perpendicular refractive index
c     parameter nscat_n ! It should be given in the param.i file
 
      include 'scatnper_nml.i'
      include 'scatnper_no_nml.i'
     
c-------------------------------------------------
c      iscat it is the switch for the n_perp scattering
c            iscat=1 the scattering swithed on, =0 swithed off
c
c      rhoscat(1:nscat_n) small radii for the scattering location
c
c      The scattering of the polar angle deltheta will be
c      deltheta=dsqrt(2.d0*scatd)*ranorm(fseed) 
c      scatd(0)         the mean square scattering angle (radians**2)
c                       for the plasma boundary reflection points   
c               
c      scatd(0:nscat_n) the mean square scattering angles (radians**2)
c                       for the interior plasma boundary points
c      rhooldsc the value of the small radius from the previous time step
c               It will be calculated in genray and output subroutines.   
c-------------------------------------------------------------------










