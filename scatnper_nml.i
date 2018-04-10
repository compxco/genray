
c     the data for the scattering of the perpendicular refractive index
c     parameter nscat_n ! It should be given in the param.i file
c     All data are from namelist 
     
      
      real*8
c-----from namelist /scatnper/
     &rhoscat,
     &scatd
c-----from namelist /scatnper/
      integer
     &iscat,
cSAP120511
     &iscat_lh_nicola

      common/scatnper_nml/
     &rhoscat(1:nscat_n),
     &scatd(0:nscat_n),
     &iscat,
cSAP120511
     &iscat_lh_nicola
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
c      iscat_lh_nicola=1 to use scattering procedure lh_scattering
c                     =0 do not use lh_scattering
c                     








