
c      the data for the profiles of the absorbed power and current drive.
c      PARAMETER (NRA) is set in param.i

      include 'onetwo_nml.i'
      include 'onetwo_no_nml.i'
c
c    NR is the maximal number of points in the small radius direction 
c      for the
c      calculations of the power and current drive radial profiles 
c      in each radial bin.
c    NRA= maximal value of NR. It is set in param.i
c    power_e is the absorbed power due to collisionless electron damping
c    power_i is the absorbed power due to collisionless ion damping
c    power_s is the absorbed power due to coll'less ion damping on each ion
c    power_cl is the absorbed power due to collisional damping
c    spower_e summed absorbed power due to coll'less electron damping
c    spower_i summed absorbed power due to collisionless ion damping
c    spower_i summed absorbed power due to coll'less ion damping on each ion
c    spower_cl summed absorbed power due to collisional damping
c    allcur     total current
c    densprof,temprof,zefprof   for output of radial profiles to .nc file
c
c050406
c    currtot is a total current along B field
c
c    cur_den_parallel is the flux surface averaged parallel current density
c                     from one ray
c                    
c   s_cur_den_parallel is the flux surface averaged parallel current density
c                     from all rays
c
c   binarea_pol is bin area of the poloidal cross-section at theta poloidal=0 
