
c      the data for the profiles of the absorbed power and current drive.
c      PARAMETER (NRA) is set in param.i
c      No namelist data
   
c-----real*8
      ! Set in p_c_prof (accumulate values from each ray element): 
      real*8 allpw_e, allpw_i, allpw_cl, allpower, allcur
      
      ! Set in dnonetwo:
      real*8 powtot_e,powtot_i,powtot_cl,
     &       currtot,
     &       parallel_cur_total,toroidal_cur_total,poloidal_cur_total
     
      ! Set in p_c_prof (add contribution from each ray element): 
      real*8,  pointer ::
     &         power(:),            !(NR)
     &         current(:),          !(NR)
     &         power_e(:),          !(NR) 
     &         power_i(:),          !(NR) 
     &         power_s(:,:),        !(NR,nbulk)
     &         power_cl(:),         !(NR)
     &         cur_den_parallel(:)  !(NR)   
      
      ! Set in other subroutines:
      real*8,  pointer ::
     &         temparr(:),            !(NR)    ! Not used
     &         zeffarr(:),            !(NR)    ! Not used/spldens?
     &	       spower(:),             !(NR)    ! sonetwo
     &         scurrent(:),           !(NR)    ! sonetwo
     &         powden(:),             !(NR)    ! dnonetwo
     &         currden(:),            !(NR)    ! dnonetwo 
     &         spower_e(:),           !(NR)    ! sonetwo
     &         spower_i(:),           !(NR)    ! sonetwo
     &         spower_s(:,:),         !(NR,nbulk)sonetwo
     5	       spower_cl(:),          !(NR)    ! sonetwo
     &         powtot_s(:),           !(nbulk) ! dnonetwo
     &         powden_e(:),           !(NR)    ! dnonetwo
     &         powden_i(:),           !(NR)    ! dnonetwo
     &         powden_s(:,:),         !(NR,nbulk)dnonetwo
     &         powden_cl(:),          !(NR)    ! dnonetwo
     &         currden_s(:),          !(NR)    ! dnonetwo
     & 	       powden_e_s(:),         !(NR)    ! Not used
     &         powden_i_s(:),         !(NR)    ! Not used
     &         powden_cl_s(:),        !(NR)    ! Not used
     & 	       densprof(:,:),         !(NR,nbulk) ! netcdfr3d
     &         temprof(:,:),          !(NR,nbulk) ! netcdfr3d
     &         zefprof(:),            !(NR)       ! netcdfr3d
     &         rhoprof(:), ![2024-01-24] rho corresponding to densprof (just for nc file)
     &         rho_bin(:),            !(NR)       ! netcdfr3d
     &         rho_bin_center(:),     !(NR-1)     ! netcdfr3d    
     &         binvol(:),             !(NR-1)  ! dnonetwo and p_c_prof
     &         binarea(:),            !(NR-1)  ! dnonetwo and p_c_prof
     &         binarea_pol(:),        !(NR-1)  ! dnonetwo and p_c_prof
     &         pollen(:),             !(NR-1)  ! dnonetwo 
     &         s_cur_den_parallel(:), !(NR-1)  ! dnonetwo and sonetwo
     &         s_cur_den_onetwo(:),   !(NR-1)  ! dnonetwo  
     &         s_cur_den_toroidal(:), !(NR)    ! dnonetwo
     &         s_cur_den_poloidal(:)  !(NR)    ! dnonetwo
  


      common/onetwo_no_nml/ powtot_e,powtot_i,powtot_cl,allpower,
     &         allpw_e,allpw_i,allpw_cl,
     &         allcur,currtot,  
     &         parallel_cur_total,toroidal_cur_total,poloidal_cur_total,
c-----pointers
     &         powtot_s,
     &         power,
     &         current,
     &         temparr,
     &         zeffarr,
     &	       spower,
     &         scurrent,
     &         powden,
     &         currden,
     &         power_e,
     &         power_i,
     &         power_s,
     &         power_cl,
     &         spower_e,
     &         spower_i,
     &         spower_s,
     &	       spower_cl,
     &         powden_e,
     &         powden_i,
     &         powden_s,
     &         powden_cl,
     &         currden_s,
     & 	       powden_e_s,
     &         powden_i_s,
     &         powden_cl_s,
     8	       densprof,
     &         temprof,
     &         zefprof,
     &         rhoprof, ![2024-01-24] rho corresponding to densprof (just for nc file)
     &         rho_bin,
     &         rho_bin_center,
     &         binvol,
     &         binarea,
     &         binarea_pol,
     &         pollen,
     &         cur_den_parallel,
     &         s_cur_den_parallel,
     &         s_cur_den_onetwo,
     &         s_cur_den_toroidal,
     &         s_cur_den_poloidal
    
c-----------------------------------------------------------------------------
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
c    spower_s summed absorbed power due to coll'less ion damping on each ion
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
