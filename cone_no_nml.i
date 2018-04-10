   
c   The data for EC cone starting conditions, obtained from namelist
c   and subroutine dinit_mr calls.

c     nraymax must be greater or equal than nray (obtained in code).

    
      real*8      
     1	          powj,zstj,rstj,
     1            phistj,betaj,alphaj,
     1            z_st_ox,r_st_ox,phi_st_ox,alpha_st_ox,beta_st_ox, 
     1            disk_power_equal,
     1            power_per_ray,
     1            disk_rho_mesh_edge,disk_rho_mesh_center,
     1            ar_z_disk_launch,  
     1            ar_r_disk_launch,  
     1            ar_phi_disk_launch,        
     1            ar_beta_disk_launch,      
     1            ar_alpha_disk_launch       

      integer
     1            iray,
     1            icone_loc

      common/cone/
     1	          powj(nraymax),zstj(nraymax),rstj(nraymax),
     1            phistj(nraymax),betaj(nraymax),alphaj(nraymax),
     1            z_st_ox,r_st_ox,phi_st_ox,alpha_st_ox,beta_st_ox,
     1            disk_power_equal(n_mesh_disk_radial_bin_a),
     1            power_per_ray,
     1            disk_rho_mesh_edge(n_mesh_disk_radial_bin_a+1),
     1            disk_rho_mesh_center(n_mesh_disk_radial_bin_a),
     1            ar_z_disk_launch(nraymax),  
     1            ar_r_disk_launch(nraymax),  
     1            ar_phi_disk_launch(nraymax),        
     1            ar_beta_disk_launch(nraymax),      
     1            ar_alpha_disk_launch(nraymax),       
     1            iray,
     1            icone_loc
c-----------------------------------------------------------------
c      ncone- number of ray starting cones
c      zst,rst,phist- are coordinates of the antenna cone vertices
c      powtot is power in each cone
c      alfast,betast,alpha1,alpha2 give ray starting angles and
c        cone angles
c      The variables dimensioned nraymax give starting locations
c        and starting angles for each of the ray.
c      alphaj(iray)- toroidal angles(radian) of ray iray
c      betaj(iray)- angles(radian) between the  horizontal
c                          plane and ray iray
c      powj(iray)-initial power flowing in the ray channel at antenna
c                    normalized Sum_{i=1,nray}powj(i)=
c                    total power in all rays(erg/sec)
c      iray         -     ordinal number of the ray
c   
c       z_st_ox,r_st_ox,phi_st_ox,alpha_st_ox,beta_st_ox -
c       the EC cone vertex coordinates
c       calculated in subroutine call antenna_vertex
c
c       icone_loc the number of EC cone. It is used in output.f
c       
c       raypatt ='genray' or 'torray' : different specification
c                of the cenral ray
c               ='diskdisk' disk to disk launch
c--------------------------------------------------------------------
c       the data for disk to disk launching at raypat='diskdisk'
c
c       n_mesh_disk_mesh_radial_bin is the number of radial bins
c                              at the launching disk
c       n_mesh_disk_angle_bin() are the numbers of angle bins 
c                               inside each radial bin
c                         on the launching disk
c       sigma_launching_disk [m]
c       part_gauss_power   persent of integral from gaussian distribution
c                          at first disk with radis  rho_launching_disk
c       rho_launching_disk is the radius of the launching disk
c       rho_focus_disk     is the radius of the second (focus)disk,
c       d_disk             is the distance between two disks
c       initial_azimuth_angle_degree(n_mesh_radial_bin_a) are angles at each
c                          radial bins
c       disk_power_equal((n_mesh_radial_bin_a) for disk to disk launch at 
c                                     raypatt='diskdisk' case
c                                     It is the power for rays
c                                     launched from different
c                                     radial bins on the disk
c       power_per_ray                 power per ray launched from first disk
c
c       disk_rho_mesh_edge(n_mesh_radial_bin_a+1) non-uniform radial mesh
c                                                 of the bin's edge points
c                                                 at the first disk
c       disk_rho_mesh_center(n_mesh_radial_bin_a) nnon-uniform radial mesh
c                                                 of the bin's center points
c                                                 at the first disk
c
c       ar_z_disk_launch(nraymax),          !space coordinates of
c       ar_r_disk_launch(nraymax),          !launch point at the
c       ar_phi_disk_launch(nraymax),        !first disk
c       ar_beta_disk_launch(nraymax),       !poloidal and toroidal angles
c       ar_alpha_disk_launch(nraymax)       !give the ray direction

