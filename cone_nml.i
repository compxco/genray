c     The data for EC cone starting conditions, obtained from namelist
c     and subroutine dinit_mr calls.
      real*8
c-----namelist /eccone/
     &            zst,rst,phist, 
     &            alfast,betast,
     &            alpha1,alpha2,
     &            powtot,
     &            cr,
     &            rho_launching_disk,sigma_launching_disk,
     &            part_gauss_power, 
     &            rho_focus_disk,d_disk,
     &            initial_azimuth_angle_degree
    
      integer
c------namelist /eccone/
     &            ncone,
     &            n_mesh_disk_radial_bin,
     &            n_mesh_disk_angle_bin,
     &            gzone,mray,nray_in,
     &            na1,na2 
     
      character*8 
c------namelist /eccone/
     &            raypatt

       common/cone_nml/
c-----real*8
     &            zst(nconea),rst(nconea),phist(nconea), 
     &            alfast(nconea),betast(nconea),
     &            alpha1(nconea),alpha2(nconea),
     &            powtot(nconea),
     &            cr(gzonemax),
     &            rho_launching_disk,sigma_launching_disk,
     &            part_gauss_power, 
     &            rho_focus_disk,d_disk,
     &           initial_azimuth_angle_degree(n_mesh_disk_radial_bin_a),
c-----integer
     &            ncone,
     &            n_mesh_disk_radial_bin,
     &            n_mesh_disk_angle_bin(n_mesh_disk_radial_bin_a),
     &            gzone,mray(gzonemax),nray_in,
     &            na1,na2,    
c-----character
     &            raypatt
