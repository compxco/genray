
 
c        ********************* GENRAY ***********************
c        *                     ------                       *
c        * GENRAY is the computer code  for  obtaining  the *
c        * trajectories of different  wave  modes  by  ray- *
c        * tracing techniques                               *
c        *                                                  *
c        * Authors: Alexander P. Smirnov (primary)          *
c        *          Moscow State University                 *
c        *          (sap@ns.cnt.ru, sap@cs.msu.su)          *
c        *          R.W. Harvey          (secondary)        *
c        *          CompX                                   *
c        *          (bobh@compxco.com)                      *
c        *                                                  *
c        * Manual: GENRAY, Report CompX-01-2000 (2000)      *
c        *                 CompX, PO Box 2672, Del Mar, CA  * 
c        *                 www.compxco.com/Genray_manual.pdf* 
c        *                                                  *
c        ****************************************************

c***********************************************************************
c
c   Copyright A.P. Smirnov and R.W. Harvey
c   CompX company, Del Mar, California, USA
c   1995-2013
c   Below, "this program" refers to the GENRAY, and associated 
c   source files and manuals.
c
c   The primary references for the code are:
c   (1) A.P. Smirnov and R.W. Harvey, Calculations of the Current Drive
c   in DIII-D with the GENRAY Ray Tracing Code, Bull. Amer. Phys. Soc.
c   Vol. 40, No. 11, p. 1837, Abstract 8P35 (1995).
c   (2) A.P. Smirnov (Moscow State University) & R.W. Harvey (CompX)
c   CompX, P.O. Box 2672, Del Mar, CA 92014 Report CompX-2000-01, Ver. 2,
c   March 17, 2003, Available at http://www.compxco.com/genray.html
c   as Genray Manual (PDF).
c   
c
c                GNU Public License Distribution Only
c   This program is free software; you can redistribute it and/or modify
c   it under the terms of the GNU General Public License as published by
c   the Free Software Foundation; either version 3 of the License, or
c   any later version (and at end of this file).
c
c   This program is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c   GNU General Public License for more details.
c
c   You should have received a copy of the GNU General Public License
c   along with this program; if not, see <http://www.gnu.org/licenses/>.
c
c         --R.W. Harvey, CompX, Del Mar, CA, USA
c
c   E-mail:  rwharvey@compxco.com
c   Address: CompX company
c            P.O. Box 2672
c            Del Mar, CA 92014-5672
c
c   It will be appreciated by CompX if useful additions to GENRAY can be
c   transmitted to the authors, for inclusion in the source distribution.
c
c***********************************************************************
c***********************************************************************
c
c
c
c
c
c
c
c
c-----------------------------------------------------------------!
c                                                  !
c           method of solution and coordinate system.              !
c                                                  !
c        the wave trajectories are obtained from the solution     !
c        of  geometrical  optics  equations  where   independent  !
c        space variables are:                                 !
c                                                  !
c           z, r, phi - cylindrical coordinates                    !
c                                                  !
c        and canonically conjugate momenta are:                      !
c                                                  !
c           k_z, k_r, r*k_phi.                                !
c                                                  !
c        Geometrical optics equation in these  variables  are     !
c        hamiltonian in form and they are  solved  by  4th-order  !
c        runge-kutta method.                                !
c        The code maintains conservation of the hamiltonian       ! 
c        function with given accuracy eps by using two different  !
c        numerical methods :                                      !
c        for isolv=1 - the solution of six ray equations          !
c         with correction of coordinates at each time step.       ! 
c        for isolv=2                                              !
c         the solution of the dispersion equation for one ray        !
c         variable and the determination the other 5 variables        !
c         from 5 geometric optics equations. The code determines  !
c         automatically which ray variable which must be          !
c         determined from dispertion relation.                    !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c                                                                 !
c        the code is written in fortran and uses open source      !
c        libraries libnetcdf and mpich (for MPI compilations).    !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c                                                                 !
c       version                                                   !
c                                                                 !
c-----------------------------------------------------------------!
c       it uses input files: equilib.dat, genray.in (genray.dat)--!


CMPIINSERTPOSITION PROGRAMSTART


      PROGRAM GENRAY

c      implicit double precision (a-h,o-z)
      implicit none

      include 'param.i' ! specifies code parameters 
      include 'commons.i'
c      include 'transport_prof.i'  [added to commons.i, SAP080629]

CMPIINSERTPOSITION DECLARATION
      
      integer qsize ! In serial run: qsize=1
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
cSm070718      dimension u(6),deru(6),aux(8,6)
      real*8  u(6),deru(6),aux(8,6),
     +energy,pitch,fdist,dfdx,dfdpitch,dfdp, !for call dskin
     +r,z,phi,fdens,fdens0,fmaxw,
     +zst1,rst1,phist1,alfast1,betast1,
     +cnteta,cnphi,      
     +win_0dsmax,ds,r0_em,z0_em,
     +psi_loc,rho_loc,clight,sum_emission,sum_emission_wall


      integer iraystop,nray,ndim,
     +initial,ihlf,n_r,n_z,ifreq,
     +i_bad_initial_conditions,n,
     +nrayelt_o_cutoff_emis,n0,igenray,
     +i_geom_optic_loc

c-----externals
      external outpt
      external rside1,rsideb1,outptb1,b,
     +fdens_fdist,fdens0_fdist,fdens_fkin2,
     +length_char,tempe,temperho,fpsi,rhopsi,r_2nd_harm,
     +dvarddz,vardens,
     +ddnsrho,densrho

      real*8 thetapol,psif,drhodz,dthetadz,dthetadr,drhodr,ddnsrho,b1
      real*8 b,dense,density_r_z_i,d_density_r_z_i_d_r,x,dxdr,dxdz,
     +d_density_r_z_i_d_z,drhopsi,
     +fdens_fdist,fdens0_fdist,fdens_fkin2,tempe,temperho,fpsi,rhopsi,
     +r_2nd_harm,zeffrho,
     +dense_no_RZ_spline

      integer length_char

      real*4 time_loop_ifreq_1,time_loop_ifreq_2,time_loop_ifreq,
     +time_before_rk,time_after_rk,time_rk,
     +time_drkgs2_1,time_drkgs2_2,
     +time_genray_1,time_genray_1a,time_before_1st_ray,time_genray_2,
     +time_emission_2,
     +time_emission_1
c------------------------------------------------------------
cfor test CD_adj_LH_efficiency
      real*8 lh_cd_efficiency,cnpar,cnper,thetapol_l,
     +u_ph_karney,u1,efficien,temp_kev,unorm

      integer n_radial
cfor interpolation_chi
      real*8 u0,theta0,
     +chi,d_chi_d_u_out,d_chi_d_theta_out
      integer n_radial0
cSAP080711
c-----genray input file: genray.in or genray.dat
      character*10 genray_in_dat

cSAP081122
      real*8 
     +length_b

cSAP090205 to check subroutine sigma_edge_n_theta_pol
      integer i,i0r,j0r,i0z,j0z

      real*8  theta_pol_radian,
     +sigma_edge_n,d_sigma_edge_n_d_theta_pol,
     +d2_sigma_edge_n_d2_theta_pol

ctest for spline density art rz mesh
      real*8 dens_rho_l,dens_rz,d_norm,diff_dens,
     +dif_d_dens_dr,dif_d_dens_dz,
     +d_norm_r,d_norm_z,d_dens_rho_r,d_dens_rho_z,
     +d_dens_spl_r,d_dens_spl_z,
     +dens_rho_theta,d_dens_rho_theta_d_rho,
     +d_dens_rho_theta_d_theta,
     +dn_drho,dro_dpsi,
     +step_rz,d_p,d_m,deriv_l,r_p,r_m,z_p,z_m,
     +x_l,x_p,x_m,dxdr_l,dxdz_l,dxdr_d,dxdz_d,
     +rho_z_p,rho_z_m,rho_r_p,rho_r_m,     
     +dens_z_p,dens_z_m,dens_r_p,dens_r_m
      integer j,k,m

c for test lsc approach
      real*8 hrho,h_v
      
      data myrank/0/  !In serial run: myrank=0; In MPI run: myrank=rank
      data qsize/1/   !In serial run: qsize=1

CMPIINSERTPOSITION INITIALIZATION      
          
      if (myrank.eq.0) then
         call cpu_time(time_genray_1) ! very beginning
      endif  !On myrank.eq.0 

c---------------------------------------------------
c     Write out the version number
c--------------------------------------------------

      write(t_,1000) version
 1000 format(15x,"GENRAY VERSION: ",a)
      if (myrank.eq.0) write(6,'(//16x,
     +"=================================================")')
      if (length_char(t_) .gt. 256) then
         if (myrank.eq.0) then
            WRITE(*,*)'GENRAY:Adjust length of t_'
         endif
         stop ! stop in all cores
      endif
      if (myrank.eq.0) write(6,*) t_(1:length_char(t_))
      if (myrank.eq.0) write(6,'(16x,
     +"=================================================",/)')
c---------------------------------------------------
c     Write out all code parameters given in param.i
c--------------------------------------------------
      call screen_print_all_parameters ! called from rank=0 only

c---------------------------------------------------
c     check the parameters in param.i
c--------------------------------------------------
      if (nrya.ne.(max(nxeqda,nyeqda)+4)) then
      if (myrank.eq.0) then
         WRITE(*,*)'it should be nrya.eq.(max(nxeqda,nyeqda)+4)'
         WRITE(*,*)'but nrya,nxeqda,nyeqda', nrya,nxeqda,nyeqd
         WRITE(*,*)'change the parameter nrya in param.i'
      endif  !On myrank.eq.0 
      stop ! stop in all cores
      endif
      
CMPIINSERTPOSITION STARTTIME



c      if (nzy.ne.(max(npsi,nteta1)+4)) then
c         write(*,*)'it should be nzy.eq.(max(npsi,nteta1)+4)'
c         write(*,*)'but nzy,npsi,nteta1', nzy,npsi,nteta1
c         write(*,*)'change the parameter nzy in param.i'
c         stop
c      endif

c      if (nzy.ne.(max0(npsi4,nteta1)+4)) then 
      if (nzy.ne.(max0(npsi,nteta1)+4)) then
         if(myrank.eq.0) then
c        write(*,*)'in param.i nzy.ne.max0(npsi4,nteta1)+4'
         WRITE(*,*)'in param.i nzy.ne.max0(npsi,nteta1)+4'
         WRITE(*,*)'nzy= ',nzy
c        write(*,*)'max0(npsi4,nteta1)',max0(npsi4,nteta1)
         WRITE(*,*)'max0(npsi,nteta1)=4',max0(npsi,nteta1)+4
         WRITE(*,*)'Change nzy in param.i'
         endif
         stop
      endif

c      if (nwka.ne.(1+3*max(jxa,iya))) then
c--------the parameters for the distribution function spline aproximation   
c         write(*,*)'it should be nwka.eq.(1+3*max(jxa,,iya))'
c         write(*,*)'but nwka,jxa,iya',nwka,jxa,iya
c         write(*,*)'change the parameter nwka in param.i'
c         stop
c      endif

      
c-----the creation of default input data like in genray.dat file
      call default_in
CWRITE       write(*,*)'genray.f after default_in'
CWRITE       write(*,*)'i_resonance_curve_integration_method=',
CWRITE      +i_resonance_curve_integration_method

c-----check if input file is genray.in, then it will change
c     default data to genray.in file (MKSA) format  
      call transform_input_data_to_MKSA
     
c-----Initialize the value of nj, in case read_transport_prof not used.
      nj=0

c---- read all namelists from input genray.in or genray.dat file
      call read_all_namelists(genray_in_dat,ndim,nray)
      if(myrank.eq.0) then
         WRITE(*,*)'genray.f after read_all_namelists ndim ',ndim
         WRITE(*,*)'in genray.f genray_in_dat=', genray_in_dat
         WRITE(*,*)'genray.f after read_all_ prmt ',prmt
         WRITE(*,*)'genray.f after read_all_ ioxm',ioxm   
         WRITE(*,*)'genray.f after read_all_ nray',nray  
         WRITE(*,*)'genray.f after read_all_namelists data in /genr/'
         WRITE(*,*)'r0x=',r0x
         WRITE(*,*)'b0=',b0
         WRITE(*,*)'outdat=',outdat
         WRITE(*,*)'stat=',stat
         WRITE(*,*)'mnemonic=', trim(mnemonic)
         WRITE(*,*)'rayop=',rayop
         WRITE(*,*)'dielectric_op=',dielectric_op
         WRITE(*,*)'partner=',partner
      endif ! myrank=0
      
      if(ilaunch.eq.2) then
CWRITE       write(*,*)'ngrill',ngrill
      do i=1,ngrill
CWRITE       write(*,*)'i,nthin(i)',i,nthin(i)
      do j=1,nthin(i)
CWRITE       write(*,*)'j,rlaunch(j,i),zlaunch(j,i)',
CWRITE      +j,rlaunch(j,i),zlaunch(j,i)
CWRITE       write(*,*)'weight_power_launch(j,i)',
CWRITE      +weight_power_launch(j,i)
      enddo
      enddo
      endif

c-----allocate data 
      call ainalloc

c-----test sigma_edge_n
      pi=4*datan(1.d0)
CWRITE       write(*,*)' i_edge_dens_anal',i_edge_dens_anal
      if(i_edge_dens_anal.eq.2) then
c--------spline data for  sigmedgn_ar
      do i=1,n_pol_edge_dens
      theta_pol_radian=2.d0*pi/( n_pol_edge_dens-1)*(i-1)
      call sigma_edge_n_theta_pol(theta_pol_radian,
     +sigma_edge_n,d_sigma_edge_n_d_theta_pol,
     +d2_sigma_edge_n_d2_theta_pol)
CWRITE       write(*,*)'i,theta_pol_radian,sigma_edge_n',
CWRITE      +i,theta_pol_radian,sigma_edge_n
CWRITE       write(*,*)'d_sigma_edge_n_d_theta_pol',
CWRITE      +d_sigma_edge_n_d_theta_pol
CWRITE       write(*,*)'d2_sigma_edge_n_d2_theta_pol',
CWRITE      +d2_sigma_edge_n_d2_theta_pol
      enddo
      endif

c-----If input file is genray.in, then it will change
c     input data to genray.dat file format 
      if (genray_in_dat.eq.'genray.in')  then
CWRITE       write(*,*)'genray.f before transform_genray_in_to_dat'
      call transform_genray_in_to_dat 
      endif
CWRITE       write(*,*)'genray.f after tramnsform_genray ',prmt
c---------------------------------------------------------------
c     reading the eqdsk data
c     and creation of the coefficients for spline approximation
c     for: feqd, psi -magnetic field functions
c          zpllim(r),zminlim (r) - limmiter boundary
c---------------------------------------------------------------
CWRITE       write(*,*)' genray.f before call equilib'
      if (myrank.eq.0) then
      call cpu_time(time_genray_1a) !read input,allocation; just before equilib
      endif  !On myrank.eq.0 

      call equilib
CWRITE       write(*,*)' genray.f after call equilib'
CWRITE       write(*,*)'NR',NR
c---------------------------------------------------------------
CWRITE       write(*,*) 'genray.f before call rhospl'
      call rhospl     
CWRITE       write(*,*) 'genray.f after call rhospl'
CWRITE       write(*,*)'arpsi',arpsi
c----------------------------------------------------
c    compare poloidal and B field lengths 
c      do n=1,npsi     
       
c        call calc_length_b(arpsi(n),length_b)
c        write(*,*)'n,arpsi(n), rhopsi(arpsi(n))',
c     &             n,arpsi(n), rhopsi(arpsi(n))
c         write(*,*)'pollength,length_b',arrho_l(n)*totlength,length_b
c         if (n.gt.1) write(*,*)'length_b/pollength',
c     &                          length_b/(arrho_l(n)*totlength)
     
c      enddo
c      stop 'genray.f calc_length_b'
c--------------------------------------------------------------------
c     reading the name of the output file
c--------------------------------------------------------------------
      i1_=91
      if(myrank.eq.0)then ! MPI write from myrank=0 only ! 
         open(i1_,file=outdat)
      endif ! myrank=0
      
      
      iraystop=0
c-------------------------------------------------------------------
c     reading input data from genray.in file
CWRITE       write(*,*) 'genray.f before call dinit_mr'
CWRITE       write(*,*)'NR=',NR
      call dinit_mr(ndim,nray)
      if(myrank.eq.0) then
         WRITE(*,*)'genray.f after dinit_mr ndim ',ndim
         WRITE(*,*)'genray.f after dinit_mr nray',nray  
         !pause
      endif ! myrank=0
c-----allocate pointers at writencdf.i
c      write(*,*)'genray.f before  ainalloc_writencdf nray',nray
c      call ainalloc_writencdf_i(nray)
c      call ainalloc_write_i(nray)

CWRITE       write(*,*)'genray.f after dinit_mr prmt=',prmt
CWRITE       write(*,*)'genray !!!!!ndim',ndim
CWRITE       write(*,*)'genray after dinit_mr nray=',nray
CWRITE       write(*,*)'genray after dinit_mr nbulk',nbulk
CWRITE       write(*,*)'v',v
CWRITE       write(*,*)'w',w
CWRITE       write(*,*)'genray.f after dinit_mr ioxm ',ioxm
CWRITE       write(*,*)'genray.f after dinit_mr freqncy0',freqncy0
      do iray=1,nray         
CWRITE       write(*,*)'1 iray,arzu0(iray),arru0(iray),arphiu0(iray)',
CWRITE      +iray,arzu0(iray),arru0(iray),arphiu0(iray)
      enddo
CWRITE       write(*,*)'genray.f before n_wall.eq.1  n_wall',n_wall 

      if (n_wall.gt.1)then
c------------------------------------------------------------------
c        create additional points at the chamber wall
c-----------------------------------------------------------------
c        interpolate wall points at the wall mesh with additinal ponts 
c        r_wall_add(n_wall_add) z_wall_add(n_wall_add)
c        and calculate the number of points of this mesh: n_wall_add)
c        Result will be in fourbe.i
      call ainalloc_fourb_i !allocate arrays for fourb.i
      call create_fine_mesh_for_chamber_wall_limiter_coordinates
      call add_horizontal_limiter_walls

      if(myrank.eq.0)then
         call MK_graph_chamber_wall
      endif ! myrank=0

c         stop 'genray,f after  MK_graph_chamber_wall'
c--------------------------------------------------------------------
c        create RZ meshes: rr_add(nxeqd_add),zz_add(nyeqd_add)
c-------------------------------------------------------------------
      call creat_fine_2D_poloidal_mesh
c--------------------------------------------------------------------
c        calculate distance from RZ mesh points to chamber wall:
c        distance_to_wall(nxeqd_add,nyeqd_add)
c-------------------------------------------------------------------
CWRITE       write(*,*)'genray.f before distance_for_wall'
      call distance_for_wall  
CWRITE       write(*,*)'genray.f after distance_for_wall'
c--------test fast subroutine: it gave wrong result 
c         call distance_for_wall_1
c         stop 'genray.f after distance_for_wall_1'

      if (i_edge_dens_rz_mesh.ne.2)then
c-------------------------------------------------------------------
c           Calculate density array:  
c           density_r_z(nxeqd_add,nyeqd_add,nbulk)
c           at (rr_add,zz_add) mesh
c           For( n_wall.gt.1) case create density fall near the chamber wall.
c-------------------------------------------------------------------
      call density_at_zr_plane  
CWRITE       write(*,*)'genray.f after density_at_zr_plane'

c            do k=1,nbulk
c              do j=1,nyeqd_add
c                do i=1,nxeqd_add
c                  write(*,*)'0k,j,i,density_r_z(i,j,k,0)',
c     &                        k,j,i,density_r_z(i,j,k,0)
c                enddo
c              enddo
c            enddo

      endif          !i_edge_dens_rz_mesh.ne.2
c-----------------------------------------------------------------------
c        write dens_rz_in file to dens_rz_out.txt for test  
c------------------------------------------------------------------------ 
      if (i_edge_dens_rz_mesh.eq.1) then
c-----------------------------------------------------------------------
c           write arays
c           dens_r_z_in(i,j,k),temperature_r_z_in(i,j, zeff_r_z_in(i,j,k)
c           and poloidal flux and small radius values in rz mesh points
c           to the output file dens_rz_out.txt
c---------------------------------------------------------------------
cSAP111002

      if(myrank.eq.0)then ! MPI write from myrank=0 only !
          call writes_dens_temp_psi_rho_rz_out_txt 
      endif ! myrank=0
      
CWRITE       write(*,*)'genray.f writes_dens_temp_psi_rho_rz_out_txt '
      endif  ! i_edge_dens_rz_mesh.eq.1) 
c-----------------------------------------------------------------------
      if (i_edge_dens_rz_mesh.eq.2) then
cSAP111002
c---------------------------------------------------------
c           reads density in [MKS] at dens_rz_in(i,j,k) at rz mesh
c           from the input file dens_rz_in.txt
c           This density will be used outside LCFS
c-------------------------------------------------------- 
      call read_density_temperature_at_rz_mesh_txt
c---------------------------------------------------------------------
c           put density dens_r_z_in,temperature_r_z_in,
c           given in the input file dens_temp_rz_in.dat
c           to arrays density_r_z, temperature_r_z
c--------------------------------------------------------------------
      call
     +put_dens_temp_rz_into_dens_temp_r_z
c           do k=1,nbulk
c              do j=1,nyeqd_add
c                do i=1,nxeqd_add
c                  write(*,*)'2 k,j,i,density_r_z(i,j,k,0)',
c     &                         k,j,i,density_r_z(i,j,k,0)
c                enddo
c              enddo
c            enddo
c-------------------------------------------------------------------
c           calculate temperature spline coefficients 
c           at RZ mesh using array temperature_r_z
c------------------------------------------------------------------
      call splcoef_temperature_r_z
c-------------------------------------------------------------------
c           calculate z_eff spline coefficients 
c           at RZ mesh using array zeff_r_z
c------------------------------------------------------------------
c           write(*,*)'genray.f before splcoef_zeff_r_z'
      call splcoef_zeff_r_z
      endif  ! i_edge_dens_rz_mesh.eq.2
c-------------------------------------------------------------------
c        calculate density spline coefficients 
c        at RZ mesh using array density_r_z
c-------------------------------------------------------------------
CWRITE       write(*,*)'genray.f before splcoef_density_r_z'
      call splcoef_density_r_z
CWRITE       write(*,*)'genray.f after splcoef_density_r_z'
c         stop 'genray.f after splcoef_density_r_z'

      endif

      goto 17
c------test of RZ density spline      for n_wall >0
      if (n_wall.ge.1) then

      phi=0.d0
      phi=1.d0
      m=0
      m=1
      do k=1,nbulk
      d_norm=0.d0
      d_norm_r=0.d0  
      d_norm_z=0.d0
      dif_d_dens_dr=0.d0
      dif_d_dens_dz=0.d0
      i0r=0
      j0r=0
      i0z=0
      j0z=0
      do i=1,nxeqd_add
      do j=1,nyeqd_add
CWRITE       write(*,*)'k,i,j,zz_add(j),rr_add(i)',
CWRITE      +k,i,j,zz_add(j),rr_add(i)
      z=zz_add(j)
      r=rr_add(i)
      bmod=b(z,r,phi) !calculate rho

c             write(*,*)'genray.f before dense_no_RZ_spline'
      dens_rho_l=dense_no_RZ_spline(z,r,phi,k)
CWRITE       write(*,*)'genray.f before dense dens_rho_l',dens_rho_l

      dens_rz=dense(z,r,phi,k)
CWRITE       write(*,*)'genray.f after dense dens_rz', dens_rz

      if(dabs(distance_to_wall(i,j,m)).le.1.d-3) then
CWRITE       write(*,*)'i,j,z,r,rhodistance_to_wall(i,j,m)',
CWRITE      +i,j,z,r,rho,distance_to_wall(i,j,m)
CWRITE       write(*,*)'dens_rho_l,dens_rz',dens_rho_l,dens_rz
      endif
 
      d_norm=d_norm+dabs(dens_rho_l-dens_rz)
      diff_dens=dens_rho_l-dens_rz
            
c------------derivatives
cSAP090408
c             write(*,*)'derivatives rho',rho

      if (rho.gt.1.d0-1.d-10) then
      thetapol_l=thetapol(z,r) ! -pi <thetapol =<pi

      if (thetapol_l.lt.0d0) then
      thetapol_l=thetapol_l+2*pi !pi< theta_pol<2pi
      endif
cSAP090408
c                write(*,*)'before dens_rho_theta_LCFS'

      call dens_rho_theta_LCFS(rho,thetapol_l,k,
     +dens_rho_theta,d_dens_rho_theta_d_rho,
     +d_dens_rho_theta_d_theta)

cSAP090408
c                write(*,*)'after dens_rho_theta_LCFS'

      d_dens_rho_r=d_dens_rho_theta_d_rho*drhodr(z,r,phi)
     ++d_dens_rho_theta_d_theta*dthetadr(z,r)

      d_dens_rho_z=d_dens_rho_theta_d_rho*drhodz(z,r,phi)
     ++d_dens_rho_theta_d_theta* dthetadz(z,r)
cSAP090408
CWRITE       write(*,*)'d_dens_rho_r,d_dens_rho_z',
CWRITE      +d_dens_rho_r,d_dens_rho_z

      goto 15
c---------------numerical derivative d_dens_rho_theta_d_rho
      step_rz=1.d-4
      step_rz=1.d-6
      call dens_rho_theta_LCFS(rho+step_rz,thetapol_l,k,
     +x_p,d_dens_rho_theta_d_rho,
     +d_dens_rho_theta_d_theta)
      call dens_rho_theta_LCFS(rho-step_rz,thetapol_l,k,
     +x_m,d_dens_rho_theta_d_rho,
     +d_dens_rho_theta_d_theta)
      dxdz_d=(x_p-x_m)/(2.d0*step_rz) !d_dens_drho numerical
CWRITE       write(*,*)'d_dens_rho_theta_d_rho,dxdz_d',
CWRITE      +d_dens_rho_theta_d_rho,dxdz_d
c---------------numerical derivative d_dens_rho_theta_d_theta
      call dens_rho_theta_LCFS(rho,thetapol_l+step_rz,k,
     +x_p,d_dens_rho_theta_d_rho,
     +d_dens_rho_theta_d_theta)
      call dens_rho_theta_LCFS(rho,thetapol_l-step_rz,k,
     +x_m,d_dens_rho_theta_d_rho,
     +d_dens_rho_theta_d_theta)
      dxdz_d=(x_p-x_m)/(2.d0*step_rz) !d_dens_drho numerical
CWRITE       write(*,*)'d_dens_rho_theta_d_theta,dxdz_d',
CWRITE      +d_dens_rho_theta_d_theta,dxdz_d
c---------------numerical derivatives drhodr(z,r,0.d0)-----------------
      bmod=b(z,r+step_rz,k) 
      x_p=rho
      bmod=b(z,r-step_rz,k)
      x_m=rho
      dxdr_d=(x_p-x_m)/(2.d0*step_rz) !d_rho/d_r numerical
CWRITE       write(*,*)'drhodr(z,r,phi),dxdr_d',
CWRITE      +drhodr(z,r,phi),dxdr_d
c---------------numerical derivatives drhodz(z,r,0.d0)-----------------
      bmod=b(z+step_rz,r,k) 
      x_p=rho
      bmod=b(z-step_rz,r,k)
      x_m=rho
      dxdz_d=(x_p-x_m)/(2.d0*step_rz) !d_rho/d_z numerical 
CWRITE       write(*,*)'drhodz(z,r,phi),dxdz_d',
CWRITE      +drhodz(z,r,phi),dxdz_d
c----------------------------------------------------------------
      else
      psi_loc=psif(z,r)
      dro_dpsi=drhopsi(psi_loc)
      dn_drho=ddnsrho(rho,k)
      d_dens_rho_r=dn_drho*dro_dpsi*dpdrd
      d_dens_rho_z=dn_drho*dro_dpsi*dpdzd
      endif
c------------density derivatives using spline
 15   continue
cSAP090409
c             write(*,*)'genray.f before d_density_r_z_i_d_r k,m',k,m
      d_dens_spl_r=d_density_r_z_i_d_r(z,r,phi,k)
      d_dens_spl_z=d_density_r_z_i_d_z(z,r,phi,k)
      goto 16
c             step_rz=1.d-5
      step_rz=1.d-4
c-------------numerical d_dens_d_z from spline rz density            
      x_p=density_r_z_i(z+step_rz,r,phi,k)
      x_m=density_r_z_i(z-step_rz,r,phi,k)
      dxdz_d=(x_p-x_m)/(2.d0*step_rz)

CWRITE       write(*,*)'numer der z spl'
CWRITE       write(*,*)'x_p,z+step_rz,x_m,z-step_rz,dxdz_d',
CWRITE      +x_p,z+step_rz,x_m,z-step_rz,dxdz_d

      x_p=density_r_z_i(z,r+step_rz,phi,k)
      x_m=density_r_z_i(z,r-step_rz,phi,k)
      dxdr_d=(x_p-x_m)/(2.d0*step_rz)

CWRITE       write(*,*)'numer der r spl'
CWRITE       write(*,*)'x_p,r+step_rz,x_m,r-step_rz,dxdr_d',
CWRITE      +x_p,r+step_rz,x_m,r-step_rz,dxdr_d

CWRITE       write(*,*)'num der spl dxdz_d,dxdr_d',dxdz_d,dxdr_d

c-------------numerical d_dens_d_z from function dense             
      bmod=b(z+step_rz,r,k) 
      x_p=dense(z+step_rz,r,phi,k)
      bmod=b(z-step_rz,r,k) 
      x_m=dense(z-step_rz,r,phi,k)
      dxdz_d=(x_p-x_m)/(2.d0*step_rz)

CWRITE       write(*,*)'numer der z dense'
CWRITE       write(*,*)'x_p,z+step_rz,x_m,z-step_rz,dxdz_d',
CWRITE      +x_p,z+step_rz,x_m,z-step_rz,dxdz_d

      bmod=b(z,r+step_rz,k) 
      x_p=dense(z,r+step_rz,phi,k)
      bmod=b(z,r-step_rz,k) 
      x_m=dense(z,r-step_rz,phi,k)
      dxdr_d=(x_p-x_m)/(2.d0*step_rz)

CWRITE       write(*,*)'numer der r dense'
CWRITE       write(*,*)'x_p,r+step_rz,x_m,r-step_rz,dxdr_d',
CWRITE      +x_p,r+step_rz,x_m,r-step_rz,dxdr_d


CWRITE       write(*,*)'num der dense dxdz_d,dxdr_d',dxdz_d,dxdr_d
c---------------------------------------------------------

 16   continue
      if (dif_d_dens_dr.lt.dabs(d_dens_rho_r-d_dens_spl_r)) then
      i0r=i
      j0r=j
      dif_d_dens_dr=dabs(d_dens_rho_r-d_dens_spl_r)
      endif

      if (dif_d_dens_dz.lt.dabs(d_dens_rho_z-d_dens_spl_z)) then
      i0z=i
      j0z=j
      dif_d_dens_dz=dabs(d_dens_rho_z-d_dens_spl_z)
      endif

      d_norm_r= d_norm_r+dabs(d_dens_rho_r-d_dens_spl_r)
      d_norm_z= d_norm_z+dabs(d_dens_rho_z-d_dens_spl_z)

      if (dabs(distance_to_wall(i,j,0)).le.1.d-3) then
CWRITE       write(*,*)'d_dens_rho_r,d_dens_spl_r',d_dens_rho_r,d_dens_spl_r
CWRITE       write(*,*)'d_dens_rho_z,d_dens_spl_z',d_dens_rho_z,d_dens_spl_z
      endif

      enddo
      enddo
CWRITE       write(*,*)'genray.f k d_norm,d_norm_r,d_norm_z',
CWRITE      +k,d_norm,d_norm_r,d_norm_z
CWRITE       write(*,*)'genray.f dif_d_dens_dr,i0r,j0r',
CWRITE      +dif_d_dens_dr,i0r,j0r
CWRITE       write(*,*)'genray.f dif_d_dens_dz,i0r,j0z',
CWRITE      +dif_d_dens_dz,i0z,j0z


      enddo
 
      endif !nwall>0
c------end test
 17   continue
c       stop 'genray.f after compare splcoef_density_r_z'

c-------------------------------------------------
c     compare two different splines for b calculations
c      call  test_b_field
c      stop 'genray.f after  test_b_field'
c-------------------------------------------------------------------
c     test of non-maxwellian distribution
c     with the 3D distribution function f written by CQL3D
c      do iray=1,nray         
c       write(*,*)'2 iray,arzu0(iray),arru0(iray),arphiu0(iray)',
c     +iray,arzu0(iray),arru0(iray),arphiu0(iray)
c      enddo
      goto 30

CWRITE       write(*,*)'genray.f i_diskf',i_diskf
     
      if (i_diskf.ne.0) then
CWRITE       write(*,*)'i_diskf= ',i_diskf
ctest
      initial=1
      call dskin(initial,energy,pitch,rho,
     +fdist,dfdx,dfdpitch,dfdp,1)
CWRITE       write(*,*)'diskf file was written'
c
      
c         call test_fdist

      r=1.50d0
      z=0.0d0
      phi=0.d0
CWRITE       write(*,*)'!!!!!!! -1'
      fdens0=fdens0_fdist(r,z,phi)
CWRITE       write(*,*)'in genray r,z,fdens0_fdist',r,z,fdens0 
      fdens=fdens_fdist(r,z,phi)
CWRITE       write(*,*)'in genray r,z,fdens_fdist',r,z,fdens

      r=1.20d0
      z=0.0d0
CWRITE       write(*,*)'!!!!!!! -2'
      fdens0=fdens0_fdist(r,z,phi)
CWRITE       write(*,*)'in genray r,z,fdens0_fdist',r,z,fdens0 
      fdens=fdens_fdist(r,z,phi)
CWRITE       write(*,*)'in genray r,z,fdens_fdist',r,z,fdens
 

      r=2.0d0
      z=0.1d0
CWRITE       write(*,*)'!!!!!!! -3'
      fdens0=fdens0_fdist(r,z,phi)
CWRITE       write(*,*)'in genray r,z,fdens0_fdist',r,z,fdens0 
      fdens=fdens_fdist(r,z,phi)
CWRITE       write(*,*)'in genray r,z,fdens_fdist',r,z,fdens
 

      fdens=fdens_fkin2(1,1,r,z)
CWRITE       write(*,*)'in genray fdens_fkin2(1,5,r,z)',fdens 
         
ctest
      energy=0.5d0
      pitch=0.5d0 
      rho=0.5d0
CWRITE       write(*,*)'!!!!energy,pitch,rho',energy,pitch,rho
      call dskin(0,energy,pitch,rho,
     +fdist,dfdx,dfdpitch,dfdp,2)
CWRITE       write(*,*)'fdist,dfdx,dfdp,dfdpitch',
CWRITE      +fdist,dfdx,dfdp,dfdpitch
CWRITE       write(*,*) 'in genray after dskin before stop'
      call fmaxw_en(energy,rho,fmaxw,dfdx)
c         fmaxw=fmaxw_en(energy,rho)
CWRITE       write(*,*)'genray.f fmaxw,dfdx',fmaxw,dfdx 
         

c         stop
cendtest
      endif
 30   continue
c--------------------------------------------------------------
 
 
c-----------------------------------------------------------------------
c     For the plotting the output Toray data.
c     READs the data from toray (old_3d.dat) and creates genray.bin file     
c     call read3d(15)
c     write(*,*)'genray.f after read3d'
c     call MK_GRAPT (1,1)
c     write(*,*)'genray.f after mk_grapt'
c-------------------------------------------------------------


c-------------------------------------------------------------------
c      calculations of contours bmod(z,r)=const 1/y_i(z,r)=1,2,..
c      call contourb
c      write(*,*)'in genray before contrb2'
c      call contrb2
c      write(*,*)'write after contourb'
c-------------------------------------------------------------------
c     set zero to arrays power and current
c     for subroutine p_c_prof
c------------------------------------------
      do iray=1,nray         
CWRITE       write(*,*)'3 iray,arzu0(iray),arru0(iray),arphiu0(iray)',
CWRITE      +iray,arzu0(iray),arru0(iray),arphiu0(iray)
      enddo
CWRITE       write(*,*)'in genray ionetwo',ionetwo
      if(ionetwo.eq.1) then
      call onetwoini
      endif

c-------------------------------------------------------------------
c     creation of the circles on horizontal plane (X,Y) with
c     radiuses rmin,rmax,ra=xma (for tokamak toroidal section
c     graphic picture).Output:XT(3,NP+1),YT(3,NP+1) in common gr.cb
      call gr3
c------------------------------------------------------------------
      call output_con1

c---------------------------------------------------------------
c     the loop on all rays
c---------------------------------------------------------------
c---------------------------------------------------------------
c     Initialize arrays for Runge-Kutta subroutine
      call arrays(ndim,deru,prmt,ihlf)

c-----------------------------------------------------------------
      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
        if (isolv.eq.1) then
        write(*,*)'the runge-kutta solution of 6 hamiltonian equations
     +  with correction which gives hamiltonian conservation'
        write(*,*)'accuracy of the hamiltonian conservation
     +  in correction procedure epscor=',prmt(4)
        endif
        if (isolv.eq.2) then
        write(*,*)'the runge-kutta solution of 5 hamiltonian equations
     +  and the solution of the dispersion relation for the sixth
     +  ray variable'
        endif
      endif ! outprint

      if (nray.gt.nraya) then
         if(myrank.eq.0)then
         WRITE(*,*)'genray.f: nray=',nray,'nraya=',nraya
         WRITE(*,*)'nray>nraya. Inrease nraya in param.i'
         endif
         stop
      endif  
      if(nray.gt.nraymax)then
         if(myrank.eq.0)then
         WRITE(*,*)'genray.f: nray=',nray,'nraymax=',nraymax
         WRITE(*,*)'nray>nraymax. Increase nraymax in param.i'
         endif
         STOP ! at all cores
      endif
      if(ncone.gt.nconea)then
         if(myrank.eq.0)then
         WRITE(*,*)'genray.f: ncone=',ncone,'nconea=',nconea
         WRITE(*,*)'ncone>nconea. Increase nconea in param.i'
         endif
         STOP ! at all cores
      endif
      if(n_mesh_disk_radial_bin.gt.n_mesh_disk_radial_bin_a)then
         if(myrank.eq.0)then
         WRITE(*,*)'genray.f: n_mesh_disk_radial_bin=',
     +                        n_mesh_disk_radial_bin,
     +                       'n_mesh_disk_radial_bin_a=',
     +                        n_mesh_disk_radial_bin_a
         WRITE(*,*)'n_mesh_disk_radial_bin>n_mesh_disk_radial_bin_a. 
     +   Increase n_mesh_disk_radial_bin_a  in param.i'
         endif
         STOP ! at all cores
      endif

c
c-----Construct names of .txt and .nc ray data files
      if( length_char(trim(mnemonic)).gt.124)
     +stop 'Adjust mnemonic in genray.f'
      write(filetxt,1001) mnemonic(1:length_char(trim(mnemonic)))
 1001 format(a,".txt")
      write(filenc,1002) mnemonic(1:length_char(trim(mnemonic)))
 1002 format(a,".nc")
c
      if( myrank.eq.0 ) then !----------------------------myrank=0
      
      if(outnetcdf.eq.'enabled')then !YuP[2018-01-17] Added
c       open text file for 3d FP code
          if(rayop.eq."text" .or. rayop.eq."both") then
             i_=92
             open(i_,file=trim(filetxt))
          endif
c-----preparing data for output mnemonic.txt and mnemonic.nc file
CWRITE       write(*,*)'in genray.f before write3d1 nray',nray
          call write3d1(nray)      
CWRITE       write(*,*)'in genray.f sub after write3d1 nrayl=',nrayl
      endif ! outnetcdf
      
      call plotinit ! open PGplot

c-----plot magnetic field contours to plot.ps file
      n_r=100
      n_z=100
      if(i_plot_b.eq.1) then ! plotting b,n,t in plot.ps file
CWRITE       write(*,*)'genray.f before map_b'
CWRITE       write(*,*)'nbulk',nbulk
CWRITE       write(*,*)'v',v
CWRITE       write(*,*)'w',w 
CWRITE       write(*,*)'genray.f before map_b'
      call map_b(n_r,n_z)   
CWRITE       write(*,*)'genray.f after map_b'
c--------------------------------------------------------------------
c       plot frequencies along the stright line of length dist_freq [m]
c       with the edge point z_freq,r_freq [m]
c       and directed by angles alpha_freq,beta_freq [degree]
c------------------------------------------------------------------
c ARIES
c      z_freq=0.0d0
c      r_freq=5.7d0
c      dist_freq=4.7d0
c      alpha_freq=180.0d0
c NSTX
c      z_freq=0.0d0
c      r_freq=1.49d0
c      dist_freq=1.28d0 !the length of the stright line

c      alpha_freq=180.0d0
c      beta_freq=90.0d0
c      nsteps_freq=780
c      n_ec_harmonics_freq=6
c      max_plot_freq200.d0
      call plot_fcefuh(z_freq,r_freq,alpha_freq,beta_freq,dist_freq,
     +nsteps_freq,n_ec_harmonics_freq,npar_freq,max_plot_freq)
c         stop 'genray.f after plot_fcefuh'
      
c---------------------------------------------------------------
c     plot fig.1 for Ehst-Karney efficiency
c
c      write(*,*)'genray.f before plot_1_Karney'
c      call plot_1_Karney
c      write(*,*)'genray.f after plot_1_Karney'
c      stop
c      call plot_1_ADJ
c
c      See also, below:call plot_1_ADJ_Karney
c---------------------------------------------------------------
c      plot hot plasma dispersion function contours
c      D(ReN_perp,ImN_perp) at given RZ N_parallel points
c      and stopped the code
c-------------------------------------------------------------
CWRITE       write(*,*)'in genray.f point_plot_disp=',point_plot_disp
         endif ! iplot_b

      If(point_plot_disp.eq.'rz_nparallel') then        
         call plot_disp_D_at_RZ_n_parallel_points     
c---------close PGgplot 
         call plotend
         stop'after plot_dispd_D_at_RZ_n_parallel_points'
      endif
c--------------------------------------------------------------
c     calculate the OX optimal directions of the EC cone central ray
CWRITE       write(*,*)'genray.f before  gr_OX_optimal_direction'
      if (i_ox.eq.1) then
cSAP080422
         ifreq_write=1 !it is used in dinit_1ray
         call gr_OX_optimal_direction(ndim)
      endif
c---------------------------------------------------------------
CWRITE       write(*,*)'genray.f before   adj_chi_function'
CWRITE       write(*,*)'genray.f i_adj=',i_adj
      if (i_adj.eq.1) then  
         if (i_calculate_or_read_adj_function.eq.1) then
c----------calculate adj chi function used for current drive efficiency
           call adj_chi_function
c          write(*,*)'genray.f ua',ua
c          write(*,*)'th0a',th0a
c          stop 'adj_chi_function'
         endif
c---------------------------------------------------------
c-------read chi function from the file iout5
         call read_iout5_chi

c-------calculate CD conductivity
c        write(*,*)'genray.f 1 ua',ua
c        write(*,*)'th0a',th0a 
         if (i_calculate_or_read_adj_function.eq.1) then  
            call DC_electric_field_adj_conductivity
         endif
         if (i_calculate_or_read_adj_function.eq.0) then  
            call read_iout3_chi
         endif
cSAP080205
c        goto 110
c------------------------------------------------------------
c       calculate numerical arrays of derivatives d_chi_du,d_chi_d_theta
c       from chi function
c------------------------------------------------------------ 
CWRITE       write(*,*)'genray.f before interpolation_chi'
         call interpolation_chi(1,n_radial0,u0,theta0,
     +             chi,d_chi_d_u_out,d_chi_d_theta_out) 
CWRITE       write(*,*)'genray.f after interpolation_chi'
c-------------------------------------------------------------
c       calculate spline coefficients for chi 2D function
c       for all radial points
c--------------------------------------------------------------
CWRITE       write(*,*)'genray.f before splcoef_chi'
         call  splcoef_chi

CWRITE       write(*,*)'genray.f after splcoef_chi'

c        stop 'stop after after splcoef_chi'
c----------------------------------------------------
c      creates data (arrays) for plot like Fig. 1
c      at Karney article Nuclear Fusion 1991 p. 1934
c      using ADJ efficiency
c--------------------------------------------------- 

CWRITE       write(*,*)'genray.f before plot_1_ADJ_Karney'
         call plot_1_ADJ_Karney
CWRITE       write(*,*)'genray.f after plot_1_ADJ_Karney'


         goto 101
ctest efficiency fo LH case
         clight=2.99792458d10     !light speed [cm/sec]
         u_ph_karney=2.d0         !LH phase velocity normalized to unorm
         u_ph_karney=4.d0 
c        u_ph_karney=1.d0 
         u1=u_ph_karney/dsqrt(2.d0) 
         thetapol_l=0.d0
c        do n_radial=1,npsi0
         do n_radial=2,2
CWRITE       write(*,*)'genray.f n_radial,psis(n_radial)',
CWRITE      +n_radial, psis(n_radial)
            rho_loc=rhopsi(psis(n_radial))
            temp_kev=temperho(rho_loc,1)
            unorm=dsqrt(temp_kev/t)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m) 
            cnpar=(clight/unorm)/u_ph_karney
            cnpar=dsqrt(1.d0+((clight/unorm)/u_ph_karney)**2)
CWRITE       write(*,*)'genray.f rho_loc,unorm,cnpar',rho_loc,unorm,cnpar
CWRITE       write(*,*)'zeffrho(rho_loc)',zeffrho(rho_loc)
            efficien=8.d0/(5.d0+zeffrho(rho_loc))*u1*u1 !test
CWRITE       write(*,*)'asimptotic efficien',efficien
            call CD_adj_LH_efficiency_anal_1(u_ph_karney,cnpar,
     +           n_radial,thetapol_l,
     +           unorm,
     +           lh_cd_efficiency)
CWRITE       write(*,*)'unorm',unorm
CWRITE       write(*,*)'lh_cd_efficiency',lh_cd_efficiency
         enddo
 101     continue
      endif ! i_adj.eq.1
c----------------------------------------------------------------

      endif  !On myrank=0    !----------------------------myrank=0

      do iray=1,nray         
CWRITE       write(*,*)'4 iray,arzu0(iray),arru0(iray),arphiu0(iray)',
CWRITE      +iray,arzu0(iray),arru0(iray),arphiu0(iray)
      enddo

      i_total_bad_initial_conditions=0
      power_launched=0.d0

CMPIINSERTPOSITION BARRIER
      if(myrank.eq.0) then  ! MPI
         call cpu_time(time_before_1st_ray)
      endif  !On myrank=0 MPI


      do 21 iray=1,nray    ! <<<<==== MAIN LOOP
c$$$         if (iray.eq.17) then
c$$$            write(*,*)'SKIPPING A RAY: Normally, remove this'
c$$$            goto 21
c$$$         endif
c$$$         if (iray.ne.17) then
c$$$            write(*,*)'SKIPPING ALL BUT 1 RAY: Normally, remove this'
c$$$            goto 21
c$$$         endif
c--------the loop over all rays        

      if(myrank.eq.0) then  ! MPI
         write(i1_,11)iray
 11      format('####  iray=',i3)
         call cpu_time(time_before_rk)
      endif  !On myrank=0 MPI

ctest for hamiltonian mapping
c         write(*,*)'genray before maphmnri'
c         ihermloc=1
c         ye=0.900751703
c         xe=2.52756889
c         write(*,*)'r0x',r0x
c         r=0.1662d+01/r0x
c         z=0.4916d-05/r0x
c         phi=0.8992d-01
c         bmod=b(z,r,phi)
c         write(*,*)'genray rho=',rho
c         te=tempe(z,r,phi,1)
c         write(*,*)'te=',te
c         stop
c         cnpar=0.762634654
c         cnper=25.3389475
c         delnper=0.1d0
c         delnperim=0.002d0
c         call maphmnri(ye,xe,te,cnpar,cnper,delnper,delnperim,ihermloc)
c         write(*,*)'genray after  maphmnri'
c         stop
cendtest

      if (i_emission.eq.0) then
      nfreq=1
CWRITE       write(*,*)'in genray i_emission=0, no emission calculation'
CWRITE       write(*,*)'genray.f put nfreq=1' 
      endif


      i_geom_optic_loc=i_geom_optic
      do 25 ifreq=1,nfreq

      i_geom_optic=i_geom_optic_loc

c for Nazikian case !070802
c           if (ifreq.gt.2) id=2

CMPIINSERTPOSITION DETERMINERANK

      write(*,*)
      write(*,*) '  iray,ifreq,nray,myrank ============= ', 
     +              iray,ifreq,nray,myrank
       
      ifreq_write=ifreq ! to set in write.i

CWRITE       write(*,*)'genray.f ifreq_write',ifreq_write

      call cpu_time(time_loop_ifreq_1)

c-----------the loop over all emission frequencies for i_emission=1
CWRITE       write(*,*)'### iray= ',iray,'ifreq= ',ifreq
c-----------set the frequency frqncy=wfreq(ifreq) and w0,v0,w0,v(i),w(i)
      call set_freq(ifreq)

      if(istart.eq.1) then
c ------------EC wave
        zst1=zstj(iray) 
        rst1=rstj(iray) 
        phist1=phistj(iray) 
        alfast1=alphaj(iray)
        betast1=betaj(iray)
CWRITE       write(*,*)'genray before dinit_1ray ,zst1,rst1,phist1',
CWRITE      +zst1,rst1,phist1,'alfast1,betast1',alfast1,betast1
CWRITE       write(*,*)'in genray iray,powj(iray)',iray,powj(iray)
        powini=powj(iray)
CWRITE       write(*,*)'genray.f before dinit_1ray ioxm ',ioxm  
        call dinit_1ray(zst1,rst1,phist1,alfast1,betast1,
     +  cnteta,cnphi,u,iraystop)
      endif !istart.eq.1

      if((istart.eq.2).or.(istart.eq.3))then
c ------------LH and FW waves >istart=2
c ------------ECR O_X mode conversion special case: istart=31 
CWRITE       write(*,*)'2 in genray before dinit_1ray iray',iray
c      write(*,*)'.......................................'
c       write(*,*)'2 in genray before dinit_1ray zu0,ru0,phiu0',
c     + rank, arzu0(iray),arru0(iray),arphiu0(iray)
c      write(*,*)'.......................................'
         powini=powinilh(iray)
CWRITE       write(*,*)'alfast1,betast1,arntheta(iray),arnphi(iray)',
CWRITE      +alfast1,betast1,arntheta(iray),arnphi(iray)
CWRITE       write(*,*)'!!!genray before dinit_1ray'
CWRITE       write(*,*)'genray.f before dinit_1ray ioxm ',ioxm
         call dinit_1ray(arzu0(iray),arru0(iray),arphiu0(iray),
     +   alfast1,betast1,arntheta(iray),arnphi(iray),u,iraystop)
      endif !istart.eq.2 .or. istart.eq.3

      i_bad_initial_conditions=0
CWRITE       write(*,*)'genray.f after dinit_1ray ',
CWRITE      +'i_bad_initial_conditions=',i_bad_initial_conditions
cBH140909: Allow for starting with (near-) zero power
cBH140909:      write(*,*)'genray.f: iraystop,powini,iray_status_one_ray=',
cBH140909:     1     iraystop,powini,iray_status_one_ray
      if (iraystop.eq.1) then
         if (iray_status_one_ray.ne.14) then  !i.e., not zero pwr start
            i_bad_initial_conditions=1
            nrayelt=0
         else
            i_bad_initial_conditions=0
            nrayelt=1 ! "short" ray: near zero initial power
         endif
         ! Ray tracing will be skipped. But send/receive data (mpi)
         ! This is needed for proper logic/control of mpi data flow.    
      else
         !irayl=irayl+1 !YuP-130531: counting 'good' rays (for data saving)
      end if

c-----------initialization for each new ray 
      prmt(6)=prmt6 
      prmt(7)=prmt(1)+prmt(6)
CWRITE       write(*,*)'genray.f after goto 24 prmt(1),prmt(6),prmt(7)',
CWRITE      +prmt(1),prmt(6),prmt(7)
 
c----------------------------------------------------------
c           call b() to calculate the small radius rho (inside b())
c           the result will be in common block  one.i
c           rhooldsc will be used in the n_perp_scattering procedure
c           rhooldsc will be in common block scatnper.i 
c-----------------------------------------------------------
      bmod=b(u(1),u(2),u(3))
      rhooldsc=rho

      if(qsize.gt.1) then
         if( myrank.eq.0 ) goto 23 
         !rank=0: Skip ray-tracing; just wait for data from other ranks
      endif

      if(iraystop.eq.1) goto 22 ! rank>0: no ray or "short" ray
         !(either bad init.condition or near-zero init.power)
         !Skip tracing, but send data (initial values)

      nstep_rk=1 ! initialize the number of Runge-Kutta time step 
            
      if (isolv.eq.1) then
c--------------------------------------------------------------
c              The Runge-Kutta solution of 6 hamiltonian equations
c              with correction which gives 
c              hamiltonian conservation with
c              accuracy epscor=prmt(4)
c--------------------------------------------------------------
      if(irkmeth.eq.0) then
c                4_th order Runge-Kutta method with constant time step
         call drkgs(prmt,u,deru,ndim,ihlf,rside1,outpt,aux)
      endif
              
      if(irkmeth.eq.1) then
c                 5-th order Runge-Kutta method with variable time step,
         call drkgs1(prmt,u,deru,ndim,ihlf,rside1,outpt,aux)
      endif
      if(irkmeth.eq.2) then
c                 4th order Runge-Kutta with variable time step,
c                          time step can be reduce or enlarge                
CWRITE       write(*,*)'genray.f before  drkgs2 ioxm.ndim',
CWRITE      +ioxm,ndim
         if(myrank.eq.0)then
            call cpu_time(time_drkgs2_1)
         endif
         !-------------------------------------------------
         call drkgs2(prmt,u,deru,ndim,ihlf,rside1,outpt,aux,
     +   i_output)
         !-------------------------------------------------
         if(myrank.eq.0)then
            call cpu_time(time_drkgs2_2)
            WRITE(*,*)'iray, time_drkgs2:', 
     +       iray,time_drkgs2_2-time_drkgs2_1
         endif
      endif ! irkmeth=2
      endif ! isolv.eq.1
          
      if (isolv.eq.2) then
c-------------------------------------------------------------------
c              The Runge-Kutta solution of 5 hamiltonian equations
c              and the solution of the dispersion relation for the sixth
c              ray variable'
c-------------------------------------------------------------------------
         call rkb1(prmt,u,deru,ndim,ihlf,rsideb1,outptb1,aux)
c------------------------------------------------------------------------
      endif ! isov.eq.2
                          
 10   format(a/a)


c         do i=1,NR-1,2         
c             write(*,*)'After RK4 i,power(i)',i,power(i)
c         enddo

 22   continue ! rank>0: skip tracing for a "short" ray (iraystop=1)

CMPIINSERTPOSITION SENDDATA
cyup At this point, other ranks are sending data to rank=0.
cyup and then skipping all the rest of the loop
cyup (SENDDATA contains goto 25 line, meaning "advance to next ray").
cyup So, skipping emission calc. and skipping write_data to nc file.
cyup But such rank>0 would not proceed to the next ray until 
cyup it gets confirmation from rank=0 that it got the data 
cyup for a given ray and recorded all data (subr.write3d).
cyup See SENDCONFIRM line just before label 25 continue:
cyup this is the "confirmation" sent by rank=0 

 23   continue ! skip handle (rank=0 skips ray tracing)

CMPIINSERTPOSITION RECVDATA
cyup Receiving data from ranks>0.
cyup But not sending "confirmation" yet to a given particular rank.
cyup Still need to do other job: write3d - record data to nc file.
cyup The subr. write3d uses data in write.i common blocks, 
cyup which will be over-written by each new ray.
cyup So, cannot receive data from any new rank until
cyup data is recorded and "confirmation" is sent (just before label 25).
cyup From this line upto label 25 - it is for rank=0 only.
      
c      write(*,*)'sums after recvdata:',
c     +     nrayelt+nrayelt_emis+ifreq0+i_ox_conversion+
c     +     ifreq_write+iray_status_one_ray,  sum(w_ires),
c     +     sum(cwexde),sum(cweyde),sum(cwezde),sum(w_ceps),
c     +     sum(delpow_e_ar),sum(delcur_par_ar)
c      ! pause
           
c------------------------------------------------------------------
c           creation of the file: genray.bin for xdraw
c           input data from common blocks gr.cb and write
c           write(*,*)'in genray before call mk_graph iray',iray
      if(outxdraw.eq.'enabled')then ! YuP[2018-01-17] Added
      if (i_emission.eq.0) then
         call mk_graph(iray,nray,ifreq)       
         call mk_gr3d(iray,nray)
c              creation of the file: npar.bin for xdraw
         call GRAPHnpr
         if (iwcntr.eq.1) then
c                 calculation of contours wb_c = consts
           call mk_grapc(iray,iwopen,iwj)
         endif
      endif ! i_emission=0
      
      if ((i_emission.eq.1)  .and.  (nrayelt.ge.2)) then 
                             ! skipping short rays
CWRITE       write(*,*)'genray.f before call mk_graph'
CWRITE       write(*,*)'genray. f (i_emission.eq.1) before mk_graph' 
         call mk_graph(iray,nray,ifreq)
CWRITE       write(*,*)'genray.f after call mk_graph'
      endif 
      endif ! outxdraw
c--------------------------------------------------------------------
c           calculation of power(array spower(NR)) and
c           current(array scurrent(NR)) radial profiles
c           as sum of profiles Sum(i=1,iray)

c---------------------------------------------------------------------
CWRITE       write(*,*)'ionetwo=',ionetwo
      if((ionetwo.eq.1)   .and.  (nrayelt.ge.1)) then
         call sonetwo
      endif
c------------------------------------------------------------
      if (i_emission.eq.1   .and.  (nrayelt.ge.2)) then
                            ! skipping short rays
         call cpu_time(time_emission_1)                
c--------------spline coeficients calculations for ray data along the ray
         call spl_ray

c--------------emission calclulation   
CWRITE       write(*,*)'genray.f before emission nrayelt_o_cutoff',
CWRITE      +nrayelt_o_cutoff

c               do n=1,nrayelt
c               write(*,*)'n,wsn(n)',n,wsn(n)
c               enddo

         call emission(tol_emis,nrayelt,nrelta,wsn,wal_emis,
     +   wj_emis,wnray,win_sn,win_0,
     +   wi_0(iray,ifreq),wi_0sn,wtaun_em,
     +   nrayelt_emis,
     +   nrayelt_o_cutoff,transm_ox,nrayelt_o_cutoff_emis)  

c               write(*,*)'genray.f aft emission ifreq,nrayelt_o_cutoff',
c     &         ifreq,nrayelt_o_cutoff

c               write(*,*)'genray.f aft emis ifreq,iray,wi_0(iray,ifreq)'
c     &         ,ifreq,iray,wi_0(iray,ifreq)

c               write(*,*)'genray.f after emission  nrayelt_emis',
c     &         nrayelt_emis

c              do n=1,nrayelt_emis
c                write(*,*)'n,wi_0sn(n)',n,wi_0sn(n)
c              enddo 
               
c               write(*,*)'in genray.f nrayelt,nrayelt_emis,
c     +         wi_0(iray,ifreq)',
c     +         nrayelt,nrayelt_emis,wi_0(iray,ifreq)

               

c              do n=1,nrayelt_emis
c                write(*,*)'n,wr_em(n),win_sn(n),win_0(n),wtaun_em(n)',
c     +          n,wr_em(n),win_sn(n),win_0(n),wtaun_em(n)
c              enddo   

c              wtau_em is the optical depth: 
c              a) from one pass, if the ray had the reflection point 
c              b) from part of the ray path, if there was full absorption
c                 before one reflection
         wtau_em(iray,ifreq)=wtaun_em(nrayelt_emis-1)

cSm070120      wi_0 was devided by wnray(1) in emission.f
c               wi_0(iray,ifreq)=wi_0(iray,ifreq)/wnray(1)**2
CWRITE       write(*,*)'genray.f wnray(1) ',wnray(1)

         wi_0t(iray,ifreq)=wi_0(iray,ifreq)/
     +   (1.d0-wallr*dexp(-wtau_em(iray,ifreq)))

c--------------write the emission output data to emis.bin file for plotting
CWRITE       write(*,*)'genray.f before mk_gremis'
c               write(*,*)'wsn',wsn
c               write(*,*)'win_0',win_0
c               write(*,*)'genray.f bef mk_gremis iray,nray,ifreq,nfreq'
c     &                                          ,iray,nray,ifreq,nfreq
         call mk_gremis (iray,nray,ifreq,nfreq)
c--------------------------------------------------------------------
c              Find the ray point P0 with the maximal value of flux
c              (at the plasma boundary from n-th bin) divided by 
c              the bin length ds : In0/ds
c              Calculate the plasma temperatute in this point P0.
c--------------------------------------------------------------------
         n0=1
         win_0dsmax=0.d0

         do n=1,nrayelt_emis-1
            ds=wsn(n+1)-wsn(n) !the bin length       
c                  write(*,*)'n,wsn(n+1),wsn(n),win_0(n),ds',
c     +            n,wsn(n+1),wsn(n),win_0(n),ds
c                  write(*,*)'genray.f n,n0,win_0dsmax,win_0(n)/ds',
c     +            n,n0,win_0dsmax,win_0(n)/ds
c                  if (ds.gt.1.d-18) then 
            if (ds.gt.1.d-12) then
            if (win_0dsmax.lt.win_0(n)/ds) then
               win_0dsmax=win_0(n)/ds
               n0=n                  
            endif
            endif
         enddo

         r0_em=wr_em(n0)
         wr0_em(ifreq)=r0_em  
         z0_em=wz_em(n0)
         wz0_em(ifreq)=z0_em  
         phi=0.d0
         bmod=b(z0_em,r0_em,phi)
         wrho0_em(ifreq)=rho
         wtemp_pl_fr(ifreq)=tempe(z0_em,r0_em,phi,1)!KeV
c               write(*,*)'genray.f ifreq,wr0_em(ifreq),wz0_em(ifreq)'
c     +         ,ifreq,wr0_em(ifreq),wz0_em(ifreq),
c     +         'wtemp_pl_fr(ifreq)',wtemp_pl_fr(ifreq)
c--------------------------------------------------------
c              calculate the major radius at EC second harmonic point
c              at z=0 
c              (to use drawemfr.in file i_r_2nd_harm=1)
c              (to use drawemf1.in file i_r_2nd_harm=0)
c-------------------------------------------------------
         if(i_r_2nd_harm.eq.1) then
            wr_2nd_harm(ifreq)=r_2nd_harm(1)
CWRITE       write(*,*)'genray.f ifreq,wr_2nd_harm(ifreq)',
CWRITE      +ifreq,wr_2nd_harm(ifreq)
         endif               

c-------------------------------------------------------
c              calculate the bulk electron temperuture at 
c              EC second harmonic point z=0
c-------------------------------------------------------
         if (wr_2nd_harm(ifreq).gt.rmax) then
            wtemp_2nd_harm(ifreq)=0.d0
         else
            psi_loc=fpsi(wr_2nd_harm(ifreq),z) 
            rho_loc=rhopsi(psi_loc)
            wtemp_2nd_harm(ifreq)=temperho(rho_loc,1)
         endif

CWRITE       write(*,*)'genray.f after wr_2nd_harm(ifreq)',
CWRITE      +wr_2nd_harm(ifreq)

         call cpu_time(time_emission_2)     
CWRITE       write(*,*)'time_emission',time_emission_2-time_emission_1     
      endif !emission    
      
CWRITE       write(*,*)'genray.f i_bad_initial_conditions',
CWRITE      +i_bad_initial_conditions

      if(i_bad_initial_conditions.eq.0) then
c--------------------------------------------------------------------
c             calculates total power launched along rays with 
c             good initial conditions
c--------------------------------------------------------------------
         power_launched=power_launched+powini
c---------------------------------------------------------------------
c             writing ray data to mnemonic.txt and/or 
c             saving data for mnemonic.nc
c-------------------------------------------------------------------  
         if(outnetcdf.eq.'enabled')then !YuP[2018-01-17]
          write(*,*)'in genray.f before write3d irayl=',irayl
          call write3d  ! save ray data; from myrank=0 only !
          write(*,*)'in genray.f  after write3d irayl=',irayl  
         endif ! outnetcdf       
              
         if (i_emission.eq.1) then
         if (nfreq.gt.1) then
c---------------------------------------------------------------------
c                   calculate the emission temperature KEV for
c                   multi-frequency case
c---------------------------------------------------------------------
         pi=4.d0*datan(1.d0)
         clight=2.9979d10
         wtemp_rad_fr(ifreq)=2.d0*pi*clight**2*
     +   wi_0(iray,ifreq)/(wfreq(ifreq)*1.d9)**2/1.6022d-9

         wtemp_rad_fr_wall(ifreq)=2.d0*pi*clight**2*
     +   wi_0t(iray,ifreq)/(wfreq(ifreq)*1.d9)**2/1.6022d-9
                  
c                   write(*,*)'genray.f iray,ifreq,wfreq(ifreq)',
c     +             iray,ifreq,wfreq(ifreq)
c                   write(*,*)' wi_0(iray,ifreq)',wi_0(iray,ifreq)
CWRITE       write(*,*)'wtemp_rad_fr(ifreq)',wtemp_rad_fr(ifreq)
CWRITE       write(*,*)'wtemp_rad_fr_wall(ifreq)',
CWRITE      +wtemp_rad_fr_wall(ifreq)
CWRITE       write(*,*)'wtemp_pl_fr(ifreq)',wtemp_pl_fr(ifreq)
         endif ! (nfreq.gt.1))
c-----------------------------------------------------------------------
c                put emission data in writencdf.i
c-----------------------------------------------------------------------
CWRITE       write(*,*)'genray.f before put_emission_in_writencdf_i'
         if(outnetcdf.eq.'enabled')then !YuP[2018-01-17] Added
         call put_emission_in_writencdf_i(ifreq,iray) 
         endif
CWRITE       write(*,*)'genray.f after put_emission_in_writencdf_i'
         endif ! (i_emission.eq.1)
c------------------------------------------------------------
         call cpu_time(time_loop_ifreq_2)
         time_loop_ifreq=time_loop_ifreq_2-time_loop_ifreq_1
CWRITE       write(*,*)'time_loop_ifreq',time_loop_ifreq
     
c------------------------------------------------------------------
c             calculate the spectrum of the emission flux
c-----------------------------------------------------------------
CWRITE       write(*,*)'genray.f i_bad_initial_conditions',
CWRITE      +i_bad_initial_conditions

CWRITE       write(*,*)'genray.f before emission_spectrum'
CWRITE       write(*,*)'ifreq,iray,nrayelt_o_cutoff_emis',
CWRITE      +ifreq,iray,nrayelt_o_cutoff_emis
           
         call emission_spectrum(ifreq,iray,
     +   nrayelt_o_cutoff_emis)
CWRITE       write(*,*)'genray.f after emission_spectrum'

      else
c-------------calculate the total number of rays having 
c             i_bad initial conditions=1
         i_total_bad_initial_conditions= 
     +   i_total_bad_initial_conditions+1 
      endif !i_bad_initial_conditions 
c--------------------------------------------------------------
           
      if(i_bad_initial_conditions.eq.1) then
      iray_status_nc(iray,ifreq)=-1
      else 
      iray_status_nc(iray,ifreq)=iray_status_one_ray
      endif  !i_bad_initial_conditions 
      
CMPIINSERTPOSITION SENDCONFIRM
cyup Very important! The data for a given ray (rank) is recorded.
cyup Confirmation is sent to the given rank: 
cyup The rank of source (src) was found in subr. recvdata(src,...).
cyup Now rank=0 is ready to receive data from other rank.

 25   continue ! ifreq   (advance to next ray)
     
          
c           do ifreq=1,nfreq
c            write(*,*)'iray,ifreq',iray,ifreq
c            write(*,*)'wi_0(iray,ifreq)',wi_0(iray,ifreq)
c            write(*,*)'wi_0_nc(iray,ifreq)',wi_0_nc(iray,ifreq)
c           enddo         
c-------------------------------------------------------------
c           if ((i_emission.eq.1).and.(nfreq.gt.1)) then
c--------------calculate the emission temperature KEV for multi-frequency case
c               pi=4.d0*datan(1.d0)
c               clight=2.9979d10
c               do ifreq=1,nfreq                  
c                  wtemp_rad_fr(ifreq)=2.d0*pi*clight**2*
c     +            wi_0(iray,ifreq)/(wfreq(ifreq)*1.d9)**2/1.6022d-9

c                  wtemp_rad_fr_wall(ifreq)=2.d0*pi*clight**2*
c     +            wi_0t(iray,ifreq)/(wfreq(ifreq)*1.d9)**2/1.6022d-9
                  
c                  write(*,*)'genray.f iray,ifreq,wfreq(ifreq)',
c     +            iray,ifreq,wfreq(ifreq)
c                  write(*,*)' wi_0(iray,ifreq)',wi_0(iray,ifreq)
c                  write(*,*)'wtemp_rad_fr(ifreq)',wtemp_rad_fr(ifreq)
c                  write(*,*)'wtemp_rad_fr_wall(ifreq)',
c     &                       wtemp_rad_fr_wall(ifreq)
c                  write(*,*)'wtemp_pl_fr(ifreq)',wtemp_pl_fr(ifreq)
c               enddo
c           endif
      if (i_emission.ne.0) call mk_gremfr(iray,nray)

      !write(*,*)'genray. before 21================ iray,rank=',iray,rank
      if(myrank.eq.0)then
         call cpu_time(time_after_rk)
         time_rk= time_after_rk-time_before_rk
         WRITE(*,*)'iray, time_rk:', iray,time_rk
      endif
 21   continue ! iray=1,nray 

         
      if ((i_emission.eq.1).and.(nfreq.gt.1))then
c-----------------------------------------------------------------------------
c     calculate the averaged temperature over all rays for multi-frequency case
c-----------------------------------------------------------------------------
      do ifreq=1,nfreq
cSAP090710
c           waveraged_temp_rad_fr_nc(nfreqa)=0.d0
      call bcast(waveraged_temp_rad_fr_nc,0.d0,
     +SIZE(waveraged_temp_rad_fr_nc)) 

      sum_emission=0.d0
      sum_emission_wall=0.d0 

      do iray=1,nray
      sum_emission=sum_emission+wi_0_nc(iray,ifreq)
      sum_emission_wall=sum_emission_wall+wi_0_nc(iray,ifreq)  
c             write(*,*)'ifreq,iray,sum_emission,sum_emission_wall',
c     &                  ifreq,iray,sum_emission,sum_emission_wall
  
      enddo

      if (sum_emission.le.1.d-100) then
      waveraged_temp_rad_fr_nc(ifreq)=0.d0
      waveraged_temp_rad_fr_wall_nc(ifreq)=0.d0
      else
      do iray=1,nray
                
      waveraged_temp_rad_fr_nc(ifreq)=
     +waveraged_temp_rad_fr_nc(ifreq)+
     +wtemp_rad_fr_nc(iray,ifreq)*
     +wi_0_nc(iray,ifreq)/sum_emission
              
      waveraged_temp_rad_fr_wall_nc(ifreq)=
     +waveraged_temp_rad_fr_wall_nc(ifreq)+
     +wtemp_rad_fr_wall_nc(iray,ifreq)*
     +wi_0_nc(iray,ifreq)/sum_emission_wall      
               
      enddo  
      endif
    
      enddo

      endif ! (i_emission.eq.1).and.(nfreq.gt.1)
   
      if( myrank.ne.0 ) then 
         goto 100
      endif ! myrank.ne.0  MPI
      

cSAP100111
c-----------------------------------------------------------------        
c     power iterations for LH_lsc approach
c-------------------------------------------------------------------
CWRITE       write(*,*)'genray.f i_lsc_approach',i_lsc_approach

      if (i_lsc_approach.eq.1)then
CWRITE       write(*,*)'genray.f before  LSC_power_iterations'
         call  LSC_power_iterations
CWRITE       write(*,*)'genray.f after LSC_power_iterations ionetwo',ionetwo
c         if (ionetwo.eq.1) then
c--------------------------------------------------------------------
c           calculation power and current density profiles at radius rho
c---------------------------------------------------------------------       
c           write(*,*)'genray.f after LSC_power_iterations bef dnonetwo'
c           call dnonetwo         
c           write(*,*)'genray.f after LSC_power_iterations aft dnonetwo'

c           write(*,1030)'i rho_bin_center powden currden '
c           do i=1,NR-1
c              write(*,1031)i,rho_bin_center(i),powden(i),currden(i)
c           enddo
c         endif
c         write(*,1030)'i,rho_bin_center_lsc(i) '
c     &      //'power_dens_watt_m3_TSC_1D(i)'
c     &      //'CD_dens_no_Edc_a_m2_TSC_1D(i)'

c         do i=1,n_psi_TSC
c            write(*,1031)i,rho_bin_center_lsc(i),
c     &      power_dens_watt_m3_TSC_1D(i),
c     &      CD_dens_no_Edc_a_m2_TSC_1D(i)

c         enddo

c1030     format(/,1x,a)
c1031     format(i3,3(1pe12.4))

c--------plot to the plot.ps file power and CD density radial profiles 
c         call calc_RF_power_CD_density_profile
c         call plot_lsc_powdens_rho
c--------create  lsc.bin file to plot delpwr(ws)
         call MK_GRAPH_LSC(nray)
c--------plot lg(fe(v_par) at shosen radii
         if (n_dim_radii_f_plots.ge.1) then
            do i=1,n_dim_radii_f_plots 
c               call plot_log_fe(n_radii_f_plots_ar(i)) 
c               call plot_log_fe_param(n_radii_f_plots_ar(i)) 
               call plot_log_fe_vpar2_param(n_radii_f_plots_ar(i)) 
            enddo   
         endif
      endif ! i_lsc_approach.eq.1
c-------end lsc ----------------------------------------

      if(outnetcdf.eq.'enabled')then !YuP[2018-01-17] Added ------------
      
        if(rayop.eq."text" .or. rayop.eq."both") then
c--------close file for 3d FP code
          close(i_)
        endif

CWRITE       write(*,*)'genray: powtot_e,powtot_i',powtot_e,powtot_i
        if(rayop.eq."netcdf" .or. rayop.eq."both") then
c         write(*,*)'read before wrtnetcdf(1)'
c          read(*,*)
c 110  continue

        call wrtnetcdf(1)
CWRITE       write(*,*)'genray.f after  wrtnetcdf(1)'
c--------it will write data to nc file in one ray point
c        at number of point is=80 and
c        at number of ray iray=1

CWRITE       write(*,*)'before wrtnetcdf`_one_ray_point(1,80,1)'

        call wrtnetcdf_one_ray_point(1,80,1)
    
        call wrtnetcdf(0)              
CWRITE       write(*,*)'genray.f after  wrtnetcdf(0)' 

c--------it will write data to nc file in one ray point
c        at number of point is=80 and
c        at number of ray iray=1
CWRITE       write(*,*)'before wrtnetcdf_one_ray_point(0,80,1)'

        call wrtnetcdf_one_ray_point(0,80,1)

CWRITE       write(*,*)'genray.f afterwrtnetcdf_one_ray_point(0,80,1)'
c-------------------------------------------------------------
c        writesw peqdsk and r,z eqdsk mesh to existing 
c           netcdf file
c-------------------------------------------------------------
        call netcdf_eqdsk_data(trim(filenc))
c-------------------------------------------------------------
        if ((n_wall.ne.0).or.(max_limiters.ne.0)) then
c-----------writes wall and limiter coordinates to existing 
c           netcdf file
CWRITE       write(*,*)'genray.f before netcdf_wall_lim_data filenc',
CWRITE      +filenc 
          call wrtnetcdf_wall_limiter_data(trim(filenc))
        endif
      
        endif ! rayop=
CWRITE       write(*,*)'genray: powtot_e,powtot_i',powtot_e,powtot_i
  
        if ((istart.eq.2).or.(istart.eq.3)) then
c--------write ray starting coordinates in filenc.nc file
CWRITE       write(*,*)'genray.f before wrtnetcdf_grill_launch filenc',
CWRITE      +filenc
          call wrtnetcdf_grill_launch(trim(filenc))
        endif
        if (istart.eq.1) then
c--------write ray starting coordinates in filenc.nc file
CWRITE       write(*,*)'genray.f before wrtnetcdf_EC_launch filenc',
CWRITE      +filenc
c         write(*,*)'read before wrtnetcdf_EC_launch'
c         read(*,*)
          call wrtnetcdf_EC_launch(trim(filenc))
c         write(*,*)'read after wrtnetcdf_EC_launch'
c         read(*,*)
        endif

        if(myrank.eq.0)then ! MPI write from rank=0 only !
           close(i1_)
        endif ! myrank=0 MPI
      
      endif ! outnetcdf ------------------------------------------------

cSAP091030
      if(ionetwo.eq.1) then
c--------------------------------------------------------------------
c        calculation power and current density profiles at radius rho
c---------------------------------------------------------------------
CWRITE       write(*,*)'genray.f before dnonetwo ionetwo',ionetwo
         if(myrank.eq.0) call dnonetwo         
CWRITE       write(*,*)'genray.f after dnonetwo ionetwo',ionetwo
      endif
 
cSAP090306
c 1010 format('total power absorbed at reflections(erg/sec)=',1pe14.6)
c      write(*,1010) w_tot_pow_absorb_at_refl
 1010 format('total power absorbed at reflections(watt)=',1pe14.6)
CWRITE       write(*,1010) w_tot_pow_absorb_at_refl*1.d-7

c--------------------------------------------------------------------
c     print out status of rays ending
c--------------------------------------------------------------------
         WRITE(*,*)'status of rays ending'
         do iray=1,nrayl
         do ifreq=1,nfreq
            WRITE(*,*)'iray,ifreq,iray_status_nc(iray,ifreq)',
     +                 iray,ifreq,iray_status_nc(iray,ifreq)
         enddo
         enddo
c-----------------------------------------------------------
cSAP100111
c-----------------------------------------------------------------        
c     power iterations for LH_lsc approach
c-------------------------------------------------------------------
CWRITE       write(*,*)'genray.f i_lsc_approach',i_lsc_approach

      if (i_lsc_approach.eq.1)then
c        write(*,*)'genray.f before  LSC_power_iterations'
c         call  LSC_power_iterations
c         write(*,*)'genray.f after LSC_power_iterations ionetwo',ionetwo
         if (ionetwo.eq.1) then
c--------------------------------------------------------------------
c           calculation power and current density profiles at radius rho
c---------------------------------------------------------------------       
c           write(*,*)'genray.f after LSC_power_iterations bef dnonetwo'
            if(myrank.eq.0) call dnonetwo         
c           write(*,*)'genray.f after LSC_power_iterations aft dnonetwo'
c           write(*,1030)'i rho_bin_center powden currden '
c           do i=1,NR-1
c              write(*,1031)i,rho_bin_center(i),powden(i),currden(i)
c           enddo
         endif
c         write(*,1030)'i,rho_bin_center_lsc(i) '
c     &      //'power_dens_watt_m3_TSC_1D(i)'
c     &      //'CD_dens_no_Edc_a_m2_TSC_1D(i)'
c         do i=1,n_psi_TSC
c            write(*,1031)i,rho_bin_center_lsc(i),
c     &      power_dens_watt_m3_TSC_1D(i),
c     &      CD_dens_no_Edc_a_m2_TSC_1D(i)
c         enddo
c1030     format(/,1x,a)
c1031     format(i3,3(1pe12.4))
c--------plot to the plot.ps file power and CD density radial profiles 
         call calc_RF_power_CD_density_profile
         call plot_lsc_powdens_rho
      endif ! i_lsc_approach.eq.1


CWRITE       write(*,*)'genray.f, partner=',partner
      if(outnetcdf.eq.'enabled')then !YuP[2018-01-17] Added
      if (partner.eq.'genray_profs_in.txt' 
     +.or. partner.eq.'genray_profs_in.nc' 
     +.or. partner.eq.'genray_profs_out.nc') then
c-----------------------------------------------------------------
c       write power density and current density profiles
c       at ONETWO small radial grid in the file: genray_profs_out
c-----------------------------------------------------------------
CWRITE       write(*,*)'genray.f before write_transport_prof'
         igenray=19         
c         call write_transport_prof(NR,nbulk,igenray,indexrho,
c     &   powden_e,powden_s,powden_i,currden,powtot_e,powtot_i,
c     &   powtot_s,powtott,currtot) 
c-------- Here: s_cur_den_onetwo is ONETWO RF current density <j_par*B>/B_0
c               currtot is the total toroidal RF current
         if(myrank.eq.0)then ! MPI write from rank=0 only !             
         call write_transport_prof(NR,nbulk,igenray,indexrho,
     +   powden_e,powden_s,powden_i,s_cur_den_onetwo,powtot_e,powtot_i,
     +   powtot_s,powtott,currtot) 
         endif ! myrank=0 MPI
      endif ! partner=
      endif ! outnetcdf

c--------------------------------------------------------------------
c     Output power and current density profiles at radius rho
c--------------------------------------------------------------------
      if(outxdraw.eq.'enabled')then ! YuP[2018-01-17] Added
      if(ionetwo.eq.1) then
CWRITE       write(*,*)'genray.f before mk_gronetwo ionetwo',ionetwo
         call mk_gronetwo
         call mk_gronetwo_1
      endif
      endif ! outxdraw

      if(outnetcdf.eq.'enabled')then !YuP[2018-01-17] Added
      if(myrank.eq.0)then
CWRITE       write(*,*)'genray.f before wrtnetcdf_prof 1'
         call wrtnetcdf_prof(trim(filenc),1)
CWRITE       write(*,*)'genray.f before wrtnetcdf_prof 0'
         call wrtnetcdf_prof(trim(filenc),0)
CWRITE       write(*,*)'genray.f after wrtnetcdf_prof 0'
      endif ! myrank=0
      endif ! outnetcdf

c--------------------------------------------------------------------
c     creates some data for drawing
c--------------------------------------------------------------------

CWRITE       write(*,*)'genray before mkgrtool' 

      if(itools.eq.1) call mkgrtool

CWRITE       write(*,*)'genray after mkgrtool' 
c--------------------------------------------------------------------
c     write dielectric tensor to mnemonic.nc file
c--------------------------------------------------------------------
      if(outnetcdf.eq.'enabled')then !YuP[2018-01-17] Added
      if (dielectric_op.eq.'enabled') then
c         write(*,*)'read before wrtnetcdf_eps'
c         read(*,*)
         call  wrtnetcdf_eps(trim(filenc))
c         write(*,*)'read after wrtnetcdf_eps'
c         read(*,*)
      endif
      endif ! outnetcdf

c--------------------------------------------------------------------
c     write emission data to mnemonic.nc file
c--------------------------------------------------------------------
      if (i_emission.eq.1) then
CWRITE       write(*,*)'genray.f before wrtnetcdf_emission filenc=',
CWRITE      +filenc
c        write(*,*)'genray.f 3'
c        do iray=1,nray
c         do ifreq=1,nfreq
c            write(*,*)'iray,ifreq',iray,ifreq
c            write(*,*)'wi_0(iray,ifreq)',wi_0(iray,ifreq)
c            write(*,*)'wi_0_nc(iray,ifreq)',wi_0_nc(iray,ifreq)
c         enddo      
c        enddo
c        write(*,*)'read before wrtnetcdf_emission'
c        read(*,*)
         if(outnetcdf.eq.'enabled')then !YuP[2018-01-17] Added
           call wrtnetcdf_emission(trim(filenc),nray)
CWRITE       write(*,*)'genray.f after wrtnetcdf_emission filenc=',
CWRITE      +filenc
           call wrtnetcdf_emission_spectrum(trim(filenc))
CWRITE       write(*,*)'genray.f after wrtnetcdf_emission_spectrum filenc='
CWRITE      +,filenc
         endif ! outnetcdf
         call read_nc(trim(filenc))
CWRITE       write(*,*)'genray.f after read_nc'
      endif
c--------------------------------------------------------------------
c     write iray_status_nc into netcdf file: mnemonic.nc file
c--------------------------------------------------------------------
c      do iray=1,nray  
c         do ifreq=1,nfreq
c            write(*,*)'genray.f iray.ifreq,iray_status_nc(iray,ifreq)',
c     &                          iray,ifreq,iray_status_nc(iray,ifreq)
c          enddo
c      enddo
      if(outnetcdf.eq.'enabled')then !YuP[2018-01-17] Added
      call writencdf_iray_status(trim(filenc))
      endif
c-------------------------------------------------------------------
CWRITE       write(*,*)'nrayelt_o_cutoff', nrayelt_o_cutoff
c                             
c  
      call cpu_time(time_genray_2)
      
      WRITE(*,'(a,1pd15.6)') 'CPU time from t_just_before_equilib=',
     +              time_genray_2-time_genray_1a
     
      WRITE(*,'(a,1pd15.6)') 'CPU time from t_just_before_1st_ray=',
     +              time_genray_2-time_before_1st_ray
     
      WRITE(*,1003) time_genray_2-time_genray_1
 1003 format('CPU TOTAL runtime [sec] ',1pd15.6)
 
      WRITE(*,1004)
 1004 format('genray.f: Normal end of program')
c-----close PGgplot 
      call plotend

 100  continue ! handle to skip plotting and writing, for myrank>0

CMPIINSERTPOSITION ENDTIME
 
CMPIINSERTPOSITION ENDMPI
     
      !!!WRITE(*,*)'AFTER MPI_FINALIZE: myrank,ierr=',myrank,ierr

      ! stop ! stop here gives error messages from MPI

      end program



c=====================================================================
c=====================================================================
CMPIINSERTPOSITION SUBS
c=====================================================================
c=====================================================================



c *****************************************************************
c ********* Tokamak data output to con1 ***************************
c *****************************************************************
      subroutine output_con1
      implicit double precision (a-h,o-z)
      include 'param.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'three.i'
      
      if(myrank.ne.0) return
      
      open(30,file='con1')
      write(30,90) nxeqd,nyeqd
90    format(8x,' the data of input file',/,
     +2x,'the numbers of points in the r direction:',i3,/,
     +2x,'the numbers of points in the z direction:',i3)
      write(30,11) xdimeqd*r0x,ydimeqd*r0x,reqd*r0x,redeqd*r0x,
     +ymideqd*r0x
11    format(2x,'the full-width of rectangle: dx=',f12.6,' m',/,
     +30x,' dy=',f12.6,' m',/,
     +2x,'the major radius of the torus:',f12.6,' m',/,
     +2x,'the major radius of the inner edge',/,
     +2x,'of rectangular grid:',f12.6,' m',/,
     +2x,'the vertical shift up-down symmetry plane:',f12.6,' m')
      write(30,12) xma*r0x,yma*r0x,psimag*b0*r0x**2,
     +psilim*b0*r0x**2,beqd*b0
12    format(2x,'the major radius of magnetic axis:',f12.6,' m',/,
     +2x,'the vertical height of magnetic axis:',f12.6,' m',/,
     +2x,'the poloidal flux function values',/,
     +2x,'at the magnetic axis:',f12.6,/,
     +2x,'and the last closed flux surface:',f12.6,/,
     +2x,'the toroidal magnetic field',/,
     +2x,'at major radius of the torus:',f12.6,' t')
      write(30,13) toteqd
13    format(2x,'the toroidal curent:',e17.8,' a')
      close(30)

      return
      end


      subroutine onetwoini
c--------------------------------------------
c     set zero to arrays power and current
c     for subroutine p_c_prof
c------------------------------------------
CSAP090630
      implicit none
c      implicit double precision (a-h,o-z)
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'onetwo.i'

c-----locals
      integer i,kk

CWRITE       write(*,*)'onetwoini NR,NRA', NR,NRA


cSAP090625
c      do i=1,NRA
      do i=1,NR
      spower(i)=0.0d0
      spower_e(i)=0.0d0
      spower_i(i)=0.0d0
      spower_cl(i)=0.0d0
      scurrent(i)=0.0d0       
      enddo
cSAP090625
c      do i=1,NRA-1
      do i=1,NR-1
      s_cur_den_parallel(i)=0.d0
      enddo
cSAP090625
c      do kk=1,nbulka
c         do i=1,NRA-1
      do kk=1,nbulk
      do i=1,NR-1
      spower_s(i,kk)=0.0d0
      enddo
      enddo


      return
      end
c*************************contourb*************************************
c  It calculates contours coordinates for the open contours R=R(z):   *
c          modb(r,z)=const and y_(r,z,phi,2)=1./n, n=1,2,...          *                              *
c          arrays ry(j,i) zpsi(i)                                *
c          j=1,100 npsi(number of countours n=1,100)                    *
c          i=0,nz (parameter nz is a number of points in z direction) *
c---------------------------------------------------------------------
c  input data are in common one.i,five.i                            *
c---------------------------------------------------------------------*
c  Output: arrays AR,AZ, zpsi,rpsi into common block 'gr              *
c          (file gr.cb)                                               *
c**********************************************************************

      subroutine contourb
      implicit double precision (a-h,o-z)
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'three.i'
      include 'five.i'
      parameter (nz=20)
      dimension       zy(nz),ry(100,nz)
      double precision ias1r,b
      external b
      character jy*2
      character outdaty*20,chj*8

      if(myrank.ne.0) return

      phi=0.d0
      btotl=b(zrmin,rmin,phi)
      bmod=btotl
      yl_i=y(zrmin,rmin,phi,2)
      btotr=b(zrmax,rmax,phi)
      bmod=btotr
      yr_i=y(zrman,rmax,phi,2)
c      write(*,*)'in contourb yl_i yr_i',yl_i,yr_i
      nll=int(1./yl_i)
      nlr=int(1./yl_i)+1
      nrl=int(1./yr_i)
      nrr=int(1./yr_i)+1
c      write(*,*)'in contourb nll,nlr,nrl,nrr',nll,nlr,nrl,nrr
      hz=(zmax-zmin)/dfloat(nz)
      epsy=1d-5 ! accuracy 
      phi=0.d0
      do i=0,nz
      z=zmin+i*hz
      zy(i)=z
      do j=nlr,nrl
c          determination of ry(j,i) where
c          y(z,ry,phi,2)=1/j (using the binary method)
      tl=rmin
      tr=rmax
      c=1.d0/float(j)
      do while ((tr-tl).gt.epsy)
      t=tl+(tr-tl)*0.5d0
      r=t
      bmod=b(z,r,phi)
      y1=y(z,r,phi,2)-c
      rtr=tr
      bmod=b(z,rtr,phi)
      y2=y(z,rtr,phi,2)-c
      if ((y1*y2).gt.0) then
      tr=t
      else
      tl=t
      end if
      end do
c          -----------------------------------------------
c          end of the binary methode
c          -----------------------------------------------
c         write(*,*)'j,i',j,i
c         write(*,*)'tr-tl,t',tr-tl,t
c         bmod=b(z,r,phi)
c         write(*,*)'y1,y(z,r,phi,2),c',y1,y(z,r,phi,2),c
c         write(*,*)'r,z',r,z
      ry(j,i)=r
      enddo !j      
      enddo !i      
c         CALL ASSIGN("assign -F f77 -N ieee u:82",ier)
CWRITE       write(*,*)'contours 1/Y_i=n will be ploted for n=',nlr,'...',nrl 
 10   format(3(1pe11.3))
 11   format(2(1pe11.3))
 12   format(34(1pe11.3))
      do j=nlr,nrl
      call charnumb(j,chj)
      outdaty=chj//'.dat'
c      write(*,*)outdaty
      open(j,file=outdaty)
      do i=0,nz
      bmod=b(zy(i),ry(j,i),phi)
     
c---------------- idx derivativs order 0.ge.idx.le.3---------------
      idx=0
      ipx=ip
      ipx4=ip+4
      rrr=ry(j,i)
      zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,rrr)
      zp=zzrp
      ipx=im
      ipx4=im+4
      zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,rrr)
      zm=zzrm
      if ((zy(i).le.zp).and.(zy(i).ge.zm)) then
c           write(*,*)'in write r,zm,z,zp',rrr,zm,zy(i),zp
c            write(1,10)ry(j,i),zy(i),dfloat(j)
c            write(2,11)ry(j,i),zy(i)
      write(j,11)ry(j,i)*100.d0,zy(i)*100.d0 ! *100 to get cm
      endif
      enddo
      close(j)
      enddo
CWRITE       write(*,*)'end of contourb'
c      close(1)
c      close(2)
c      close(3)
      return
      end
***********charnumb*******************************************
*     The transformation of the integer j to character chj
*     It is assumed that the number of the decimal numbers in j<=8
******************************************************
c     input : j   is integer 
c     ATTENTION: the number of decimal positions in j should be<=8
c     output: chj is character*8 
c-------------------------------------------------------------
      subroutine charnumb(j,chj)
      integer nj
      parameter (nj=8) ! the max number of the decimal numbers in j 
      character chj*8
      character ikch
      dimension ikch(nj)
      ich0=ichar('0')

      k=1
      jk=j/10

      do while(jk.ne.0) 
      k=k+1
      jk=jk/10
      enddo       
c      write(*,*)'k is the number of the decimal nimbers in j',k
      if(k.gt.nj)then
CWRITE       write(*,*)'in charnumb k>nj'
      stop
      endif
      do i=1,nj
      ikch(i)='0'
      enddo

      jk=j
      do i=1,k
      jkn=jk/10
      ik=jk-jkn*10
      jk=jkn
c       ikch(k) is the decimal number in the j in the k position  
c        The position numeration is from the right side
      ikch(i)=char(ik+ich0)
      enddo
      chj=ikch(8)//ikch(7)//ikch(6)//ikch(5)//ikch(4)//ikch(3)//
     +ikch(2)//ikch(1)
c      write(*,*)'in charnumb chj=',chj
      return
      end

c*************************contrb1 *************************************
c  It calculates contours coordinates for the open contours R=R(z):   *
c          modb(r,z)=const and y_(r,z,phi,2)=1./n, n=1,2,...          *                              *
c          arrays ry(j,i) zpsi(i)                                *
c          j=1,100 npsi(number of countours n=1,100)                    *
c          i=0,20 ( a number of points in z direction)                *
c---------------------------------------------------------------------
c  input data are in common one.i,five.i                            *
c---------------------------------------------------------------------*
c* Output: arrays (r,z,) for the given values 1/Y_i=n to               *
c  files: outputy=n.dat                                       *
c**********************************************************************

      subroutine contrb1
      implicit double precision (a-h,o-z)
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'gr.i'
      include 'write.i'
      dimension       zy(20),ry(100,20)
      double precision ias1r
      character chj*8
      character outdaty*20
      
      if(myrank.ne.0) return
      
      iw_j=iwj
c      write(*,*)'in contrb1 iw_j',iw_j
      phi=0.d0

      nnr=40
      nz=40

      hr=xdimeqd/(nnr-1.d0)
      hz=ydimeqd/(nz-1.d0)
      bmin=b(zrmin,rmin,phi)
c      write(*,*)'0 bmin',bmin
     
      open(5,file='wdwci.doc')
      do i=1,nnr+1
      r=4.d0+0.05d0*(dfloat(i)-1.d0)
      do j=1,nz
      z=-2.d0+0.1d0*(dfloat(j)-1.d0)
      bmod=b(z,r,phi)
      dyci=1.d0/y(z,r,phi,iw_j)
      write(5,10)r,z,dyci
      enddo
      enddo
      close(5)
     
      nnr=40
      nz=30
      nz=71
      hr=xdimeqd/(nnr-1.d0)
      hz=ydimeqd/(nz-1.d0)
      hz=14.d0/(nz-1.d0)
c      hr=0.05d0
c      nz=0.1d0
      open(4,file='btot.doc')
      do i=1,nnr
      r=redeqd +hr*(i-1.d0)
c      r=4.d0+hr*(i-1.d0)
      do j=1,nz
      z=ymideqd+hz*(j-1)-ydimeqd*0.5
      z=-7.0d0+hz*(j-1)
c        write(*,*)'in contrb1 i,z,j,r',i,z,j,r
      btot=b(z,r,phi)
c        write(*,*)'in contrb1 i,z,j,r, btot',i,z,j,r,btot
      write(4,10)r,z,btot
      enddo
      enddo
      close(4)

      btotl=b(zrmin,rmin,phi)
      bmod=btotl
      yl_i=y(zrmin,rmin,phi,iw_j)
      btotr=b(zrmax,rmax,phi)
      bmod=btotr
      yr_i=y(zrman,rmax,phi,iw_j)
CWRITE       write(*,*)'in contourb btotl btotr',btotl,btotr
CWRITE       write(*,*)'in contourb yl_i yr_i',yl_i,yr_i
      nll=int(1./yl_i)
      nlr=int(1./yl_i)+1
      nrl=int(1./yr_i)
      nrr=int(1./yr_i)+1
CWRITE       write(*,*)'in contrb1 nll,nlr,nrl,nrr',nll,nlr,nrl,nrr
      nz=20
      hz=(zmax-zmin)/dfloat(nz)
      epsy=1d-5 ! accuracy 
      phi=0.d0
      do i=1,nz
      z=zmin+i*hz
      zy(i)=z
      do j=nlr,nrl
c          determination of ry(j,i) where
c          y(z,ry,phi,iw_j)=1/j (using the binary method)
      tl=rmin
      tr=rmax
      c=1.d0/float(j)
      do while ((tr-tl).gt.epsy)
      t=tl+(tr-tl)*0.5d0
      r=t
      bmod=b(z,r,phi)
      y1=y(z,r,phi,iw_j)-c
      rtr=tr
      bmod=b(z,rtr,phi)
      y2=y(z,rtr,phi,iw_j)-c
      if ((y1*y2).gt.0) then
      tr=t
      else
      tl=t
      end if
      end do
c          -----------------------------------------------
c          end of the binary method
c          -----------------------------------------------

      ry(j,i)=r
      enddo !j      
      enddo !i      
c         CALL ASSIGN("assign -F f77 -N ieee u:84",ier)
CWRITE       write(*,*)'contours 1/Y_i=n will be plotted for n=',nlr,'...',nrl 
 10   format(3(1pe11.3))
 11   format(2(1pe11.3))
 12   format(34(1pe11.3))
      do j=nlr,nrl
      call charnumb(j,chj)
      outdaty=chj//'.dat'
CWRITE       write(*,*)outdaty
      open(1,file=outdaty)
      do i=1,20
      bmod=b(zy(i),ry(j,i),phi)
c        write(*,*)'j,y=1/j,y(zy(i),ry(j,i),phi,iw_j),zy(i),ry(j,i)',
c     6        j,1.d0/dfloat(j),y(zy(i),ry(j,i),phi,iw_j),zy(i),ry(j,i)
     
c---------------- idx derivativs order 0.ge.idx.le.3---------------
      idx=0
      ipx=ip
      ipx4=ip+4
      rrr=ry(j,i)
      zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,rrr)
      zp=zzrp
      ipx=im
      ipx4=im+4
      zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,rrr)
      zm=zzrm
c         write(*,*)'in bound r,zm,z,zp',rrr,zm,zy(i),zp
      if ((zy(i).le.zp).and.(zy(i).ge.zm)) then
c           write(*,*)'in write r,zm,z,zp',rrr,zm,zy(i),zp
c           write(*,*)'ry(j,i)*100,zy(i)*100)'
c           write(*,*)ry(j,i)*100,zy(i)*100

c            write(1,10)ry(j,i),zy(i),dfloat(j)
c            write(2,11)ry(j,i),zy(i)
      write(1,11)ry(j,i)*100.d0,zy(i)*100.d0 ! *100 to get cm
      WRITE(84) REAL(ry(j,i)*100),REAL(zy(i)*100),
     +REAL(XT(3,NP+1)),REAL(YT(3,NP+1)),
     +REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     +REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     +REAL(salphal(nrayelt))
      WRITE(83,5)ry(j,i)*100,zy(i)*100,
     +XT(3,NP+1),YT(3,NP+1),
     +ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     +spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     +salphal(nrayelt)
      endif
      enddo
      WRITE(84)
      WRITE(83,2)
 2    format(/)
      close(1)
      enddo
 5    format(11(1pe10.3))
      return
      end





c*************************contrb2 *************************************
c  It calculates contours coordinastes for contours:                  *
c          modb(r,z)=const and y_(r,z,phi,2)=1./n, n=1,2,...          *
c          for the case then the closed contours exist.                *
c          arrays ry(j,i) thetac(i)                                *
c          j=1,100 npsi(number of countours n=1,100)                    *
c          i=0,nthetac (number of points in thetac direction)                  *
c---------------------------------------------------------------------
c  input data are in common one.i,five.i                            *
c---------------------------------------------------------------------*
c* Output: arrays (r,z,) for the give values 1/Y_i=n to               *
c  files: outputy=n.dat                                       *
c**********************************************************************

      subroutine contrb2
      implicit double precision (a-h,o-z)
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'gr.i'
      include 'write.i'
      parameter (nthetac=20, nthetac1=nthetac+1)
      dimension       thetac(nthetac),zy(100,nthetac1),ry(100,nthetac1)
      double precision ias1r
      character outdaty*20,chj*8
      
      if(myrank.ne.0) return
      
      iw_j=iwj
      phi=0.d0
c------------------------------------------------------
c     creation of the btot.doc file with the coordinates (r,z) and
c     the values of mod(b(r,z))
c     calculation the coordinates (rbmin zbmin) of the point
c     inside the plasma where btot(rbmin,zbmin)=bmin
c------------------------------------------------------
c      write(*,*)'in contrb2 before nz=40'
      nnr=40
      nz=40
      hr=xdimeqd/(nnr-1.d0)
      hz=ydimeqd/(nz-1.d0)
      bmin=b(zrmin,rmin,phi)
c      write(*,*)'0 bmin',bmin
      open(4,file='btot.doc')
      open(5,file='wdwci.doc')
      do i=1,41
      r=4.d0+0.05d0*(dfloat(i)-1.d0)
      do j=1,40
      z=-2.d0+0.1d0*(dfloat(j)-1.d0)
      bmod=b(z,r,phi)
      dyci=1.d0/y(z,r,phi,iw_j)
      write(5,10)r,z,dyci
      enddo
      enddo
      close(5)
      close(4)
c      stop
      do i=1,nnr
      r=redeqd +hr*(i-1.d0)
      do j=1,nz
      z=ymideqd+hz*(j-1)-ydimeqd*0.5
      btot=b(z,r,phi)
cc          bmod=b(z,r,phi)
cc        dyci=1.d0/y(z,r,phi,iw_j)
c        write(*,*)'in contrb2 i,z,j,r,btot,bmin',i,z,j,r,btot      ,bmin
c---------------- idx derivativs order 0.ge.idx.le.3---------------
      idx=0
      ipx=ip
      ipx4=ip+4
      rrr=r
      zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,rrr)
      zp=zzrp
      ipx=im
      ipx4=im+4
      zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,rrr)
      zm=zzrm
      if ((z.le.zp).and.(z.ge.zm)) then
      if(bmin.gt.btot) then
      bmin=btot
      rbmin=r
      zbmin=z
      endif
      endif
      write(4,10)r,z,btot
cc        write(5,10)r,z,dyci
      enddo
      enddo
      close(4)
c      close(5)
CWRITE       write(*,*)'in contrb2 bmin, zbmin,rbmin',bmin, zbmin,rbmin
c      stop
      btotmin=b(zbmin,rbmin,phi)
      bmod=btotmin
      ybmin_i=y(zbmin,rbin,phi,iw_j)
c      write(*,*)'in contrb2 ybmin_i',ybmin_i
c------------------------------------------------------
c     determination of yl_i(on the left side of plasma)
c     determination of yr_i(on the right side of plasma)
      btotl=b(zrmin,rmin,phi)
      bmod=btotl
      yl_i=y(zrmin,rmin,phi,iw_j)
      btotr=b(zrmax,rmax,phi)
      bmod=btotr
      yr_i=y(zrman,rmax,phi,iw_j)
c      write(*,*)'in contrb2 btotl btotr',btotl,btotr
c      write(*,*)'in contrb2 yl_i yr_i',yl_i,yr_i
c---------------------------------------------------------
c     determination of the ysep=1/Y_ci on the the last closed contour  
      if(yl_i.gt.yr_i) then
      ysep=yl_i
      else
      ysep=yr_i
      endif
c      write(*,*)'in contrb2 ysep',ysep
      ny_min=int(1.d0/ysep+2.d0)
      ny_max=int(1.d0/ybmin_i)
c      write(*,*)'in contrb2 ny_min,ny_max',ny_min,ny_max
c---------------------------------------------------------
c  Calculations coordinates of contours 1/Yci=n                              *
c          r(n,thetac)   z(n,thetac)                             
c          arrays ry(j,i) zy(j,i)                                
c          j=ny_min,ny_max(number of contours )            
c          i=1,nthetac+1(number of points in the poloidal angle)  
      pi=4.d0*datan(1.d0)
      hteta=2.d0*pi/dble(nthetac)
      epsy=1d-3 ! accuracy 
c------------------------------------------------------
c     
c     2)we will create the limiter points using the close flux
c      surface psi(r,z)=psilim*psifactr, here
c      psifactr is a parameter (it must be .le.1) to avoide the
c      problems with the nonmonotonic psi function near the separatrix.
c      psifactr is given in genray.in file (It is in common/one/)
c     ------------------------------------
c      psilim=psimag+(psilim-psimag)*psifactr
c     
c-------------------------------------------------------
c        ipsi=1 !  to calculate contours
c        ipsi=0 !  to read contours data from file:psi.bin
c       if (ipsi.eq.0) then
c         open(1,file='psi.bin',form='unformatted',status='old')
c         do i=1,npsi
c        do j=1,nteta1
c             read(1)zpsi(i,j),rpsi(i,j)
c        end do
c         end do
c       close(1)
c       go to 200
c       endif
c------------------------------------------------------------------
      do 101 i=1,nthetac
      theta=hteta*(dble(i)-0.5d0)
      sintet=dsin(theta)
      costet=dcos(theta)
      tini=0.d0
      tm=3.5d0
      htini=tm*0.02d0
      maxiter=10
      htmin=1.d-3
      do 20 j=ny_max,ny_min,-1
CWRITE       write(*,*)'genray.f contrb2 number of contour Y=1/n j=',j      
      n0=j 
c------------------------------------------------
c          binary method for solution of equation 1/Y(r(t),z(t))=n0
      call zrcntr2(tini,n0,costet,sintet,tm,htini,
     +ierr,zt0,rt0,t0,zbmin,rbmin,epsy)
c------------------------------------------------
      tini=t0
      if (ierr.eq.1) then
      zy(j,i)=zt0
      ry(j,i)=rt0
CWRITE       write(*,*)'Yj,Thi,zy(j,i),ry(j,i)',j,i,zy(j,i),ry(j,i)
      else
CWRITE       write(*,*)'zrcontor gave ierr ',ierr
CWRITE       write(*,*)'it is impossible to find the 1/y surface '
CWRITE       write(*,*)' with n',n0,'j=',j
CWRITE       write(*,*)' i=',i,'theta',theta
CWRITE       write(*,*)' with n0-1',n0-1
CWRITE       write(*,*)' it is possible to change 
CWRITE      +ny_min=arpsi(j-1) or reduse the factor
CWRITE      +psifactr in the subroutine equilib'
      stop
      endif
20    continue
101   continue

      do 40 j=ny_max,ny_min,-1
      zy(j,nthetac1)=zy(j,1)
      ry(j,nthetac1)=ry(j,1)
40    continue
c----------------------------------------------------------
c     if ipsi=1 then continue,write file y.bin
      open(3,file='y.bin',form='unformatted')
      do i=1,npsi
      do j=1,nthetac1
      write(3)zy(i,j),ry(i,j)
      end do
      end do
      close(3)
c------------------------------------------------------------
c     if ipsi=1 then continue,ipsi=0 then read file psi.bin
200   continue
c-------------------------------------------------------------
 10   format(3(1pe11.3))
 11   format(2(1pe11.3))
 12   format(34(1pe11.3))
      do j=ny_min,ny_max
      call charnumb(j,chj)
      outdaty=chj//'.dat'
c      write(*,*)outdaty
      open(1,file=outdaty)
c      write(j,*)'r ',j
c      goto 17 !!!!
      do i=1,nthetac1
      bmod=b(zy(j,i),ry(j,i),phi)
CWRITE       write(*,*)'j,y=1/j,y(zy(i),ry(j,i),phi,iw_j),zy(j,i),ry(j,i)',
CWRITE      +j,1.d0/dfloat(j),y(zy(j,i),ry(j,i),phi,iw_j),zy(j,i),ry(j,i)
     
      WRITE(84) REAL(ry(j,i)*100.),REAL(zy(j,i)*100.),
     +REAL(XT(3,NP+1)),REAL(YT(3,NP+1)),
     +REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     +REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     +REAL(salphal(nrayelt))
      WRITE(83,5)ry(j,i)*100,zy(j,i)*100,
     +XT(3,NP+1),YT(3,NP+1),
     +ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     +spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     +salphal(nrayelt)
      enddo !i
      WRITE(84)
      WRITE(83,2)
 2    format(/)
 17   continue
      close(1)
      enddo !j
 5    format(11(1pe10.3))
c      write(*,*)'end of contrb2'
c      stop
      return
      end
c*************************zrcntr2*********************************** *
c  It calculates  contours coordinates of the countor point              *
c   r(n,teta)   z(n,teta) for the given:                         *
c   1/y_ci(z,r)=n0 and poloidal angle teta0(in radians)                *
c   It using the binary methode                                        *
c---------------------------------------------------------------------
c  input data are in common blocks /three/,/five/,/gr.cb/             *
c  tini -initial (left) value of the t( ray parameter )                   *
c  n0   -given value of the 1/y_ci                                *
c  costet0,sintet0 -for the given poloidal angle                   *
c  rbmin major radius of point mod(b)=min (inside the plasma)             *
c  zbmin Z coordinate of point mod(b)=min (inside the plasma)             *
c  tm   -maximal value of t (right)                               *
c  htini   ininial step of the t                                 *
c  epsy-accuracy of equation solution (=max(abs( 1/y_n-1/y_n+1))
c----------------------------------------------------------------------*
c  Output: ,zt0(n0,teta0),rt0(n0,teta0) and parameter t0               *
c           rt0=xma+t0*costet0 ,  zt0=yma+t0*sintet0                   *
c  ierr -index if the solution was obtained =1      else=0                   *
c**********************************************************************
      subroutine zrcntr2(tini,n0,costet0,sintet0,tm,htini,
     +ierr,zt0,rt0,t0,zbmin,rbmin,epsy)
      implicit real*8 (a-h,o-z)
      include 'param.i'
      include 'three.i'
      include 'one.i'
      iw_j=iwj
      ierr=1
      ht=htini
      phi=0.d0
      n0=-n0
c      write(*,*)'in zrcntr2 n0,costet0,sintet0',
c     1 n0,costet0,sintet0
c      write(*,*)'tini,tm,htini,epsy'
c     1 ,tini,tm,htini,epsy
      t1=tini
      rt1=rbmin+t1*costet0
      zt1=zbmin+t1*sintet0
      bmod=b(zt1,rt1,phi)
      rn1=-1.d0/y(zt1,rt1,phi,iw_j)
 10   t2=t1+ht
c      write(*,*)'10 t2,t1,rt1,zt1,rn1',t2,t1,rt1,zt1,rn1
      t2=dmin1(t2,tm)
c      write(*,*)'t2',t2
      if((t2.eq.tm.and.t1.eq.tm).or.
     +(t2.eq.0.d0.and.t1.eq.0.d0)) then
      ierr=0
c       write(*,*)'in zrcntr2 error exit ierr',ierr
c         write(*,*)'t2,t1,rt1,zt1,rt2,zt2',t2,t1,rt1,zt1,rt2,zt2
c       write(*,*)'in zrcntr2 error rn1,rn2,n0',rn1,rn2,n0
c        error exit
      goto30
      endif

      rt2=rbmin+t2*costet0
      zt2=zbmin+t2*sintet0
      bmod=b(zt2,rt2,phi)
      rn2=-1.d0/y(zt2,rt2,phi,iw_j)
      dn=(dfloat(n0)-rn1)*(dfloat(n0)-rn2)
c      write(*,*)'10 t2,rt2,zt2,rn2,dn',t2,rt2,zt2,rn2,dn
      if(dn.le.0.d0)go to 20
      t1=t2
      rn1=rn2
      goto 10
c-----------------------------------------------------------
c     n0 is between rn1 and rn2
c     binary iteration methode
c-----------------------------------------------------------
 20   continue
      tr=t2
      tl=t1
      do while ((tr-tl).gt.epsy)
      t=tl+(tr-tl)*0.5d0
      r=rbmin+t*costet0
      z=zbmin+t*sintet0
      bmod=b(z,r,phi)
      rn1=-1.d0/y(z,r,phi,iw_j)-n0
      rtr=rbmin+tr*costet0
      ztr=zbmin+tr*sintet0
      bmod=b(ztr,rtr,phi)
      rn2=-1.d0/y(ztr,rtr,phi,iw_j)-n0
      if ((rn1*rn2).gt.0) then
      tr=t
      else
      tl=t
      end if
      end do
c     -----------------------------------------------
c          end of the binary methode
c     -----------------------------------------------
      t0=t
      zt0=z
      rt0=r
 30   continue
c      write(*,*)'the end of zrctr2 t0',t0
      return
      end


      subroutine check_param(i_op)
c-----check parameters in param.i
c     i_op=1 check the parameters for eqdsk
c     i_op=2 check the parameters for grill in genray.in

      implicit double precision (a-h,o-z)
      include 'param.i' ! gives the input parameters 
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'commons.i'
   
      if ((i_op.lt.1).or.(i_op.gt.2)) then
         if(myrank.eq.0)then
         WRITE(*,*)'in genray.f,in check_param: wrong value i_op=',i_op
         WRITE(*,*)'It shoud be (i_op.eq.1).or.(i_op.eq.2)'
         WRITE(*,*)'Change i_op in call check_param(i_op)'
         endif
         stop
      endif

      if (nraymax.lt.0) then
         if(myrank.eq.0)then
         WRITE(*,*)'in param.i nraymax <0,it should be >0'
         WRITE(*,*)'Change nraymax in param.i'
         endif
         stop
      endif 

      if (i_op.eq.1) then
c--------check the parameters for eqdsk 
              
      if (nveqd.gt.nxeqd) then
         if(myrank.eq.0) WRITE(*,10110)
         stop
      endif
10110 format("subroutine equilib-nveqd > nxeqd")
cSm030224
      if ((nxeqd.gt.nxeqda).or.(nyeqd.gt.nyeqda)) then
      if(myrank.eq.0) WRITE(*,1000) nxeqd,nxeqda,nyeqd,nyeqda
 1000 format('in equilib.dat in input',/,
     +'the dimensions of eqdsk (in eqilib.dat) nxeqd or nyeqd',/,
     +'are bigger then the parameters nxeqda or nyeqda in param.i'
     +,/,'nxeqd=',I5,'nxeqda=',I5
     +,/,'nyeqd=',I5,'nyeqda=',I5
     +,/,'Please change nxeqda or nyeqda in param.i')
      stop
      endif

      if (nrya.ne.(max0(nxeqda,nyeqda)+4)) then
         if(myrank.eq.0) then
         WRITE(*,*)'in param.i nry.ne.(max0(nxeqda,nyeqda)+4)'
         WRITE(*,*)'nrya=',nrya
         WRITE(*,*)'max0(nxeqda,nyeqda)',max0(nxeqda,nyeqda)
         WRITE(*,*)'Change nrya in param.i'
         endif
         stop
      endif 
               
cSAP091124      
c         if (nzy.ne.(max0(npsi4,nteta1)+4)) then 
      if (nzy.ne.(max0(npsi,nteta1)+4)) then
         if(myrank.eq.0) then
cSAP091124
c            WRITE(*,*)'in param.i nzy.ne.max0(npsi4,nteta1)+4'
         WRITE(*,*)'in param.i nzy.ne.max0(npsi,nteta1)+4'
         WRITE(*,*)'nzy= ',nzy
cSAP091124
c            WRITE(*,*)'max0(npsi4,nteta1)',max0(npsi4,nteta1)
         WRITE(*,*)'max0(npsi,nteta1)=4',max0(npsi,nteta1)+4
         WRITE(*,*)'Change nzy in param.i'
         endif
         stop
      endif

      goto 10

      endif !i_op=1

      if (i_op.eq.2) then
c--------check the parameters for grill in genray.in
      if (ngrill.gt.ngrilla) then
 20   format('equilib in check_param ngrill>ngrilla',/,
     +'it should be ngrilla.ge,ngrill',/,
     +'ngrilla= ',I4,'ngrill= ',I4,/,
     +'Change ngrilla in param.i or grilld in genray.in')
      if(myrank.eq.0) WRITE(*,20)ngrilla,ngrill
      stop
      endif
 
      nmax=0
      do i=1,ngrill
      if (nnkpar(i).gt.nmax) nmax=nnkpar(i)
      enddo

      if (nnkprmax.ne.nmax)then
      if(myrank.eq.0) WRITE(*,30)nnkprmax,nmax
 30   format('genray.f in check_param: it should be',/,
     +'nnkprmax.ge.max{i=1,ngrill}nnkpar(i)',/,
     +'nnkprmax= ',I4,'max{i=1,ngrill}nnkpar(i)= ',I4)
      goto 10             
      endif 
          
      nmax=0
      do i=1,ngrill
      if (nthin(i).gt.nmax) nmax=nthin(i)
      enddo

      if (nthinmax.ne.nmax)then
      if(myrank.eq.0) WRITE(*,40)nthinmax,nmax
 40   format('genray.f in check_param: it should be',/,
     +'nthinmax.ge.max{i=1,ngrill}nthin(i)',/,
     +'nthinmax= ',I4,'max{i=1,ngrill}nthin(i)= ',I4)
      goto 10             
      endif 

      endif !i_op=2

 10   continue

      return
      end
              


              
      subroutine screen_print_all_parameters ! called from rank=0 only
c---------------------------------------------------
c     Write out all code parameters given in param.i
c---------------------------------------------------
      implicit none
      include 'param.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      if(myrank.eq.0) then
      WRITE(*,*)'*********** code parameteres ************' 
c--------------------------------------------------------------
      WRITE(*,*)'number of rays in common/cone/ EC-cone: nraymax = ',
     +nraymax

      WRITE(*,*)'maximum number of EC cones: nconea = ',nconea

      WRITE(*,*)'maximum number of radial bins at the launching disk'//
     +' :  n_mesh_disk_radial_bin_a = ', n_mesh_disk_radial_bin_a

c************************************************************
c     for common/five/
      WRITE(*,*)'nxeqda = ',nxeqda
      WRITE(*,*)'nyeqda = ',nyeqda
      WRITE(*,*)'nlimit = ',nlimit
      WRITE(*,*)'nx4a=nxeqda+4 = ',nx4a
      WRITE(*,*)'ny4a=nyeqda+4 = ',ny4a
      WRITE(*,*)'nrya=ny4a = ',nrya
      WRITE(*,*)'nlim4=nlimit+4 = ',nlim4
      WRITE(*,*)'It should be nrya=max(nxeqda,nyeqda)+4'
      endif

      if(nrya.ne.(max(nxeqda,nyeqda)+4)) then
      if(myrank.eq.0) then
      WRITE(*,*)'nrya.ne.(max(nxeqda,nyeqda)+4))'
      WRITE(*,*)'It should be nrya=max(nxeqda,nyeqda)+4'
      WRITE(*,*)'Please change nrya for nrya=nx4a or nx4a in param.i'
      WRITE(*,*)'giving nrya=max(nxeqda,nyeqda)+4'
      WRITE(*,*)'and recompile the code'
      endif
      stop ! stop in all cores
      endif  

      if(myrank.eq.0) then
c************************************************************
c     for common/fourb/
      WRITE(*,*)'nves = ',nves 
c************************************************************
c     for common gr.i

      WRITE(*,*)'number of points along each contour '//
     +'in the poloidal direction: nteta = ',nteta

      WRITE(*,*)'number of contours: npsi = ',npsi

      WRITE(*,*)'accuracy for the determination of contours'//
     +'poins coordinates(zpsi,rpsi): epspsi = ',epspsi

      WRITE(*,*)'nteta1=nteta+1 = ',nteta1

      WRITE(*,*)'number of contours for the ploting of tokamak'//
     +'cross section: NL = ',NL

      WRITE(*,*)'number of points along each contours in poloidal'//
     +' direction for plotting the tokamak cross section: NP = ',NP

c************************************************************
c     for common/grill/
      
      WRITE(*,*)'for common/grill/'
      WRITE(*,*)'ngrilla = ',ngrilla
      WRITE(*,*)'nnkprmax = ',nnkprmax
      WRITE(*,*)'nraymaxl = ',nraymax
      WRITE(*,*)',nthinmax = ',nthinmax

      WRITE(*,*)'for N_toroidal N_poloidal initial condition'//
     +' at i_n_poloidal=4' 
      WRITE(*,*)'max{i=1,ngril1}nnktor(i): nnktormax = ',nnktormax
      WRITE(*,*)'max{i=1,ngril1}nnkpol(i): nnkpolmax = ',nnkpolmax 
     
c************************************************************
c     for common/ions/
      WRITE(*,*)'max(nbulk): nbulka = ',nbulka
      WRITE(*,*)'nbulkma=nbulka-1   = ',nbulkma

c************************************************************
c     for common/onetwo/
      WRITE(*,*)'maximal number of bin boundaries in the small radius'
     +//' direction'
      WRITE(*,*)'for power and CD radial profiles'//
     +'Power and current are tabulated at (NR-1) bin centers :'//
     +'NRA = ',NRA
       
c************************************************************
c     for common/rho/
      
      WRITE(*,*)'npsi4=npsi+4 = ',npsi4
      WRITE(*,*)'nteta14=nteta1+4 = ',nteta14
      WRITE(*,*)'nzy=nteta1+4 = ',nzy
      endif
      
      if (nzy.ne.(max(npsi,nteta1)+4)) then
      if(myrank.eq.0) then
      WRITE(*,*)'nzy.ne.(max(npsi,nteta1)+4)'
      WRITE(*,*)'It should be nzy=max(npsi,nteta1)+4'
      WRITE(*,*)'npsi = ',npsi
      WRITE(*,*)'nteta1 = ',nteta1
      WRITE(*,*)'Please change  nzy in param.i '//
     +'for nzy=nteta1+4 or nzy=npsi+4 '//
     +'and recompile the code'
      endif
      stop ! stop in all cores
      endif

      if(myrank.eq.0) then
c************************************************************
c     for common/six/
      WRITE(*,*)'max number of points in arrays with '// 
     +'plasma density, temperature, zeff. tpop, and vflow: '//
     +' ndensa = ',ndensa
      WRITE(*,*)'ndens4a=ndensa+4 = ',ndens4a

c************************************************************
c     for common/write/

      WRITE(*,*)'max number of output points along each ray: '//
     +'nrelta =',nrelta

      WRITE(*,*)'n_relt_harma = ',n_relt_harma
      WRITE(*,*)'n_relt_harm1a = ',n_relt_harm1a
      WRITE(*,*)'n_relt_harm2a = ',n_relt_harm2a
c************************************************************
c     for common/scatnper/ (the data for the n_perp scattering )

      WRITE(*,*)'namber of points in rhoscat(nscat_n) and '//
     +'scatd(0:nscat_n): nscat_n = ',nscat_n
     
c*************************************************************
c     for common/dskincomm/ the data for reading dskin file

      WRITE(*,*)'to read distribution function f(iya,jx,lrz,ngen) '//
     +'written by CQL3D code'
      
      WRITE(*,*)'number of the points in the pitch angle mesh: iya = ',
     +iya

      WRITE(*,*)'iya1=iya+1 = ',iya1
      
c***************************************************************
c     for common/emissa/
     
      WRITE(*,*)'nrelta4=nrelta+4 = ', nrelta4
 
c*************************************************************
c     for common /output/
     
      WRITE(*,*)'n_plot_dispa =',n_plot_dispa
      WRITE(*,*)'n_plot_disp_colda = ',n_plot_disp_colda

      WRITE(*,*)'m_r_nperp_a = ',m_r_nperp_a 
      WRITE(*,*)'m_i_nperp_a = ',m_i_nperp_a 
      WRITE(*,*)'n_contour_plot_disp_a = ',n_contour_plot_disp_a 

c*************************************************************
c-----for hot plasma roots n_hot_roots_a max number of hot plasma root
      WRITE(*,*)'n_hot_roots_a = ',n_hot_roots_a
c******************************************************************
c     for small output step near EC resonance points for power calculation 
c     It will be used for data at namelist /numercl/ and common /one/
      WRITE(*,*)'n_power_switch_resonance_a = ',
     +n_power_switch_resonance_a

c******************************************************************
c     for toray EC launch
      WRITE(*,*)'gzonemax = ',gzonemax 


c********************************************************************
c     to read wall and limiter coordinates (r,z)
c     for writencdf.i:
 
      WRITE(*,*)'maximal numbes of wall points: n_wall_a = ',n_wall_a
 
      WRITE(*,*)'maximal numbes of limiter points: n_limiter_a = ',
     +n_limiter_a

      WRITE(*,*)' maximal number of limiters: max_limiters_a = ',
     +max_limiters_a

c     maximal number of wall points having the given poloidal angle theta
      
      WRITE(*,*)'maximal number of wall points having '//
     +' the given poloidal '//
     +' angle theta: n_rho_wall_a = ', n_rho_wall_a

c***********************************************************************
c     to read normalized exponential density falls at poloidal mesh
c     for edge_prof_nml,i
c     n_pol_edge_dens_a is a maximal number of mesh points
    
      WRITE(*,*)'to read normalized exponential density falls'//
     +' at poloidal mesh for edge_prof_nml,i '
      WRITE(*,*)'maximal number of mesh points: n_pol_edge_dens_a = ',
     +n_pol_edge_dens_a
c********************************************************************
c     to create fine meshes of wall coordinates with additional points
      
      WRITE(*,*)'to create fine meshes of wall coordinates '//
     +'with additional points: n_wall_add_a = ',n_wall_add_a
c*********************************************************************
c     for (R,Z) meshes rr_add, zz_add at the poloidal plane
c     for the density fall near the wall
      
c      write(*,*)'max number of r mesh points at the poloidal plane: '//
c     & 'nxeqd_add_a = ', nxeqd_add_a 
c      write(*,*)'max number of z mesh points at the poloidal plane: '// 
c     & 'nyeqd_add_a = ', nyeqd_add_a 
c--------------------------------------------------------------------
      WRITE(*,*)'*********** END print screen code parameters ********'
      endif ! rank=0
      
      return
      end







c***********************************************************************
c***********************************************************************

c                       GNU GENERAL PUBLIC LICENSE
c                          Version 3, 29 June 2007
c   
c    Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
c    Everyone is permitted to copy and distribute verbatim copies
c    of this license document, but changing it is not allowed.
c   
c                               Preamble
c   
c     The GNU General Public License is a free, copyleft license for
c   software and other kinds of works.
c   
c     The licenses for most software and other practical works are designed
c   to take away your freedom to share and change the works.  By contrast,
c   the GNU General Public License is intended to guarantee your freedom to
c   share and change all versions of a program--to make sure it remains free
c   software for all its users.  We, the Free Software Foundation, use the
c   GNU General Public License for most of our software; it applies also to
c   any other work released this way by its authors.  You can apply it to
c   your programs, too.
c   
c     When we speak of free software, we are referring to freedom, not
c   price.  Our General Public Licenses are designed to make sure that you
c   have the freedom to distribute copies of free software (and charge for
c   them if you wish), that you receive source code or can get it if you
c   want it, that you can change the software or use pieces of it in new
c   free programs, and that you know you can do these things.
c   
c     To protect your rights, we need to prevent others from denying you
c   these rights or asking you to surrender the rights.  Therefore, you have
c   certain responsibilities if you distribute copies of the software, or if
c   you modify it: responsibilities to respect the freedom of others.
c   
c     For example, if you distribute copies of such a program, whether
c   gratis or for a fee, you must pass on to the recipients the same
c   freedoms that you received.  You must make sure that they, too, receive
c   or can get the source code.  And you must show them these terms so they
c   know their rights.
c   
c     Developers that use the GNU GPL protect your rights with two steps:
c   (1) assert copyright on the software, and (2) offer you this License
c   giving you legal permission to copy, distribute and/or modify it.
c   
c     For the developers' and authors' protection, the GPL clearly explains
c   that there is no warranty for this free software.  For both users' and
c   authors' sake, the GPL requires that modified versions be marked as
c   changed, so that their problems will not be attributed erroneously to
c   authors of previous versions.
c   
c     Some devices are designed to deny users access to install or run
c   modified versions of the software inside them, although the manufacturer
c   can do so.  This is fundamentally incompatible with the aim of
c   protecting users' freedom to change the software.  The systematic
c   pattern of such abuse occurs in the area of products for individuals to
c   use, which is precisely where it is most unacceptable.  Therefore, we
c   have designed this version of the GPL to prohibit the practice for those
c   products.  If such problems arise substantially in other domains, we
c   stand ready to extend this provision to those domains in future versions
c   of the GPL, as needed to protect the freedom of users.
c   
c     Finally, every program is threatened constantly by software patents.
c   States should not allow patents to restrict development and use of
c   software on general-purpose computers, but in those that do, we wish to
c   avoid the special danger that patents applied to a free program could
c   make it effectively proprietary.  To prevent this, the GPL assures that
c   patents cannot be used to render the program non-free.
c   
c     The precise terms and conditions for copying, distribution and
c   modification follow.
c   
c                          TERMS AND CONDITIONS
c   
c     0. Definitions.
c   
c     "This License" refers to version 3 of the GNU General Public License.
c   
c     "Copyright" also means copyright-like laws that apply to other kinds of
c   works, such as semiconductor masks.
c   
c     "The Program" refers to any copyrightable work licensed under this
c   License.  Each licensee is addressed as "you".  "Licensees" and
c   "recipients" may be individuals or organizations.
c   
c     To "modify" a work means to copy from or adapt all or part of the work
c   in a fashion requiring copyright permission, other than the making of an
c   exact copy.  The resulting work is called a "modified version" of the
c   earlier work or a work "based on" the earlier work.
c   
c     A "covered work" means either the unmodified Program or a work based
c   on the Program.
c   
c     To "propagate" a work means to do anything with it that, without
c   permission, would make you directly or secondarily liable for
c   infringement under applicable copyright law, except executing it on a
c   computer or modifying a private copy.  Propagation includes copying,
c   distribution (with or without modification), making available to the
c   public, and in some countries other activities as well.
c   
c     To "convey" a work means any kind of propagation that enables other
c   parties to make or receive copies.  Mere interaction with a user through
c   a computer network, with no transfer of a copy, is not conveying.
c   
c     An interactive user interface displays "Appropriate Legal Notices"
c   to the extent that it includes a convenient and prominently visible
c   feature that (1) displays an appropriate copyright notice, and (2)
c   tells the user that there is no warranty for the work (except to the
c   extent that warranties are provided), that licensees may convey the
c   work under this License, and how to view a copy of this License.  If
c   the interface presents a list of user commands or options, such as a
c   menu, a prominent item in the list meets this criterion.
c   
c     1. Source Code.
c   
c     The "source code" for a work means the preferred form of the work
c   for making modifications to it.  "Object code" means any non-source
c   form of a work.
c   
c     A "Standard Interface" means an interface that either is an official
c   standard defined by a recognized standards body, or, in the case of
c   interfaces specified for a particular programming language, one that
c   is widely used among developers working in that language.
c   
c     The "System Libraries" of an executable work include anything, other
c   than the work as a whole, that (a) is included in the normal form of
c   packaging a Major Component, but which is not part of that Major
c   Component, and (b) serves only to enable use of the work with that
c   Major Component, or to implement a Standard Interface for which an
c   implementation is available to the public in source code form.  A
c   "Major Component", in this context, means a major essential component
c   (kernel, window system, and so on) of the specific operating system
c   (if any) on which the executable work runs, or a compiler used to
c   produce the work, or an object code interpreter used to run it.
c   
c     The "Corresponding Source" for a work in object code form means all
c   the source code needed to generate, install, and (for an executable
c   work) run the object code and to modify the work, including scripts to
c   control those activities.  However, it does not include the work's
c   System Libraries, or general-purpose tools or generally available free
c   programs which are used unmodified in performing those activities but
c   which are not part of the work.  For example, Corresponding Source
c   includes interface definition files associated with source files for
c   the work, and the source code for shared libraries and dynamically
c   linked subprograms that the work is specifically designed to require,
c   such as by intimate data communication or control flow between those
c   subprograms and other parts of the work.
c   
c     The Corresponding Source need not include anything that users
c   can regenerate automatically from other parts of the Corresponding
c   Source.
c   
c     The Corresponding Source for a work in source code form is that
c   same work.
c   
c     2. Basic Permissions.
c   
c     All rights granted under this License are granted for the term of
c   copyright on the Program, and are irrevocable provided the stated
c   conditions are met.  This License explicitly affirms your unlimited
c   permission to run the unmodified Program.  The output from running a
c   covered work is covered by this License only if the output, given its
c   content, constitutes a covered work.  This License acknowledges your
c   rights of fair use or other equivalent, as provided by copyright law.
c   
c     You may make, run and propagate covered works that you do not
c   convey, without conditions so long as your license otherwise remains
c   in force.  You may convey covered works to others for the sole purpose
c   of having them make modifications exclusively for you, or provide you
c   with facilities for running those works, provided that you comply with
c   the terms of this License in conveying all material for which you do
c   not control copyright.  Those thus making or running the covered works
c   for you must do so exclusively on your behalf, under your direction
c   and control, on terms that prohibit them from making any copies of
c   your copyrighted material outside their relationship with you.
c   
c     Conveying under any other circumstances is permitted solely under
c   the conditions stated below.  Sublicensing is not allowed; section 10
c   makes it unnecessary.
c   
c     3. Protecting Users' Legal Rights From Anti-Circumvention Law.
c   
c     No covered work shall be deemed part of an effective technological
c   measure under any applicable law fulfilling obligations under article
c   11 of the WIPO copyright treaty adopted on 20 December 1996, or
c   similar laws prohibiting or restricting circumvention of such
c   measures.
c   
c     When you convey a covered work, you waive any legal power to forbid
c   circumvention of technological measures to the extent such circumvention
c   is effected by exercising rights under this License with respect to
c   the covered work, and you disclaim any intention to limit operation or
c   modification of the work as a means of enforcing, against the work's
c   users, your or third parties' legal rights to forbid circumvention of
c   technological measures.
c   
c     4. Conveying Verbatim Copies.
c   
c     You may convey verbatim copies of the Program's source code as you
c   receive it, in any medium, provided that you conspicuously and
c   appropriately publish on each copy an appropriate copyright notice;
c   keep intact all notices stating that this License and any
c   non-permissive terms added in accord with section 7 apply to the code;
c   keep intact all notices of the absence of any warranty; and give all
c   recipients a copy of this License along with the Program.
c   
c     You may charge any price or no price for each copy that you convey,
c   and you may offer support or warranty protection for a fee.
c   
c     5. Conveying Modified Source Versions.
c   
c     You may convey a work based on the Program, or the modifications to
c   produce it from the Program, in the form of source code under the
c   terms of section 4, provided that you also meet all of these conditions:
c   
c       a) The work must carry prominent notices stating that you modified
c       it, and giving a relevant date.
c   
c       b) The work must carry prominent notices stating that it is
c       released under this License and any conditions added under section
c       7.  This requirement modifies the requirement in section 4 to
c       "keep intact all notices".
c   
c       c) You must license the entire work, as a whole, under this
c       License to anyone who comes into possession of a copy.  This
c       License will therefore apply, along with any applicable section 7
c       additional terms, to the whole of the work, and all its parts,
c       regardless of how they are packaged.  This License gives no
c       permission to license the work in any other way, but it does not
c       invalidate such permission if you have separately received it.
c   
c       d) If the work has interactive user interfaces, each must display
c       Appropriate Legal Notices; however, if the Program has interactive
c       interfaces that do not display Appropriate Legal Notices, your
c       work need not make them do so.
c   
c     A compilation of a covered work with other separate and independent
c   works, which are not by their nature extensions of the covered work,
c   and which are not combined with it such as to form a larger program,
c   in or on a volume of a storage or distribution medium, is called an
c   "aggregate" if the compilation and its resulting copyright are not
c   used to limit the access or legal rights of the compilation's users
c   beyond what the individual works permit.  Inclusion of a covered work
c   in an aggregate does not cause this License to apply to the other
c   parts of the aggregate.
c   
c     6. Conveying Non-Source Forms.
c   
c     You may convey a covered work in object code form under the terms
c   of sections 4 and 5, provided that you also convey the
c   machine-readable Corresponding Source under the terms of this License,
c   in one of these ways:
c   
c       a) Convey the object code in, or embodied in, a physical product
c       (including a physical distribution medium), accompanied by the
c       Corresponding Source fixed on a durable physical medium
c       customarily used for software interchange.
c   
c       b) Convey the object code in, or embodied in, a physical product
c       (including a physical distribution medium), accompanied by a
c       written offer, valid for at least three years and valid for as
c       long as you offer spare parts or customer support for that product
c       model, to give anyone who possesses the object code either (1) a
c       copy of the Corresponding Source for all the software in the
c       product that is covered by this License, on a durable physical
c       medium customarily used for software interchange, for a price no
c       more than your reasonable cost of physically performing this
c       conveying of source, or (2) access to copy the
c       Corresponding Source from a network server at no charge.
c   
c       c) Convey individual copies of the object code with a copy of the
c       written offer to provide the Corresponding Source.  This
c       alternative is allowed only occasionally and noncommercially, and
c       only if you received the object code with such an offer, in accord
c       with subsection 6b.
c   
c       d) Convey the object code by offering access from a designated
c       place (gratis or for a charge), and offer equivalent access to the
c       Corresponding Source in the same way through the same place at no
c       further charge.  You need not require recipients to copy the
c       Corresponding Source along with the object code.  If the place to
c       copy the object code is a network server, the Corresponding Source
c       may be on a different server (operated by you or a third party)
c       that supports equivalent copying facilities, provided you maintain
c       clear directions next to the object code saying where to find the
c       Corresponding Source.  Regardless of what server hosts the
c       Corresponding Source, you remain obligated to ensure that it is
c       available for as long as needed to satisfy these requirements.
c   
c       e) Convey the object code using peer-to-peer transmission, provided
c       you inform other peers where the object code and Corresponding
c       Source of the work are being offered to the general public at no
c       charge under subsection 6d.
c   
c     A separable portion of the object code, whose source code is excluded
c   from the Corresponding Source as a System Library, need not be
c   included in conveying the object code work.
c   
c     A "User Product" is either (1) a "consumer product", which means any
c   tangible personal property which is normally used for personal, family,
c   or household purposes, or (2) anything designed or sold for incorporation
c   into a dwelling.  In determining whether a product is a consumer product,
c   doubtful cases shall be resolved in favor of coverage.  For a particular
c   product received by a particular user, "normally used" refers to a
c   typical or common use of that class of product, regardless of the status
c   of the particular user or of the way in which the particular user
c   actually uses, or expects or is expected to use, the product.  A product
c   is a consumer product regardless of whether the product has substantial
c   commercial, industrial or non-consumer uses, unless such uses represent
c   the only significant mode of use of the product.
c   
c     "Installation Information" for a User Product means any methods,
c   procedures, authorization keys, or other information required to install
c   and execute modified versions of a covered work in that User Product from
c   a modified version of its Corresponding Source.  The information must
c   suffice to ensure that the continued functioning of the modified object
c   code is in no case prevented or interfered with solely because
c   modification has been made.
c   
c     If you convey an object code work under this section in, or with, or
c   specifically for use in, a User Product, and the conveying occurs as
c   part of a transaction in which the right of possession and use of the
c   User Product is transferred to the recipient in perpetuity or for a
c   fixed term (regardless of how the transaction is characterized), the
c   Corresponding Source conveyed under this section must be accompanied
c   by the Installation Information.  But this requirement does not apply
c   if neither you nor any third party retains the ability to install
c   modified object code on the User Product (for example, the work has
c   been installed in ROM).
c   
c     The requirement to provide Installation Information does not include a
c   requirement to continue to provide support service, warranty, or updates
c   for a work that has been modified or installed by the recipient, or for
c   the User Product in which it has been modified or installed.  Access to a
c   network may be denied when the modification itself materially and
c   adversely affects the operation of the network or violates the rules and
c   protocols for communication across the network.
c   
c     Corresponding Source conveyed, and Installation Information provided,
c   in accord with this section must be in a format that is publicly
c   documented (and with an implementation available to the public in
c   source code form), and must require no special password or key for
c   unpacking, reading or copying.
c   
c     7. Additional Terms.
c   
c     "Additional permissions" are terms that supplement the terms of this
c   License by making exceptions from one or more of its conditions.
c   Additional permissions that are applicable to the entire Program shall
c   be treated as though they were included in this License, to the extent
c   that they are valid under applicable law.  If additional permissions
c   apply only to part of the Program, that part may be used separately
c   under those permissions, but the entire Program remains governed by
c   this License without regard to the additional permissions.
c   
c     When you convey a copy of a covered work, you may at your option
c   remove any additional permissions from that copy, or from any part of
c   it.  (Additional permissions may be written to require their own
c   removal in certain cases when you modify the work.)  You may place
c   additional permissions on material, added by you to a covered work,
c   for which you have or can give appropriate copyright permission.
c   
c     Notwithstanding any other provision of this License, for material you
c   add to a covered work, you may (if authorized by the copyright holders of
c   that material) supplement the terms of this License with terms:
c   
c       a) Disclaiming warranty or limiting liability differently from the
c       terms of sections 15 and 16 of this License; or
c   
c       b) Requiring preservation of specified reasonable legal notices or
c       author attributions in that material or in the Appropriate Legal
c       Notices displayed by works containing it; or
c   
c       c) Prohibiting misrepresentation of the origin of that material, or
c       requiring that modified versions of such material be marked in
c       reasonable ways as different from the original version; or
c   
c       d) Limiting the use for publicity purposes of names of licensors or
c       authors of the material; or
c   
c       e) Declining to grant rights under trademark law for use of some
c       trade names, trademarks, or service marks; or
c   
c       f) Requiring indemnification of licensors and authors of that
c       material by anyone who conveys the material (or modified versions of
c       it) with contractual assumptions of liability to the recipient, for
c       any liability that these contractual assumptions directly impose on
c       those licensors and authors.
c   
c     All other non-permissive additional terms are considered "further
c   restrictions" within the meaning of section 10.  If the Program as you
c   received it, or any part of it, contains a notice stating that it is
c   governed by this License along with a term that is a further
c   restriction, you may remove that term.  If a license document contains
c   a further restriction but permits relicensing or conveying under this
c   License, you may add to a covered work material governed by the terms
c   of that license document, provided that the further restriction does
c   not survive such relicensing or conveying.
c   
c     If you add terms to a covered work in accord with this section, you
c   must place, in the relevant source files, a statement of the
c   additional terms that apply to those files, or a notice indicating
c   where to find the applicable terms.
c   
c     Additional terms, permissive or non-permissive, may be stated in the
c   form of a separately written license, or stated as exceptions;
c   the above requirements apply either way.
c   
c     8. Termination.
c   
c     You may not propagate or modify a covered work except as expressly
c   provided under this License.  Any attempt otherwise to propagate or
c   modify it is void, and will automatically terminate your rights under
c   this License (including any patent licenses granted under the third
c   paragraph of section 11).
c   
c     However, if you cease all violation of this License, then your
c   license from a particular copyright holder is reinstated (a)
c   provisionally, unless and until the copyright holder explicitly and
c   finally terminates your license, and (b) permanently, if the copyright
c   holder fails to notify you of the violation by some reasonable means
c   prior to 60 days after the cessation.
c   
c     Moreover, your license from a particular copyright holder is
c   reinstated permanently if the copyright holder notifies you of the
c   violation by some reasonable means, this is the first time you have
c   received notice of violation of this License (for any work) from that
c   copyright holder, and you cure the violation prior to 30 days after
c   your receipt of the notice.
c   
c     Termination of your rights under this section does not terminate the
c   licenses of parties who have received copies or rights from you under
c   this License.  If your rights have been terminated and not permanently
c   reinstated, you do not qualify to receive new licenses for the same
c   material under section 10.
c   
c     9. Acceptance Not Required for Having Copies.
c   
c     You are not required to accept this License in order to receive or
c   run a copy of the Program.  Ancillary propagation of a covered work
c   occurring solely as a consequence of using peer-to-peer transmission
c   to receive a copy likewise does not require acceptance.  However,
c   nothing other than this License grants you permission to propagate or
c   modify any covered work.  These actions infringe copyright if you do
c   not accept this License.  Therefore, by modifying or propagating a
c   covered work, you indicate your acceptance of this License to do so.
c   
c     10. Automatic Licensing of Downstream Recipients.
c   
c     Each time you convey a covered work, the recipient automatically
c   receives a license from the original licensors, to run, modify and
c   propagate that work, subject to this License.  You are not responsible
c   for enforcing compliance by third parties with this License.
c   
c     An "entity transaction" is a transaction transferring control of an
c   organization, or substantially all assets of one, or subdividing an
c   organization, or merging organizations.  If propagation of a covered
c   work results from an entity transaction, each party to that
c   transaction who receives a copy of the work also receives whatever
c   licenses to the work the party's predecessor in interest had or could
c   give under the previous paragraph, plus a right to possession of the
c   Corresponding Source of the work from the predecessor in interest, if
c   the predecessor has it or can get it with reasonable efforts.
c   
c     You may not impose any further restrictions on the exercise of the
c   rights granted or affirmed under this License.  For example, you may
c   not impose a license fee, royalty, or other charge for exercise of
c   rights granted under this License, and you may not initiate litigation
c   (including a cross-claim or counterclaim in a lawsuit) alleging that
c   any patent claim is infringed by making, using, selling, offering for
c   sale, or importing the Program or any portion of it.
c   
c     11. Patents.
c   
c     A "contributor" is a copyright holder who authorizes use under this
c   License of the Program or a work on which the Program is based.  The
c   work thus licensed is called the contributor's "contributor version".
c   
c     A contributor's "essential patent claims" are all patent claims
c   owned or controlled by the contributor, whether already acquired or
c   hereafter acquired, that would be infringed by some manner, permitted
c   by this License, of making, using, or selling its contributor version,
c   but do not include claims that would be infringed only as a
c   consequence of further modification of the contributor version.  For
c   purposes of this definition, "control" includes the right to grant
c   patent sublicenses in a manner consistent with the requirements of
c   this License.
c   
c     Each contributor grants you a non-exclusive, worldwide, royalty-free
c   patent license under the contributor's essential patent claims, to
c   make, use, sell, offer for sale, import and otherwise run, modify and
c   propagate the contents of its contributor version.
c   
c     In the following three paragraphs, a "patent license" is any express
c   agreement or commitment, however denominated, not to enforce a patent
c   (such as an express permission to practice a patent or covenant not to
c   sue for patent infringement).  To "grant" such a patent license to a
c   party means to make such an agreement or commitment not to enforce a
c   patent against the party.
c   
c     If you convey a covered work, knowingly relying on a patent license,
c   and the Corresponding Source of the work is not available for anyone
c   to copy, free of charge and under the terms of this License, through a
c   publicly available network server or other readily accessible means,
c   then you must either (1) cause the Corresponding Source to be so
c   available, or (2) arrange to deprive yourself of the benefit of the
c   patent license for this particular work, or (3) arrange, in a manner
c   consistent with the requirements of this License, to extend the patent
c   license to downstream recipients.  "Knowingly relying" means you have
c   actual knowledge that, but for the patent license, your conveying the
c   covered work in a country, or your recipient's use of the covered work
c   in a country, would infringe one or more identifiable patents in that
c   country that you have reason to believe are valid.
c   
c     If, pursuant to or in connection with a single transaction or
c   arrangement, you convey, or propagate by procuring conveyance of, a
c   covered work, and grant a patent license to some of the parties
c   receiving the covered work authorizing them to use, propagate, modify
c   or convey a specific copy of the covered work, then the patent license
c   you grant is automatically extended to all recipients of the covered
c   work and works based on it.
c   
c     A patent license is "discriminatory" if it does not include within
c   the scope of its coverage, prohibits the exercise of, or is
c   conditioned on the non-exercise of one or more of the rights that are
c   specifically granted under this License.  You may not convey a covered
c   work if you are a party to an arrangement with a third party that is
c   in the business of distributing software, under which you make payment
c   to the third party based on the extent of your activity of conveying
c   the work, and under which the third party grants, to any of the
c   parties who would receive the covered work from you, a discriminatory
c   patent license (a) in connection with copies of the covered work
c   conveyed by you (or copies made from those copies), or (b) primarily
c   for and in connection with specific products or compilations that
c   contain the covered work, unless you entered into that arrangement,
c   or that patent license was granted, prior to 28 March 2007.
c   
c     Nothing in this License shall be construed as excluding or limiting
c   any implied license or other defenses to infringement that may
c   otherwise be available to you under applicable patent law.
c   
c     12. No Surrender of Others' Freedom.
c   
c     If conditions are imposed on you (whether by court order, agreement or
c   otherwise) that contradict the conditions of this License, they do not
c   excuse you from the conditions of this License.  If you cannot convey a
c   covered work so as to satisfy simultaneously your obligations under this
c   License and any other pertinent obligations, then as a consequence you may
c   not convey it at all.  For example, if you agree to terms that obligate you
c   to collect a royalty for further conveying from those to whom you convey
c   the Program, the only way you could satisfy both those terms and this
c   License would be to refrain entirely from conveying the Program.
c   
c     13. Use with the GNU Affero General Public License.
c   
c     Notwithstanding any other provision of this License, you have
c   permission to link or combine any covered work with a work licensed
c   under version 3 of the GNU Affero General Public License into a single
c   combined work, and to convey the resulting work.  The terms of this
c   License will continue to apply to the part which is the covered work,
c   but the special requirements of the GNU Affero General Public License,
c   section 13, concerning interaction through a network will apply to the
c   combination as such.
c   
c     14. Revised Versions of this License.
c   
c     The Free Software Foundation may publish revised and/or new versions of
c   the GNU General Public License from time to time.  Such new versions will
c   be similar in spirit to the present version, but may differ in detail to
c   address new problems or concerns.
c   
c     Each version is given a distinguishing version number.  If the
c   Program specifies that a certain numbered version of the GNU General
c   Public License "or any later version" applies to it, you have the
c   option of following the terms and conditions either of that numbered
c   version or of any later version published by the Free Software
c   Foundation.  If the Program does not specify a version number of the
c   GNU General Public License, you may choose any version ever published
c   by the Free Software Foundation.
c   
c     If the Program specifies that a proxy can decide which future
c   versions of the GNU General Public License can be used, that proxy's
c   public statement of acceptance of a version permanently authorizes you
c   to choose that version for the Program.
c   
c     Later license versions may give you additional or different
c   permissions.  However, no additional obligations are imposed on any
c   author or copyright holder as a result of your choosing to follow a
c   later version.
c   
c     15. Disclaimer of Warranty.
c   
c     THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
c   APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
c   HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
c   OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
c   THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
c   PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
c   IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
c   ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
c   
c     16. Limitation of Liability.
c   
c     IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
c   WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS
c   THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY
c   GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE
c   USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF
c   DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD
c   PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
c   EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
c   SUCH DAMAGES.
c   
c     17. Interpretation of Sections 15 and 16.
c   
c     If the disclaimer of warranty and limitation of liability provided
c   above cannot be given local legal effect according to their terms,
c   reviewing courts shall apply local law that most closely approximates
c   an absolute waiver of all civil liability in connection with the
c   Program, unless a warranty or assumption of liability accompanies a
c   copy of the Program in return for a fee.
c   
c                        END OF TERMS AND CONDITIONS
c   
c               How to Apply These Terms to Your New Programs
c   
c     If you develop a new program, and you want it to be of the greatest
c   possible use to the public, the best way to achieve this is to make it
c   free software which everyone can redistribute and change under these terms.
c   
c     To do so, attach the following notices to the program.  It is safest
c   to attach them to the start of each source file to most effectively
c   state the exclusion of warranty; and each file should have at least
c   the "copyright" line and a pointer to where the full notice is found.
c   
c       <one line to give the program's name and a brief idea of what it does.>
c       Copyright (C) <year>  <name of author>
c   
c       This program is free software: you can redistribute it and/or modify
c       it under the terms of the GNU General Public License as published by
c       the Free Software Foundation, either version 3 of the License, or
c       (at your option) any later version.
c   
c       This program is distributed in the hope that it will be useful,
c       but WITHOUT ANY WARRANTY; without even the implied warranty of
c       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c       GNU General Public License for more details.
c   
c       You should have received a copy of the GNU General Public License
c       along with this program.  If not, see <http://www.gnu.org/licenses/>.
c   
c   Also add information on how to contact you by electronic and paper mail.
c   
c     If the program does terminal interaction, make it output a short
c   notice like this when it starts in an interactive mode:
c   
c       <program>  Copyright (C) <year>  <name of author>
c       This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
c       This is free software, and you are welcome to redistribute it
c       under certain conditions; type `show c' for details.
c   
c   The hypothetical commands `show w' and `show c' should show the appropriate
c   parts of the General Public License.  Of course, your program's commands
c   might be different; for a GUI interface, you would use an "about box".
c   
c     You should also get your employer (if you work as a programmer) or school,
c   if any, to sign a "copyright disclaimer" for the program, if necessary.
c   For more information on this, and how to apply and follow the GNU GPL, see
c   <http://www.gnu.org/licenses/>.
c   
c     The GNU General Public License does not permit incorporating your program
c   into proprietary programs.  If your program is a subroutine library, you
c   may consider it more useful to permit linking proprietary applications with
c   the library.  If this is what you want to do, use the GNU Lesser General
c   Public License instead of this License.  But first, please read
c   <http://www.gnu.org/philosophy/why-not-lgpl.html>.
c   
c***********************************************************************
c***********************************************************************
