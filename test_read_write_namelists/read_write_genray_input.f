

c======================================================================
c NOTE: When using with prepare_genray_input (for Plasma State),comment 
c       out if(partner..)-clause surrounding 'call read_transport_prof'.
c       This gives read_write_genray_input_prep.f, and it should be
c       kept in sync with read_write_genray_input.f.
c================s======================================================

     
c      program prepare_genray_input
      subroutine prepare_genray_input
!     (~84 lines to end)
!     Prepares new genray.in namelist file for genray, using
!     several files from the genray source, and the old 
!     genray.in input file.

      implicit none

      include 'param.i'
c----------------------------------------
c     declaration and storage of all namelists variables
c-----------------------------------
      include 'cone_nml.i'
      include 'emissa_nml.i'
      include 'grill_nml.i'
      include 'ions_nml.i'
      include 'one_nml.i'
      include 'onetwo_nml.i'
      include 'output_nml.i'
      include 'scatnper_nml.i'    
      include 'six_nml.i'   
      include 'adj_nml.i'
cSAP090203
      include 'edge_prof_nml.i'

cSAP111028
c      include 'lsc_approach.i'
      include 'lsc_approach_nml.i'
c--------------------------------------------------------------
c     declaration of namelist variables which
c     are used in subroutine dinit_mr and
c     are not used in any common block
c--------------------------------------------------------------
      include 'dinit_nml.i' 
c--------------------------------------------------------------
cBH111027      include 'rkutta.i' ! contains common /rkutta/ prmt(9)
c--------------------------------------------------------------
c     list of all namelists
c--------------------------------------------------------------
      include 'name.i'    ! list of all namelists

c--------------------------------------------------------------
      integer ndim,nray

      integer i,j,k

cSAP080711
c-----genray input file: genray.in or genray.dat
      character*10 genray_in_dat

      save ndim  ! to give ndim to
                 ! subroutine write_all_namelists
c--------------------------------------------------------------
c     the creation of default input data like in
c     genray.dat files
c--------------------------------------------------------------    
      call default_in
c--------------------------------------------------------------
c     check if input file is genray.in, in which case change
c     default data to genray.in file (MKSA) format
c--------------------------------------------------------------  
      call transform_input_data_to_MKSA
c---------------------------------------------------------------
c     Read all namelists from input genray.in or genray.dat file
c     Here genray_in_dat is an output argument,
c     As prepare_genray_input uses only genray.in, it will have
c     genray_in_dat='genray.in'
c---------------------------------------------------------------
      call read_all_namelists(genray_in_dat,ndim,nray)
      write(*,*)'after first read_all_namelists, ndim,nray=',ndim,nray
c---------------------------------------------------------------------------
c     writes all namelists to genray.in file
c----------------------------------------------------------------------------
      call write_all_namelists(ndim)
cSAP080730
       write(*,*)'###after write_all_namelists'
c----------------------------------------------------------------------
c     To check the new created genray.in file it reads new genray.in
c-----------------------------------------------------------------------  
      call read_all_namelists(genray_in_dat,ndim,nray)
      write(*,*)'after new read_all_namelists, ndim,nray=',ndim,nray
 
21       format(5e16.9)

      stop
      end    !!!!!!!prepare_genray_input!!!!!!!!!!!!!!!!!!!!!!!!!!!





      subroutine read_all_namelists(genray_in_dat,ndim,nray)
!     [end at line ~1483]
c--------------------------------------------------------------
c     reads all namelists from genray.in or genray.dat input file
c--------------------------------------------------------------
c-----output:
c     genray_in_dat is the name of nml file 'genray.in' or 'genray.dat'
c     found in the PWD.
c     ndim is a number of ray-tracing equations (=6 at isolv=1)
c                                                (=5 at isolv=2) 
c     nray is a number of all rays. 
c          It will be calculated for
c          1)istartr=1 and
c            raypatt='diskdisk' or 'diskbeam' 
      implicit none

      include 'param.i'
c----------------------------------------
c     declaration and storage of all namelists variables
c-----------------------------------
      include 'cone_nml.i'
      include 'emissa_nml.i'
      include 'grill_nml.i'
      include 'ions_nml.i'
      include 'one_nml.i'
      include 'onetwo_nml.i'
      include 'output_nml.i'
      include 'scatnper_nml.i'    
      include 'six_nml.i'  
      include 'adj_nml.i'
cSAP090203
      include 'edge_prof_nml.i'
cSAP091222
      include 'lsc_approach_nml.i'
c--------------------------------------------------------------
c     declaration of namelist variables which
c     are used in subroutine dinit_mr and
c     are not used in any common block
c--------------------------------------------------------------  
      include 'dinit_nml.i' 
c--------------------------------------------------------------
      include 'rkutta.i' ! contains common /rkutta/ prmt(9)
c--------------------------------------------------------------
c     list of all namelists
c--------------------------------------------------------------
      include 'name.i'    ! list of all namelists
c---------------------------------------------------------------------------
c-----output
      integer ndim,nray
      character*10 genray_in_dat
c-----local
      integer kode,icheck,i,k,j,nbulk1,i1,i_unit,
     &i_genray_in_transformation
      real*8 powtott, nnkprmxn, nthinmxn, psi0, rho0,
     &rhom(ndensa),h,rho,
     &te0(nbulka),teb(nbulka),
     &prof2_uniform(ndensa,nbulka)
c-----external
      integer length_char
      real*8 fpsi, rhopsi

      write(*,*)'in subroutine read_all_namelists 1'

      i_unit=1

      i_genray_in_transformation=0

cSAP100914
c      open(i_unit,file = 'genray.in',status='old',iostat=kode)
      open(i_unit,file = 'genray.dat',status='old',iostat=kode)
cSAP100914
c      write(*,*)'in read_all_namelists after open file=genray.in kode',
c     &kode
      write(*,*)'in read_all_namelists after open file=genray.dat kode',
     &kode
cSAP101004  
      if (kode.eq.0) genray_in_dat='genray.dat'  !the name of used input file

cSAP100914
c      if (kode.eq.0) then
      if (kode.ne.0) then
        i_genray_in_transformation=1 !input data file genray.in
                                     !in MKSA system
cSAP100914                                     
        open(i_unit,file='genray.in',status='old',iostat=kode)
        if (kode.ne.0) then
            write(*,*)' read_all_namelists:'
            write(*,*)' Neither genray.in or genray.dat r present'
            stop
        endif
cSAP080711
        genray_in_dat='genray.in'        !the name of used input file
c---------------------------------------------------------------
c       transform default data to MKSA system
c---------------------------------------------------------------
c        call default_in(i_genray_in_transformation) 
      endif

cSAP100914
c      if (kode.ne.0) then
c         open(i_unit,file='genray.dat',status='old',iostat=kode)
ccSAP080711
c         genray_in_dat='genray.dat'      !the name of used input file
c         open(i_unit,file='genray.in',status='old',iostat=kode)
c         genray_in_dat='genray.in'      !the name of used input fil
c         if (kode.ne.0) then
c            write(*,*)' prepare_genray_input:'
c            write(*,*)' Neither genray.in or genray.dat r present'
c            stop
c         endif
c      endif

      rewind(unit=i_unit)
      read(i_unit,genr,iostat=kode)
      write(*,*)'in read_all_namelists after read genr kode=',kode
      call check_read(kode,'genr')
      write(*,*)'print data from /genr/ obtained from input'

cSAP091030
      if (partner.ne.'disabled') ionetwo=1
c----------------------------------------------------------------------
      call bcast(r_wall(1),0.d0,n_wall_a)
      call bcast(z_wall(1),0.d0,n_wall_a)
      call ibcast(n_limiter,0,max_limiters)
      call bcast(r_limiter(1,1),0.d0,n_limiter_a*max_limiters_a)
      call bcast(z_limiter(1,1),0.d0,n_limiter_a*max_limiters_a)
      
      rewind(unit=i_unit)
      read(i_unit,tokamak,iostat=kode)
      write(*,*)' prepare_genray_input after read tokamak kode=',kode
      write(*,*)'Tokamak data'
c      write(*,tokamak)   
      call check_read(kode,'tokamak')

c-----check the input data in namelist /tokamak/
      if ((indexrho.lt.1).or.(indexrho.gt.6)) then
         write(*,*)'prepare_genray_input in reading namelist /tokamak/'
         write(*,*)'It should be 0<indexrho<7, but indexrho =',indexrho
         write(*,*)'Change indexrho in genray.in file'
         stop
      endif

      if ((ipsi.lt.0).or.(ipsi.gt.1)) then
         write(*,*)'prepare_genray_inputi n reading namelist /tokamak/'
         write(*,*)'It should be -1<ipsi<2, but ipsi =',ipsi
         write(*,*)'Change ipsi in genray.in file'
         stop
      endif

      if ((ionetwo.lt.0).or.(ionetwo.gt.1)) then
         write(*,*)'prepare_genray_input in reading namelist /tokamak/'
         write(*,*)'It should be -1<ionetwo<2, but ionetwo =',ionetwo
         write(*,*)'Change ionetwo in genray.in file'
         stop
      endif

      if ((ieffic.lt.1).or.(ieffic.gt.6)) then
         write(*,*)'prepare_genray_input in reading namelist /tokamak/'
         write(*,*)'It should be 1 =<ieffic =< 6, but ieffic =',ieffic
         write(*,*)'Change ieffic in genray.in or genray.dat file'
         stop
      endif
      
      if ((psifactr.le.0).or.(psifactr.gt.1)) then
         write(*,*)'prepare_genray_input in reading namelitst /tokamak/'
         write(*,*)'It should be 0<psifactr=<1, but psifactr=',psifactr
         write(*,*)'Change psifactr in genray.in file'
         stop
      endif

      write(*,*)'deltripl,nloop,i_ripple',deltripl,nloop,i_ripple
      if((i_ripple.ne.1).and.(i_ripple.ne.2))then
        write(*,*)'prepare_genray_input i_riplle should be equal 1 or 2'
        write(*,*)'i_ripple=',i_ripple
        write(*,*)'Please change i_ripple in /tokamak/'
        write(*,*)'in genray.in.file'
        stop
      endif

      if (length_char(eqdskin).gt.512) then
         write(*,1001)
 1001    format('STOP: eqdskin spec too long')
         STOP
      endif
    
      if (NR.gt.NRA) then
        write(*,*)'NR > NRA'
        write(*,*)'it should be NR.le.NRA'
        write(*,*)'please change NR in genray.dat or NRA in param.i' 
        STOP 
      endif 

      if (n_wall.gt.n_wall_a) then
        write(*,*)'n_wall > n_wall_a'
        write(*,*)'it should be n_wall.le.n_wall_a'
        write(*,*)'Please change n_wall in genray.dat or
     &  n_wall_a in param.i' 
        STOP 
      endif

      if (max_limiters.gt.max_limiters_a) then
        write(*,*)'max_limiters > max_limiters_a'
        write(*,*)'it should be max_limiters.le.max_limiters_a'
        write(*,*)'Please change max_limiter_a in param.i' 
        write(*,*)'and recomple the code' 
        STOP 
      endif

      do i=1,max_limiters
         if (n_limiter(i).gt.n_limiter_a) then
           write(*,*)'n_limiter(i).gt.n_limiter_a) at i=',i
           write(*,*)'n_limiter(i),n_limiter_a',n_limiter(i),n_limiter_a
           write(*,*)'Please change n_limiter_a in param.i'
           write(*,*)'and recompile the code'
           STOP
         endif
      enddo

      if (n_wall.gt.0) then
         if((r_wall(1).ne.r_wall(n_wall)).or.
     &      (z_wall(1).ne.z_wall(n_wall))) then
            write(*,*)'The first wall point should coinside
     &                 with n_wall wall point'
            write(*,*)'But they do not coinside'
            write(*,*)'r_wall(1)=', r_wall(1)
            write(*,*)'r_wall(n_wall)=', r_wall(n_wall)
            write(*,*)'z_wall(1)=', z_wall(1)
            write(*,*)'z_wall(n_wall)=', z_wall(n_wall)
            write(*,*)'Please correct wall cordinates in genray.dat'
            write(*,*)'or genray.in file'
        endif
      endif
 
      h_add_wall=h_add_wall/r0x 
      if (n_wall.gt.0) then
         do i=1,n_wall
            z_wall(i)=z_wall(i)/r0x
            r_wall(i)=r_wall(i)/r0x
         enddo
      endif

      if (max_limiters.gt.0) then
        do i=1,max_limiters
           do j=1,n_limiter(i) 
              r_limiter(j,i)=r_limiter(j,i)/r0x 
              z_limiter(j,i)=z_limiter(j,i)/r0x
           enddo
        enddo     
      endif
      
c------------------------------------------------------------------------
      rewind(unit=i_unit)
      read(i_unit,wave,iostat=kode)
      call check_read(kode,'wave')
     
      if(((istart.eq.1).or.(istart.eq.3)).and.(ibw.eq.1)) then
        write(*,*)'ibw=1 but istart.ne.2'
        write(*,*)'Please change the input data in genray.in file'
        stop
      endif

c--------------------------------------------------------
c     i_vgr_ini determines the direction of the wave in the initial point
c     i_vg_ini=+1 wave group velocity is directed into the plasma
c             =-1 ouside the plasma
c     It was proposed that the poloidal flux has a minimum at the plasma
C        write(*,*)'i_vgr_inname_uniform_mesh_profiles.i'     i=',i_vgr_ini
C------------------------------------------------------------------------
      if (cnperp_plot_max.gt.cN_perp_root_max) then
        write(*,*)'================WARNING============'
        write(*,*)'cnperp_plot_max.gt.cN_perp_root_max'
        write(*,*)'cnperp_plot_max,cN_perp_root_max',
     &             cnperp_plot_max,cN_perp_root_max
        write(*,*)'The code will set cN_perp_root_max=cnperp_plot_max'
        cN_perp_root_max=cnperp_plot_max
        write(*,*)'==================================='
      endif

c---------------------------------------------------------
c     The data for  the n_perpendicular angle scattering
      rewind(unit=i_unit)
      read(i_unit,scatnper,iostat=kode)
      call check_read(kode,'scatnper')
c---------------------------------------------------------

      if (istart.eq.1) then
c---------------------------------------------------------
c        start point is outside the plasma, ECR case
c        the reading of the data for EC cone
c----------------------------------------------------------        
         rewind(unit=i_unit)             
         read(i_unit,eccone,iostat=kode)
         call check_read(kode,'eccone')

        if (raypatt.eq.'toray') then
      
           if (gzone.gt.gzonemax) then
              write(*,*)'number of ray pattern zones, gzone=',gzone
              write(*,*)'max value of gzone: gzonemax=',gzonemax
              write(*,*) 'gzone.gt.gzonemax' 
              write(*,*)'it should b egzone.le.gzonemax'
              write(*,*)'please decrease gzone in genray.dat'
              write(*,*)'or increase gzonemax in param.i and recomplie'
              stop! 'in prepare_genray_input.f gzone.gt.gzonemax'
           endif

           if (gzone.eq.0)then
              write(*,*)'WARNING!!!  nray_in must be set for gzone=0'
              write(*,*)'nray_in=',nray_in
           endif
          
        endif ! raypatt.eq.'toray'
   
        if ((raypatt.eq.'diskdisk').or.(raypatt.eq.'diskbeam')) then
c --------------------------------------------------------------
c          calculate the radius of the first disk using
c          sigma_launching_disk,part_gauss_power
c-----------------------------------------------------------------
           if( part_gauss_power.le.0.d0) then
              write(*,*)'part_gauss_power.le.0.d0 in input file'
              write(*,*)'part_gauss_power should be positive'
              write(*,*)'Please change part_gauss_power'
              write(*,*)'in input files genray.dat or genray.in'
              stop 'part_gauss_power'
           endif

           if (part_gauss_power.lt.1.d0) then
              sigma_launching_disk=rho_launching_disk/
     &        dsqrt(dlog(1.d0/(1.d0-part_gauss_power)))                  
           else
              part_gauss_power=1.d0-
     &        dexp(-(rho_launching_disk/sigma_launching_disk)**2)             
           endif
        
           ncone=1
           write(*,*)'ncone set =1 in read_write_genray_input eccone'

           nray=0
           do j=1,n_mesh_disk_radial_bin
             do i=1,n_mesh_disk_angle_bin(j)
                nray=nray+1
             enddo
           enddo
c----------nray is a number of all rays.It is the output argument.
       
           if (nray.gt.nraymax) then 
             write(*,*)'number of rays launched from the disk'
             write(*,*)'nray.gt.nraymax'
             write(*,*)'nray,nramax',nray,nraymax
             write(*,*)'Please increase nraymax in param.i'
             write(*,*)'and recompile the code'
             stop 'in prepare_genray_input.f'
           endif
         endif !raypatt = 'diskdisk'.or.'diskbeam'
    

      else 
         !istart= 2 or 3     
c-------------------------------------------------
c        start point is inside the plasma, LH or FW case
c        the reading of the data for LH grill
c--------------------------------------------------
c         call inigrill  
c         call check_param(2)
c         if (ngrilld.ne.0) then
c          pause  'Attention !!! genray.in file contains old ngrilld'
c 1000     format('Attention!!! genray.in file contains old ngrilld',/,
c     &   'It can be if the old genray.in version is used.',/,
c     &   'The new genray.in file uses ngrill instead of ngrilld',/,
c     &   'The code will put ngrill=ngrilld')
c          write(*,1000)
c          ngrill=ngrilld
c      endif !istart
    
         rewind(unit=i_unit)
         read(i_unit,grill,iostat=kode)
         call check_read(kode,'grill')
         write(*,grill)
              
         if (ngrilld.ne.0) then
         pause'Attention!!!genray.in file contains old variable ngrilld'
 1000       format('Attention!!! genray.in file contains old ngrilld',/,
     &      'It can be if the old genray.in version is used.',/,
     &      'The new genray.in file ues ngrill instead of ngrilld',/,
     &      'The code will put ngrill=ngrilld')
            write(*,1000)
            ngrill=ngrilld
         endif

         if (ilaunch.eq.1) then
           write(*,*)
           write(*,*)'grill_lh'
           write(*,*)' ilaunch=1 launches one ray in plasma from'
           write(*,*)'       the R0launch,Phi0launch,Z0Launch location.'
           write(*,*)'       Setting ngrill=1,nthin(1)=1,nnkpar(1)=1,'
           write(*,*)'       Setting nnktor(1)=1,nnkpol(1)=1'
           write(*,*)
           ngrill=1
           nthin(1)=1
           nnkpar(1)=1
           nnktor(1)=1
           nnkpol(1)=1
         endif

cSAP111028
         if (ilaunch.eq.2) then
           write(*,*)           
           do i=1,ngrill
              if(nthin(i).gt.nthinmax) then
                write(*,*)'grill_lh ilaunch=2'
                write(*,*)'i=','nthin(i).gt.nthinmax'
                write(*,*)'It should be nthin(i).le.nthinmax'
                write(*,*)'nthin(i),nthinmax',nthin(i),nthinmax
                write(*,*)'Please change nthinmax in param.i file'
                write(*,*)'and recompile the code'
                stop 'in read_write_genray_input'
              endif
           enddo
  
           call read_rz_launch_grill_nml(i_unit,kode)

         endif

         if ((i_grill_pol_mesh.ne.1).and.(i_grill_pol_mesh.ne.2)) then  
           write(*,*)'read_write genray_input.f in inigrill:'
           write(*,*)'     it should be i_grill_pol_mesh=1 or =2 but'
           write(*,*)'     i_grill_pol_mesh=',i_grill_pol_mesh
           write(*,*)'     Please change i_grill_pol_mesh in input file'
           write(*,*)'     genray.dat or genray.ini'
           stop 'after read grill'
         endif                

         write(*,*)'from grill i_n_poloidal,ngrill',i_n_poloidal,ngrill  
         if ((i_n_poloidal.lt.1).or.(i_n_poloidal.gt.4)) then
            write(*,*)'in inigrill i_n_poloidal<1 or i_n_poloidal>4' 
            write(*,*)'please change i_n_poloidal in genray.in'
            stop
         endif

         if((i_rho_cutoff.lt.0).or.(i_rho_cutoff.gt.1)) then
           write(*,*)'i_rho_cutoff= ',i_rho_cutoff
           write(*,*)'i_rho_cutoff<0 or i_rho_cutoff>1'
           write(*,*)'it should be i_rho_cutoff=0 or =1'
           write(*,*)'psease change i_rho_cutoff in genray.in file'
         endif
    
         if(ngrill.gt.ngrilla) then
           write(*,*)'ngrill.gt.ngrilla'
           write(*,*)'ngrill=',ngrill,'ngrilla=',ngrilla
           write(*,*)'please change these parameters'
           write(*,*)'in genray.in or in param.i'
           stop
         endif        
c-------------------------------
        powtott=0.d0
        do i=1,ngrill
c----------powers(1:ngrilla)  power in one grill (MWatts)
c          (total input power to grills(in MWatts) will 
c           be powtott=sum{powers})           
           powtott=powtott+powers(i)
        enddo
        write(*,*)'powtott=sum{input powers}=',powtott
c-------------------------------
c       control that parameters nnkprmax, nraymax,ntinmax in param.i
c       are suitable
c-------------------------------
c       nnkprmax must be equal max{i=1,ngrill}nnkpar(i)
c       nnkprmxn=max{i=1,ngrill}nnkpar(i)
        nnkprmxn=0
        do i=1,ngrill
          if(nnkpar(i).gt.nnkprmxn) then
            nnkprmxn=nnkpar(i)
          endif
        enddo     
        if(nnkprmxn.gt.nnkprmax) then
          write(*,*)'nnkprmxn=',nnkprmxn,'nnkprmax=',nnkprmax
          write(*,*)'nnkprmxn.gt.nnkprmax'
          write(*,*)'please change nnkprmax in param.i'
          stop
        endif
c-------------------------------
c       nthinmax must be equal max{i=1,ngrill}nthin(i)
c       nthinmxn=max{i=1,ngrill}nthin(i)
        nthinmxn=0
        do i=1,ngrill
          if(nthin(i).gt.nthinmxn) then
            nthinmxn=nthin(i)
          endif
        enddo
        if(nthinmxn.gt.nthinmax) then
          write(*,*)'nthinmxn=',nthinmxn,'nthinmax=',nthinmax
          write(*,*)'nthinmxn.gt.nthinmax'
          write(*,*)'please change nthinmax in param.i'
          stop
        endif        
c------------------------------------------------
c     end of the reading of the data for LH grill
c------------------------------------------------
      endif !istart
c--------------------------------------------------------------------
   
      rewind(unit=i_unit)
      read(i_unit,dispers,iostat=kode)
      call check_read(kode,'dispers')
          
      if(n_relt_harm.lt.1)then
        write(*,*)'n_relt_harm.lt.1'
        write(*,*)'n_relt_harm=',n_relt_harm
        write(*,*)'it should be n_relt_harm.ge.1'
        write(*,*)'Please change n_relt_harm in genray.dat'
        stop 'in read_write_genray_input.f check /dispers/'
      endif

      
      do i=1,nbulk 
        if ((i_salphal(i).ne.1).and.(i_salphal(i).ne.0)) then
          write(*,*)'(i_salphal(i).ne.1).or.(i_salphal(i).ne.0)'
          write(*,*)'It should i_salphal(i) =0 or =1'
          write(*,*)'i=',i,'i_salphal(i)',i_salphal(i)
          write(*,*)'Please chagne i_saplhal(i)'
          write(*,*)'in input genray.in or genray.dat file'
          stop 'reaf_write_genray_input.f i_salphal problem'
        endif
      enddo 

      if ((ion_absorption.ne.'enabled').and.
     &    (ion_absorption.ne.'disabled')) then
           write(*,*)'It should be ion_absorption=enabled'
           write(*,*)'or .and.ion_absorption=disabled'
           write(*,*)'ion_absorption=',ion_absorption
           write(*,*)'Please change ion_absorption'
           write(*,*)'in input genray.in or genray.dat file'
           stop 'reaf_write_genray_input.f ion_absorption problem'
       endif
c----------------------------------------------------------------------
      rewind(unit=i_unit)
      read(i_unit,numercl,iostat=kode)
      call check_read(kode,'numercl')
      
      write(*,*)'in read_write_genray_input.f after read(numercl)'
      write(*,numercl)

      if (nrelt.gt.nrelta) then
         write(*,*)'nrelt must be .le. nrelta (a parameter in param.i)'
         write(*,*)'nrel,nrelta',nrelt,nrelta
         stop
      endif

      if((i_resonance_curve_integration_method.lt.1).or.
     &   (i_resonance_curve_integration_method.gt.4)) then
         write(*,*)'it should be' 
         write(*,*)'i_resonance_curve_integration_method =1,2,3,4'
         write(*,*)'check in genray.dat '
         write(*,*)'i_resonance_curve_integration_method =',
     &   i_resonance_curve_integration_method
         stop
      endif

      ndim=ndim1
c -------------------------------------------------
c     'Runge-Kutta method parameters '
c-------------------------------------------------
      prmt(1)=prmt1   ! tau initial= prmt(1)
      prmt(2)=prmt2   ! tau final=prmt(2)
      prmt(3)=prmt3   ! initial tau step=prmt(3)
      prmt(4)=prmt4   ! required accuracy=prmt(4)
      prmt6=prmt6/r0x ! normalization of the output step
      prmt(6)=prmt6   ! hprint=prmt(6) time period for results output
c      ihlf=prmt9
      prmt(9)=prmt9   ! accuracy of hamiltonian
      write(*,*)'in read_all prmt=',prmt
c-----check numerical or analytical differentiation
      
      if((id.eq.1).or.(id.eq.2).or.(id.eq.3).or.(id.eq.6)) then
        icheck=1 !possible to use the analytical differentiation
      else    
        icheck=2 !impossible to use the analytical differentiation
      endif

      if ((idif.eq.1).and.(icheck.eq.2)) then
        write(*,*)'impossible to use analytical differentiation'
        write(*,*)'idif=1 for used dispersion function id= ',id
        write(*,*)'please set idif=2 in genray.dat or genray.in'
        stop 'read_write_genray_input.f idif=1' 
      endif

      if((i_output.ne.1).and.(i_output.ne.2)) then
         write(*,*)'i_output should be =1 or 2'
         write(*,*)'i_output=',i_output
         write(*,*)'please set i_output in genray.dat or genray.in'
         stop 'read_write_genray_input.f i_output'       
      endif
        
c------------------------------------------------------------------
      rewind(unit=i_unit)
      read(i_unit,output,iostat=kode)      
      call check_read(kode,'output')

      if((n_plot_disp.lt.0).or.(n_plot_disp.gt.n_plot_dispa)) then
        write(*,*)'in  prepare_genray_input'
        write(*,*)'(n_plot_disp.lt.0).or.(n_plot_disp.gt.n_plot_dispa)'
        write(*,*)'n_plot_disp,n_plot_dispa',n_plot_disp,n_plot_dispa
        write(*,*)'please change n_plot_disp in genray.dat or'
        write(*,*)'increase n_plot_dispa in param.i'
        write(*,*)'and recomplile the code'
        stop 'output'
      endif

      if((n_plot_disp_cold.lt.0).or.
     &   (n_plot_disp_cold.gt.n_plot_disp_colda)) then
        write(*,*)'in in  prepare_genray_inputd'
        write(*,*)'(n_plot_disp_cold.lt.0).or.'
        write(*,*)'(n_plot_disp_cold.gt.n_plot_disp_colda)'
        write(*,*)'n_plot_disp_cold,n_plot_disp_colda=',
     &  n_plot_disp_cold,n_plot_disp_colda
        write(*,*)'please change n_plot_disp_cold in genray.dat'
        write(*,*)'or increase n_plot_disp_colda in param.i'
        write(*,*)'and recomplile the code'
        stop 'read_write__genray_input.f output'
      endif
     
      if ((i_plot_wave_normal_cold.ne.1).and.
     &    (i_plot_wave_normal_cold.ne.0)) then
         write(*,*)'in  prepare_genray_input'
         write(*,*)'i_plot_wave_normal_cold .ne.1 and .ne.0'
         write(*,*)'i_plot_wave_normal_cold=',i_plot_wave_normal_cold
         write(*,*)'Please change i_plot_wave_normal_cold in genray.dat'
         stop 'prepare_genray_input.f'
      endif

      if (n_contour_plot_disp.gt.n_contour_plot_disp_a)then
         write(*,*)'in  prepare_genray_input'
         write(*,*)'n_contour_plot_disp.gt.n_contour_plot_disp_a'
         write(*,*)'n_contour_plot_disp',n_contour_plot_disp
         write(*,*)'n_contour_plot_disp_a',n_contour_plot_disp_a
         write(*,*)'Please change n_contour_plot_disp_a in param.i'
         write(*,*)'and recompile the code'
         write(*,*)'or change n_contour_plot_disp in genray.dat'
         stop 'prepare_genray_input.f'
      endif
c---------------------------------------------------------
      rewind(unit=i_unit)
      read(i_unit,plasma,iostat=kode)
      call check_read(kode,'plasma')
      
c     nbulk>=1 number of plasma components
      if (nbulk.gt.nbulka) then
         write(*,*)'nbulka=',nbulka
         write(*,*)'nbulk=',nbulk
	 write(*,*)'nbulka.lt.nbulk'
	 write(*,*)'change the parameter nbulka in param.i'
	 stop
      endif
     
      if (ndens.gt.ndensa) then
         write(*,*)'ndens=',ndens,'ndensa=',ndensa
	 write(*,*)'ndensa.lt.ndens'
         write(*,*)'ndensa is given in param.i'
         write(*,*)'ndens  is given in genray.dat'
	 write(*,*)'change the paramter ndensa in param.i'         
	 stop
      endif

      if((izeff.eq.3).and.(nbulk.eq.1)) then
        write(*,*)'for izeff=1 it is necessary that  nbulk > 1'
        write(*,*)'izeff=',izeff,'nbulk=',nbulk
        write(*,*)'change izeff or nbulk and profiles in genray.in'
        stop 
      endif

c-----------------------------------------------------
      rewind(unit=i_unit)
      read(i_unit,species,iostat=kode)
      call check_read(kode,'species')
     
c   plasma component charges charge(i)=mod(charge(i)/charge_electron)
c----------------------------------------------------
c     electron charge charge(1) should be =1
      if( charge(1).ne.1.d0) then 
        write(*,*)'Warning in dinit: charge(1) should be equal=1'
        write(*,*)'but charge(1)=',charge(1),'control it'
	stop
      endif
      do i=2,nbulk
         if(charge(i).lt.charge(i-1). and. izeff.ne.2) then
      write(*,*)'Warning in dinit:it should be charge(i).ge.charge(i-1)'
      write(*,*)'But in i=',i,'charge(i).lt.charge(i-1)'
           write(*,*)'Please correct genray.in file'
	   stop
	 endif
      enddo

c$$$      do i=1,nbulk
c$$$         write(*,*)'i, charge(i)',i,charge(i)
c$$$      enddo
c-----------------------------------------------------
c     plasma components mass dmas(i)=Mass(i)/Mass_electron
c-----------------------------------------------------
c$$$      do i=1,nbulk
c$$$         write(*,*)'i, dmas(i)',i,dmas(i)
c$$$      enddo
c     calculation of nbulk1
      if(((izeff.eq.0).or.(izeff.eq.2)).or.(izeff.eq.3)) then
c        izeff=0, zeff will be calculated using the given ions;
c                 electron density will be calculated using ion's densities;
c             =1  ion densities nbulk and nbulk-1 will be calculated
c                 in dinit.f  using
c                 Zeff, electon density and ion's densities(i), i=2,nbulk-1;
c        izeff=2, zeff will not coincide with the plasma components
c             =3  it uses eqdsk pres (pressure) and ions densities_i
c                 for i=2,... nbulk
c                 Let temperature T_E=T_i
c                 pres=dens1(k,1)*temp1(k,1)+
c                      Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
c                 In this case we will calculate Zeff(rho),
c                 dens_electron(rho) and T_e(rho)=T_i(rho)
c             =4  it uses eqdsk pres (pressure), zeff,ions densities
c                 for i=2,... nbulk-2 (nbulk>2) and T_e(rho),T_i(rho)
c                 pres=dens1(k,1)*temp1(k,1)+
c                      Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
c                 In this case we will calculate dens_electron(rho) and
c                 ion densities for i=nbulk and i=nbulk-1)
         nbulk1=nbulk         
      else
c        izeff=1, zeff is given, the ions component will be calculated
         if (nbulk.eq.1) nbulk1=1
         if (nbulk.eq.2) then
	    nbulk1=2
	    write(*,*)'nbulk=2 Zeff must be equal charge(2) control it'
	    write(*,*)'use the option izeff=0'
	    stop
	 endif
         if (nbulk.gt.2) nbulk1=nbulk-2
      endif !izeff
      write(*,*)'nbulk1=',nbulk1
c------------------------------------------------------------------
c     small radius uniform mesh for plasma profiles
c     It will be recalculated in dinit.f for common  six_no_nml.i
c-----------------------------------------------------
      h=1.d0/(ndens-1)
      do i=1,ndens
        rhom(i)=h*(i-1)
      enddo
c------------------------------------------------------------------
c     The parameters for the density fluctuations
      rewind(unit=i_unit)
      read(i_unit,varden,iostat=kode)
      call check_read(kode,'varden')  
      
c------------------------------------------------------------------
c      call bcast(dens1(1,1),0.d0,nbulka*ndensa)
c      call bcast(temp1(1,1),0.d0,nbulka*ndensa)
c      call bcast(tpop1(1,1),1.d0,nbulka*ndensa)
c      call bcast(vflow1(1,1),0.d0,nbulka*ndensa)

      if(idens.eq.0) then
c--------Analytical radial profiles        
         rewind(unit=i_unit)
	 read(i_unit,denprof,iostat=kode)
         call check_read(kode,'denprof')
  
         rewind(unit=i_unit)
         read(i_unit,tpopprof,iostat=kode)
         call check_read(kode,'tpopprof')

         rewind(unit=i_unit)
         read(i_unit,vflprof,iostat=kode)
         call check_read(kode,'vfloprof')

         rewind(unit=i_unit)
	 read(i_unit,zprof,iostat=kode)
         call check_read(kode,'zprof')

         rewind(unit=i_unit)
	 read(i_unit,tprof,iostat=kode)
         call check_read(kode,'tprof')

      endif ! idens analytical

      write(*,*)'prepare_genray_input: partner=',partner



c     Following if(partner....)-clause to be commented out giving
c     read_write_genray_input_prep.f, for use with 
c     prepare_genray_input.f90.
      if (partner.eq.'genray_profs_in.txt' .or. 
     1    partner.eq.'genray_profs_in.nc') then
c-----------------------------------------------------------
c        read density ,temperature and zeff  profiles from  
c        the file created by ONETWO : genray_profs_in
c
c        Read in from external file
c        1) the density, temperature and zeff profiles
c          at genray small radial mesh:
c          dens1,temp1,zeff1,eqdskin
c        2) nbulk - the number of plasma species 
c        3) set izeff =2 (for nbulk.ge.1)
c-----------------------------------------------------------
         call read_transport_prof(ndens,charge,dmas,
     &   nbulk,dens1,temp1,tpop1,vflow1,zeff1,izeff,eqdskin,partner)

cSAP080711
         if (genray_in_dat.eq.'genray.in') then
c-----------All data should be for this case in MKSA (genray.in form). 
c           However read_transport_prof gives out dens1 in genray.dat form.
c           Transform dens1 from genray.dat to MKSA form.
            do i=1,nbulk
	      do k=1,ndens                            
                 dens1(k,i)=dens1(k,i)*1.d19
              enddo
            enddo
         endif ! genray_in_dat='genray.in' 

         go to 20    
      endif         ! partner

c----------------------------------
      if (idens.eq.1) then
c-----------------------------------------------------------------------
c        spline approximation of the density, temperature, zeff,
c        tpop=T_perp/T_parallel, vflow
c        radial profiles
c---------------------------------------------------------------------	 
         if(ndens.gt.ndensa) then
            write(*,*)'ndensa,ndens',ndensa,ndens
            write(*,*)'ndens > ndensa'
            Write(*,*)'it should be ndens.le.ndensa'
            stop 'read_write_genray_input.f'
         endif

         if(nbulk.gt.nbulka) then
            write(*,*)'nbulka,nbulk',nbulka,nbulk
            write(*,*)'nbulk > nbulka'
            Write(*,*)'it should be nbulk.le.nbulka'
            stop 'read_write_genray_input.f'
         endif

c=====================================================================
c        read density profiles  
c------------------------------------------------------------------ 
         if(nonuniform_profile_mesh.eq.'disabled') then 
c--------------------------------------------------------------------
c          read density profiles 'dentab' at
c          uniform mesh in column form
c------------------------------------------------------------------    
           call read_uniform_column_profile(i_unit,
     &     'dentab',nbulk,
     &      ndens,prof2_uniform,prof1_uniform,kode)
         else 
c--------------------------------------------------------------------
c          read density 'profiles dentab_nonuniform_line' at
c          nonuniform mesh in row form
c------------------------------------------------------------------
           call read_nonuniform_line_profile (i_unit,
     &     'dentab_nonuniform_line',nbulk,ndens,
     &     prof2_uniform,
     &     dens1_nonuniform,radii_nonuniform_dens1,nj_tab_dens1,kode)

         endif
c-------------------------------------------------------------------        
	 if ((izeff.eq.0).or.(izeff.eq.3)) then
	   i1=2
	 else
	   i1=1
	 endif

cSAP090315
         write(*,*)'at density reading kode=',kode
         if(kode.eq.0) then
c----------density data reading has complited succefully
           do i=i1,nbulk1
	    do k=1,ndens
               dens1(k,i)=prof2_uniform(k,i)
            enddo
           enddo
         else
           if (kode.lt.0) then             
              write(*,*)'end of input file was detected at reading'
              write(*,*)'dentab or dentab_nonuniform_line'
              write(*,*)'data will be set in default_in'
           else 
              write(*,*)'an error has occurred at reading'
              write(*,*)'dentab or dentab_nonuniform_line' 
              write(*,*)'Please change the namelist for density' 
              stop 'dentab or dentab_nonuniform_line'
           endif
         endif

21       format(5e16.9)

c        end the density profiles reading
c=====================================================================
c        read temperature profiles  
c------------------------------------------------------------------ 
         if(nonuniform_profile_mesh.eq.'disabled') then 
c--------------------------------------------------------------------
c          read temperature profiles 'temtab' at
c          uniform mesh in column form
c------------------------------------------------------------------    
           call read_uniform_column_profile(i_unit,
     &     'temtab',nbulk,
     &     ndens,prof2_uniform,prof1_uniform,kode)
         else 
c--------------------------------------------------------------------
c          read temperature profiles 'temtab_nonuniform_line' at
c          nonuniform mesh in row form
c------------------------------------------------------------------
           call read_nonuniform_line_profile(i_unit,
     &     'temtab_nonuniform_line',nbulk,ndens,
     &     prof2_uniform,
     &     temp1_nonuniform,radii_nonuniform_temp1,nj_tab_temp1,kode)           
         endif
c------------------------------------------------------------------     
cSAP090315
         write(*,*)'at temperature reading kode=',kode
         if(kode.eq.0) then
c----------temperature data reading has complited succefully
           do i=1,nbulk
	    do k=1,ndens
               temp1(k,i)=prof2_uniform(k,i)
            enddo
           enddo
         else
           if (kode.lt.0) then             
              write(*,*)'end of input file was detected at reading'
              write(*,*)'temtab or temtab_nonuniform_line' 
              write(*,*)'data will be set in default_in'
           else 
              write(*,*)'an error has occurred at reading'
              write(*,*)'temtab or temtab_nonuniform_line' 
              write(*,*)'Please change the namelist for temperature' 
              stop 'temtab or temtab_nonuniform_line'
           endif
         endif  
         

c        end temperature reading
c=====================================================================
c        read tpop profiles  
c------------------------------------------------------------------ 
         if(nonuniform_profile_mesh.eq.'disabled') then 
c--------------------------------------------------------------------
c          read tpop profiles 'tpotab' at
c          uniform mesh in column form
c------------------------------------------------------------------    
           call read_uniform_column_profile(i_unit,
     &     'tpoptab',nbulk,
     &     ndens,prof2_uniform,prof1_uniform,kode)
         else 
c--------------------------------------------------------------------
c          read tpop profiles 'tpoptab_nonuniform_line' at
c          nonuniform mesh in row form
c------------------------------------------------------------------
cSAP090311
           write(*,*)'before call read_nonuniform_line_profile'
           write(*,*)'tpoptab_nonuniform_line'

           call read_nonuniform_line_profile(i_unit,
     &     'tpoptab_nonuniform_line',nbulk,ndens,
     &     prof2_uniform,
     &     tpop1_nonuniform,radii_nonuniform_tpop1,nj_tab_tpop1,kode)

cSAP090311
         write(*,*)'after call read_nonuniform_line_profile'
         write(*,*)'tpoptab_nonuniform_line'
         write(*,*)'nj_tab_tpop1',nj_tab_tpop1

         endif

cSAP090315
         write(*,*)'at tpop reading kode=',kode
         if(kode.eq.0) then
c----------tpop data reading has complited succefully
           do i=1,nbulk
	    do k=1,ndens
               tpop1(k,i)=prof2_uniform(k,i)
            enddo
           enddo
         else
           if (kode.lt.0) then             
              write(*,*)'end of input file was detected at reading'
              write(*,*)'tpoptab or tpoptab_nonuniform_line' 
              write(*,*)'data will be set in default_in'

              do i=1,nbulk
	        do k=1,ndens
                   write(*,*)'k,i,tpop1(k,i)',k,i,tpop1(k,i)
                enddo
              enddo

           else 
              write(*,*)'an error has occurred at reading'
              write(*,*)'tpoptab or tpoptab_nonuniform_line'
              write(*,*)'Please change the namelist for tpop' 
              stop 'tpoptab or tpoptab_nonuniform_line' 
           endif
         endif  

c------------------------------------------------------------------
c        end tpop reading
c=====================================================================
c        read vflow profiles  
c------------------------------------------------------------------ 
         if(nonuniform_profile_mesh.eq.'disabled') then 
c--------------------------------------------------------------------
c          read tpop profiles 'vflowtab' at
c          uniform mesh in column form
c------------------------------------------------------------------    
           call read_uniform_column_profile(i_unit,
     &     'vflowtab',nbulk,
     &     ndens,prof2_uniform,prof1_uniform,kode)
         else 
c--------------------------------------------------------------------
c          read tpop profiles 'vflowtab_nonuniform_line' at
c          nonuniform mesh in row form
c------------------------------------------------------------------
           call read_nonuniform_line_profile(i_unit,
     &     'vflowtab_nonuniform_line',nbulk,ndens,
     &     prof2_uniform,
     &     vflow1_nonuniform,radii_nonuniform_vflow1,nj_tab_vflow1,kode)
         endif
c------------------------------------------------------------------
cSAP090315
         write(*,*)'at vflow reading kode=',kode
         if(kode.eq.0) then
c----------vflow data reading has complited succefully
           do i=1,nbulk
	    do k=1,ndens
               vflow1(k,i)=prof2_uniform(k,i)
            enddo
           enddo
         else
           if (kode.lt.0) then                           
              write(*,*)'end of input file was detected at reading'
              write(*,*)'vflowtab or vflowtab_nonuniform_line' 
              write(*,*)'data will be set in default_in'
           else 
              write(*,*)'an error has occurred at reading'
              write(*,*)'vflowtab or vflowtab_nonuniform_line'
              write(*,*)'Please change the namelist for vflow' 
              stop 'vflowtab or vflow_nonuniform_line' 
           endif
         endif  
         

c        end vflow reading 
c=====================================================================
         
c        read zeff profiles  
c------------------------------------------------------------------ 
         if(nonuniform_profile_mesh.eq.'disabled') then 
c--------------------------------------------------------------------
c          read zeff profile 'zeftab' at
c          uniform mesh in column form
c------------------------------------------------------------------    
           call read_uniform_column_profile(i_unit,
     &     'zeftab',nbulk,
     &     ndens,prof2_uniform,zeff1,kode) 
           write(*,zeftab)
           write(*,*)'after read zeff1',zeff1         
         else 
c--------------------------------------------------------------------
c          read zeff profile 'zeftab_nonuniform_line' at
c          nonuniform mesh in row form
c------------------------------------------------------------------
cSAP090311
           write(*,*)'in read_all_namelists before'
           write(*,*)'read_nonuniform_line_profile'

           call read_nonuniform_line_profile(i_unit,
     &     'zeftab_nonuniform_line',nbulk,ndens,
     &     prof2_uniform,
     &     zeff1_nonuniform,radii_nonuniform_zeff1,nj_tab_zeff1,kode)

	   do k=1,ndens
             zeff1(k)=prof2_uniform(k,1)
             write(*,*)'k,zeff1(k)',k,zeff1(k)
           enddo
          
           
         endif

cSAP090315
         write(*,*)'at zeff reading kode=',kode
         if(kode.eq.0) then
c----------zeff data reading has complited succefully
           
         else
           if (kode.lt.0) then                           
              write(*,*)'end of input file was detected at reading'
              write(*,*)'zeftab or zeftab_nonuniform_line' 
              write(*,*)'data will be set in default_in'
           else 
              write(*,*)'an error has occurred at reading'
              write(*,*)'zeftab or zeftab_nonuniform_line'
              write(*,*)'Please change the namelist for zeff' 
              stop 'zeftab or zeftab_nonuniform_line'
           endif
         endif  
         
cSAP090311
           write(*,*)'in read_all_namelists after'
           write(*,*)'read_nonuniform_line_profile'
           write(*,*)'nj_tab_zeff1',nj_tab_zeff1
           write(*,*)'zeff1_nonuniform',zeff1_nonuniform
           write(*,*)'radii_nonuniform_zeff1',radii_nonuniform_zeff1

c------------------------------------------------------------------            
c        end zeff reading
c=================================================================
      endif! idens=1
c-----------------------------------------------------------------
 20   continue !from:  if (partner.eq.
c-----------------------------------------------------------
c     reads the data for reading 3D distribution function from diskf      
c     written by CQL3D
      rewind(unit=i_unit) 
      read(i_unit,read_diskf,iostat=kode)
      call check_read(kode,'read_diskf')
 
c------------------------------------------------------------
c     read the data for emission calculations
      rewind(unit=i_unit)
      read(i_unit,emission,iostat=kode)
      call check_read(kode,'emission')
      
      if (i_emission.eq.1) then

c       if (nfreqa.lt.nfreq) then
c        write(*,*)'dinit_mr nfreqa<nfreq'
c        write(*,*)'it should be nfreqa.ge.nfreq'
c        write(*,*)'nfreq,nfreqa',nfreq,nfreqa
c        write(*,*)'change nfreqa in param.i and recompile the code
c     +  or reduce nfreq in genray.in'
c        stop
c       endif


c      if ((i_emission_spectrum.eq.1).and.(jx_kin_a.lt.jx_kin)) then
c       write(*,*)'dinit_mr jx_kina<jx_kin'
c       write(*,*)'it should be jx_kin_a.ge.jx'
c       write(*,*)'change jx_kin_a in param.i and recompile the code
c    +  or reduce jx_kin in genray.in'
c       stop
c      endif
 
       if ((wallr.lt.0.d0).or.(wallr.gt.1.d0)) then
        write(*,*)'it should be {0=< wallr =<1}'
        write(*,*)'but wallr=',wallr
        write(*,*)'change wallr in genray.in'
        stop
       endif
      endif ! i_emission=1

c---------------------------------------------------------
c     read the data for for EC cone vertex coordinates calculation.
c     This case is used for the optimal OX mode conversion.
      rewind(unit=i_unit)
      read(i_unit,ox,iostat=kode)
      call check_read(kode,'ox')

      if(i_ox.eq.1) then
        istart=3
        prmt(3)=-prmt3 !to create the negative time
        i_vgr_ini=+1
        ireflm=1   
      endif

      if(((i_ox.ne.0).and.(i_ox.ne.1)).and.(i_ox.ne.2)) then
         write(*,*)'i_ox can  =0 or =1 or =2'
         write(*,*)'in namelist /ox/ i_ox=',i_ox
         write(*,*)'please change i_ox in input file'
         stop 'in prepare_genray_input.f  /ox/'
      endif
c---------------------------------------------------------
c     read the data for for adj calculations
    
      rewind(unit=i_unit)
      read(i_unit,adj_nml,iostat=kode)
      call check_read(kode,'adj_nml')      

c      if ( nmax_chi.gt.nmax_chi_a) then
c         write(*,*)'in read_write_genray_input.f'
c         write(*,*)'nmax_chi.gt.nmax_chi_a'
c         write(*,*)'it should be nmax_chi.le.nmax_chi_a'
c         write(*,*)'nmax_chi,nmax_chi_a',nmax_chi,nmax_chi_a
c         write(*,*)'Please increase nmax_chi_a in param.i'
c         write(*,*)'and recomplile the code'
c         stop 'in read_write_genray_input.f'
c      endif

c      if ( imax_chi.gt.nmax_chi_a) then
c         write(*,*)'in read_write_genray_input.f'
c         write(*,*)'imax_chi.gt.imax_chi_a'
c         write(*,*)'it should be imax_chi.le.imax_chi_a'
c         write(*,*)'imax_chi,imax_chi_a',imax_chi,imax_chi_a
c         write(*,*)'Please increase imax_chi_a in param.i'
c         write(*,*)'and recomplile the code'
c         stop 'in read_write_genray_input.f'
c      endif

c      if ( npsi0.gt.npsi0_a) then
c         write(*,*)'in read_write_genray_input.f'
c         write(*,*)'npsi0.gt.npsi0_a'
c         write(*,*)'it should be npsi0.le.npsi0_a'
c         write(*,*)'npsi0,npsi0_a',npsi0,npsi0_a
c         write(*,*)'Please increase ,npsi0_a in param.i'
c         write(*,*)'and recomplile the code'
c         stop 'in read_write_genray_input.f'
c      endif 

c      if ( nthp0.gt.nthp0_a) then
c         write(*,*)'in read_write_genray_input.f'
c         write(*,*)'nthp0.gt.nthp0_a'
c         write(*,*)'it should be nthp0.le.nthp0_a'
c         write(*,*)'nthp0,nthp0_a',nthp0,nthp0_a
c         write(*,*)'Please increase ,nthp0_a in param.i'
c         write(*,*)'and recomplile the code'
c         stop 'in read_write_genray_input.f'
c      endif

      if ((ieffic.eq.5).and.(i_adj.ne.1)) then
         write(*,*)'in read_write_genray_input.f'
         write(*,*)'in reading namelist /adj_nml/'
         write(*,*)'ieffic=5 but i_adj.ne.1'
         write(*,*)'It should be i_adj=1 for ieffic=5'
         write(*,*)'ieffic=',ieffic,'i_adj=',i_adj
         write(*,*)'Please change ieffic or i_adj'
         write(*,*)'in genray.in or genray.dat file'
         stop
      endif

      if ((ieffic.eq.6).and.(i_adj.ne.1)) then
         write(*,*)'in read_write_genray_input.f'
         write(*,*)'in reading namelist /adj_nml/'
         write(*,*)'ieffic=6 but i_adj.ne.1'
         write(*,*)'It should be i_adj=1 for ieffic=6'
         write(*,*)'ieffic=',ieffic,'i_adj=',i_adj
         write(*,*)'Please change ieffic or i_adj'
         write(*,*)'in genray.in or genray.dat file'
         stop
      endif

      if ((i_chi_interpolation.ne.1).and.
     &    (i_chi_interpolation.ne.2)) then 
         write(*,*)'in read_write_genray_input.f'
         write(*,*)'in reading namelist /adj_nml/'
         write(*,*)'it shoud be i_chi_interpolation=1 or =2'
         write(*,*)'but i_chi_interpolation=',i_chi_interpolation
         write(*,*)'Please change i_chi_interpolation'
         write(*,*)'in genray.in or genray.dat file'
         stop 'ir in reading namelist /adj_nml/'
      endif
c---------------------------------------------------------
cSAP090203
c     read the data for for calculations density profile
c     outside LCFS 
    
      rewind(unit=i_unit)
      write(*,*)'before edge_prof_nml'
      read(i_unit,edge_prof_nml,iostat=kode)
      call check_read(kode,'edge_prof_nml')
      write(*,*)'after  edge_prof_nml'
      write(*,edge_prof_nml)

      if ((i_edge_dens_anal.ne.0).and.
     &    ((i_edge_dens_anal.ne.1).and.(i_edge_dens_anal.ne.2))) then
         write(*,*)'in read_write_genray_input.f'
         write(*,*)'i_edge_dens_anal.ne.0 and .ne.0, 1 or2'
         write(*,*)'it should be i_edge_dens_anal =0 or =1 or =2 '
         write(*,*)'i_edge_dens_anal=',i_edge_dens_anal
         write(*,*)'Please change i_edge_dens_anal'
         write(*,*)'in genray.dat or in genray.in file'
         stop 'in read_write_genray_input.f'
      endif

      if (((i_edge_dens_rz_mesh.ne.0).and.
     &  (i_edge_dens_rz_mesh.ne.1)).and.(i_edge_dens_rz_mesh.ne.2)) then
         write(*,*)'in read_write_genray_input.f'
         write(*,*)'i_edge_dens_rz_mesh .ne.0 and .ne.1 and .ne.2'
         write(*,*)'it should be i_edge_dens_rz_mesh =0 or =1  or =2'
         write(*,*)'i_edge_dens_rz_mesh=',i_edge_dens_rz_mesh 
         write(*,*)'Please change i_edge_dens_rz_mesh'
         write(*,*)'in genray.dat or in genray.in file'
         stop 'in read_write_genray_input.f'
      endif

      if (n_pol_edge_dens.gt.n_pol_edge_dens_a) then
         write(*,*)'in read_write_genray_input.f'
         write(*,*)'n_pol_edge_dens.gt.n_pol_edge_dens_a'
         write(*,*)'it should be n_pol_edge_dens.le.n_pol_edge_dens_a'
         write(*,*)'n_pol_edge_dens=',n_pol_edge_dens
         write(*,*)'n_pol_edge_dens_a=',n_pol_edge_dens_a
         write(*,*)'Please increase n_pol_edge_dens_a in param.i'
         write(*,*)'and recomplile the code'
         stop 'in read_write_genray_input.f'
      endif

c      if (nxeqd_add.gt. nxeqd_add) then
c         write(*,*)'in read_write_genray_input.f'
c         write(*,*)'nxeqd_add.gt. nxeqd_add_a'
c         write(*,*)'it should be nxeqd_add.le.nxeqd_add_a'
c         write(*,*)'nxeqd_add=',nxeqd_add
c         write(*,*)'nxeqd_add_a=',nxeqd_add_a
c         write(*,*)'Please increase nxeqd_add_a in param.i'
c         write(*,*)'and recomplile the code' 
c         stop 'in read_write_genray_input.f'
c      endif

c      if (nyeqd_add.gt.nyeqd_add_a) then
c         write(*,*)'in read_write_genray_input.f'
c         write(*,*)'nyeqd_add.gt.nyeqd_add_a'
c         write(*,*)'it should be nyeqd_add.le.nyeqd_add_a'
c         write(*,*)'nyeqd_add=',nyeqd_add
c         write(*,*)'nyeqd_add_a=',nyeqd_add_a
c         write(*,*)'Please increase nyeqd_add_a in param.i'
c         write(*,*)'and recomplile the code' 
c         stop 'in read_write_genray_input.f'
c      endif


c-----test edge profiles spline
c-----calculte tables for edge_prof
c     theta_pol_edge_dens_ar_degree(i=1,n_pol_edge_dens)
c     sigmedgn_ar(i=1,n_pol_edge_dens)
c
c     from the input data for analytical profile
c     using analytical profile  like in
c     subroutine sigma_edge_n_theta_pol(theta_pol_radian
c
c     calculated profiles will be in edge_prof_nml.i

c      write(*,*)'WARNINIG the aed table was created using test tables'
c      write(*,*)'using subroutine create_edge_prof_table'
c      call create_edge_prof_table
c-----end    test edge profiles spline  
 
c      stop 'read_read_all_namelists '
c-----------------------------------------------------------
      rewind(unit=i_unit)
      read(i_unit,lsc_approach_nml,iostat=kode)
      call check_read(kode,'lsc_approach_nml')
      write(*,lsc_approach_nml)    
    
      if(i_lsc_approach.eq.1) then
         if (n_psi_TSC.gt.ndensa) then
            write(*,*)'n_psi_TSC.gt.ndensa'
            write(*,*)'n_psi_TSC=',n_psi_TSC
            write(*,*)'ndensa=',ndensa
            write(*,*)'it should be n_psi_TSC.le.ndensa'
            write(*,*)'please change n_psi_TS in genray.dat'
            write(*,*)'or change ndensda in param.i'           
            STOP 
         endif

         iabsorp=2 !to use part with iabsorp.eq.2 in prep3d.f
         write(*,*)'WARNNING: it was set iabsorp=2'

c--------read E_DC TSC radial profile at uniform mesh
         call read_uniform_column_profile(i_unit,'EdcTSCtab',1,
     &   n_psi_TSC,prof2_uniform,EdcTSC_1D,kode) 
         write(*,EdcTSCtab)
         write(*,*)'EdcTSC_1D',EdcTSC_1D
c--------read J parallel TSC radial profile at uniform mesh
c         call read_uniform_column_profile(i_unit,'JparTSCtab',1,
c     &   n_psi_TSC,prof2_uniform,JparTSC_1D,kode)
c         write(*,EdcTSCtab)
c         write(*,*)'Jpar_1D',JparTSC_1D
      endif   
 
c-------------------------------------------------------------
c     end of reading genray.dat or genray.in file
         
      close(i_unit)      
      write(*,*)'in prepare_genray_input.f end reading genray.in file'  

c     end of reading genray.in file

      return

      end     !!!!!!!!!!!!!!!!!read_all_namelists!!!!!!!!!!!!

c============================================================

      subroutine write_all_namelists(ndim)
!     [end at line ~1903]
c--------------------------------------------------------------
c     writes all namelists to genray.in or genray.dat input file
c--------------------------------------------------------------
      implicit none

      include 'param.i'
c---------------------------------------------------------
c-----input 
      integer ndim !the number of ray-tracing equations
c----------------------------------------
c     declaration and storage of all namelists variables
c-----------------------------------
      include 'cone_nml.i'
      include 'emissa_nml.i'
      include 'grill_nml.i'
      include 'ions_nml.i'
      include 'one_nml.i'
      include 'onetwo_nml.i'
      include 'output_nml.i'
      include 'scatnper_nml.i'    
      include 'six_nml.i'  !dens1(ndensa,nbulka),temp1(ndensa,nbulka)
                           !tpop1(ndensa,nbulka),vflow1(ndensa,nbulka)
                           !zeff1_nonuniform(ndensa,nbulka)
      include 'adj_nml.i'
cSAP090203
      include 'edge_prof_nml.i'
cSAP091222
      include 'lsc_approach_nml.i'
c--------------------------------------------------------------
c     declaration of namelist variables which
c     are used in subroutine dinit_mr and
c     are not used in any common block
c--------------------------------------------------------------  
      include 'dinit_nml.i' 
c--------------------------------------------------------------
      include 'rkutta.i' ! contains common /rkutta/ prmt(9)
c--------------------------------------------------------------
c     list of all namelists
c--------------------------------------------------------------
      include 'name.i'    ! list of all namelists

c      integer nray
c-----local
      integer kode,icheck,i,k,j,nbulk1,i1,i_unit,
     &i_genray_in_transformation
      real*8 powtott, nnkprmxn, nthinmxn, psi0, rho0,
     &rhom(ndensa),h,rho,
     &te0(nbulka),teb(nbulka)
c-----external
      integer length_char
      real*8 fpsi, rhopsi
    
      i_unit=1
cSAP100914
c      i_genray_in_transformation=1    ! to rewrite data into genray.in 
c                                      ! in MKSA system
c      open(i_unit,file='genray.in',delim='apostrophe',
c     &     status='old',iostat=kode)

      i_genray_in_transformation=0    ! do not rewrite data into genray.in 
                                      ! in MKSA system
      open(i_unit,file='genray.dat',delim='apostrophe',
     &     status='old',iostat=kode)

      if (kode.ne.0) then
cSAP100914
c         open(i_unit,file='genray.dat',delim='apostrophe',
c     &        status='old',iostat=kode)
c         i_genray_in_transformation=0 ! to rewrite data into genray.dat
         open(i_unit,file='genray.in',delim='apostrophe',
     &        status='old',iostat=kode)
         i_genray_in_transformation=1! to rewrite data into genray.in

         if (kode.ne.0) then
            write(*,*)' prepare_genray_input:'
            write(*,*)' Neither genray.in or genray.dat r present'
            stop
         endif
      endif

      write(i_unit,genr,iostat=kode)     
      call check_read(kode,'genr')
c----------------------------------------------------------------------

      write(i_unit,tokamak,iostat=kode)     
      call check_read(kode,'tokamak')

      if (n_wall.gt.0) then !normalizatiom
         do i=1,n_wall
            z_wall(i)=z_wall(i)*r0x
            r_wall(i)=r_wall(i)*r0x
         enddo
      endif

      if (max_limiters.gt.0) then
        do i=1,max_limiters
           do j=1,n_limiter(i) 
              r_limiter(j,i)=r_limiter(j,i)*r0x 
              z_limiter(j,i)=z_limiter(j,i)*r0x
           enddo
        enddo     
      endif

c------------------------------------------------------------------------
c$$$      if (i_genray_in_transformation.eq.1) then 
c$$$c--------transformation to MKSA system for genray.in file
c$$$         frqncy=frqncy*1.d+9   !from GHZ to HZ
c$$$      endif
 
      write(i_unit,wave,iostat=kode)
      call check_read(kode,'wave')
c---------------------------------------------------------
c     The data for  the n_perpendicular angle scattering
c--------------------------------------------------------
      write(i_unit,scatnper,iostat=kode)
      call check_read(kode,'scatnper')
c---------------------------------------------------------
      if (istart.eq.1) then
c---------------------------------------------------------
c        start point is outside the plasma, ECR case,
c        writing of the data for EC cone
c----------------------------------------------------------                
c$$$         if (i_genray_in_transformation.eq.1) then 
c$$$c----------transformation to MKSA system for genray.in file
c$$$            do i=1,nconea
c$$$              powtot(i)=powtot(i)*1.0d+6   !from MWatt  to Watt
c$$$            enddo
c$$$          endif

         write(i_unit,eccone,iostat=kode) 
         call check_read(kode,'eccone')
      else 
         !istart= 2 or 3     
c-------------------------------------------------
c        start point is inside the plasma, LH or FW case
c        the reading of the data for LH grill
c---------------------------------------  
c$$$         if (i_genray_in_transformation.eq.1) then 
c$$$c----------transformation to MKSA system for genray.in file
c$$$            do i=1,ngrilla
c$$$              powers(i)= powers(i)*1.0d+6   !from MWatt  to Watt
c$$$            enddo
c$$$          endif

         write(i_unit,grill,iostat=kode) 
         call check_read(kode,'grill')

         if (ilaunch.eq.2) then
c-------------------------------------------------------
c           write rz_launch_grill_tab
c------------------------------------------------
            call write_rz_launch_grill_nml(i_unit,kode)  
            call check_read(kode,'rz_launch_grill_tab')
         endif
   
      endif !istart
c--------------------------------------------------------------------
      write(i_unit,dispers,iostat=kode)
      call check_read(kode,'dispers')     
c----------------------------------------------------------------------
      prmt1=prmt(1)
      prmt2=prmt(2)
      prmt3=prmt(3)
      prmt4=prmt(4)
      prmt6=prmt(6)
      prmt9=prmt(9)
      ndim1=ndim    !put argument ndim to namelist paramer ndim1

      write(i_unit,numercl,iostat=kode)
      call check_read(kode,'numercl')
c------------------------------------------------------------------
c$$$      if (i_genray_in_transformation.eq.1) then 
c$$$c--------transformation to MKSA system for genray.in file
c$$$         max_plot_freq=max_plot_freq*1.d9 !from GHZ to HZ
c$$$      endif

      write(i_unit,output,iostat=kode)
      call check_read(kode,'output')
c--------------------------------------------------------
      write(i_unit,plasma,iostat=kode)
      call check_read(kode,'plasma') 
c-----------------------------------------------------
      write(i_unit,species,iostat=kode)
      call check_read(kode,'species')
c---------------------------------------------------------------------
c     plasma component charges charge(i)=mod(charge(i)/charge_electron)
c-----------------------------------------------------
c     plasma components mass dmas(i)=Mass(i)/Mass_electron
c----------------------------------------------------- 
c     calculation of nbulk1
      if(((izeff.eq.0).or.(izeff.eq.2)).or.(izeff.eq.3)) then
c        izeff=0, zeff will be calculated using the given ions;
c                 electron density will be calculated using ion's densities;
c             =1  ion densities nbulk and nbulk-1 will be calculated  using
c                 Zeff, electon density and ion's densities(i), i=2,nbulk-1;
c        izeff=2, zeff will not coincide with the plasma components
c             =3  it uses eqdsk pres (pressure) and ions densities_i
c                 for i=2,... nbulk
c                 Let temperature T_E=T_i
c                 pres=dens1(k,1)*temp1(k,1)+
c                      Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
c                 In this case we will calculate Zeff(rho),
c                 dens_electron(rho) and T_e(rho)=T_i(rho)
c             =4  it uses eqdsk pres (pressure), zeff,ions densities
c                 for i=2,... nbulk-2 (nbulk>2) and T_e(rho),T_i(rho)
c                 pres=dens1(k,1)*temp1(k,1)+
c                      Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
c                 In this case we will calculate dens_electron(rho) and
c                 ion densities for i=nbulk and i=nbulk-1)
         nbulk1=nbulk         
      else
c        izeff=1, zeff is given, the ions component will be calculated
         if (nbulk.eq.1) nbulk1=1
         if (nbulk.eq.2) then
	    nbulk1=2
	    write(*,*)'nbulk=2 Zeff must be equal charge(2) control it'
	    write(*,*)'use the option izeff=0'
	    stop
	 endif
         if (nbulk.gt.2) nbulk1=nbulk-2
      endif !izeff
c------------------------------------------------------------------
c     small radius uniform mesh for plasma profiles
c     It will be recalculated in dinit.f for common  six_no_nml.i
c-----------------------------------------------------
      h=1.d0/(ndens-1)
      do i=1,ndens
        rhom(i)=h*(i-1)
      enddo
c------------------------------------------------------------------
c     The parameters for the density fluctuations
      write(i_unit,varden,iostat=kode)
      call check_read(kode,'varden')
c------------------------------------------------------------------
c     The parameters for the density fluctuations
c------------------------------------------------------------------
      if(idens.eq.0) then
c------------------------------------------------------------------
c         Analytical radial profiles
c------------------------------------------------------------------
c$$$         if (i_genray_in_transformation.eq.1) then 
c$$$c-----------transformation to MKSA system for genray.in file
c$$$            do i=1,nbulk
c$$$              dense0(i)=dense0(i)*1.d+9   !from 10**19/m*3  to 1/m**3
c$$$              denseb(i)=dense0(i)*1.d+9   !from 10**19/m*3  to 1/m**3
c$$$            enddo
c$$$         endif

	 write(i_unit,denprof,iostat=kode)
         call check_read(kode,'denprof')

         write(i_unit,tpopprof,iostat=kode)
         call check_read(kode,'tpopprof')

         write(i_unit,vflprof,iostat=kode)
         call check_read(kode,'vfloprof')

	 write(i_unit,zprof,iostat=kode)
         call check_read(kode,'zprof')

	 write(i_unit,tprof,iostat=kode)
         call check_read(kode,'tprof')

      endif ! idens analytical


      if (idens.eq.1) then
c-----------------------------------------------------------------------
c        spline approximation of the density, temperature, zeff,
c        tpop=T_perp/T_parallel, vflow
c        radial profiles
c        input of the arrays on the radial mesh	from dtzprof.dat
c---------------------------------------------------------------------
c        set partner
c---------------------------------------------------------------------
         partner='disabled'
c---------------------------------------------------------------------
c        set den_scale and temp_scale =1
c----------------------------------------------------------------------
         do i=1,nbulk
            den_scale(i)=1.d0
            temp_scale(i)=1.d0 
         enddo
c--------------------------------------------------------------------         
         if(nonuniform_profile_mesh.eq.'enabled')then
c---------------------------------------------------------------------
c          nonuniform mesh profiles
c--------------------------------------------------------------------         
c$$$           if (i_genray_in_transformation.eq.1) then
c$$$c------------------------------------------------------
c$$$c             transformation of the density at nonuniform_profile_mesh 
c$$$c             from genray.in to genray.dat form 
c$$$c-------------------------------------------------------------------
c$$$              do k=1,ndensa
c$$$                do j=1,nbulka  
c$$$                   dens1_nonuniform(k,j)=dens1_nonuniform(k,j)*1.d19  !from 10**19/m**3 to 1/m**3
c$$$                enddo
c$$$              enddo
c$$$           endif

           call write_nonuniform_line_profile(i_unit,
     &     'dentab_nonuniform_line',nbulk,ndens,
     &     dens1_nonuniform,radii_nonuniform_dens1,nj_tab_dens1)
           
           call write_nonuniform_line_profile(i_unit,
     &     'temtab_nonuniform_line',nbulk,ndens,
     &     temp1_nonuniform,radii_nonuniform_temp1,nj_tab_temp1)
          
            call write_nonuniform_line_profile(i_unit,
     &     'tpoptab_nonuniform_line',nbulk,ndens,
     &     tpop1_nonuniform,radii_nonuniform_tpop1,nj_tab_tpop1)

           call write_nonuniform_line_profile(i_unit,
     &     'vflowtab_nonuniform_line',nbulk,ndens,
     &     vflow1_nonuniform,radii_nonuniform_vflow1,nj_tab_vflow1)

           call write_nonuniform_line_profile(i_unit,
     &     'zeftab_nonuniform_line',nbulk,ndens,
     &     zeff1_nonuniform,radii_nonuniform_zeff1,nj_tab_zeff1)
   
         else
c---------------------------------------------------------------------
c          uniform mesh profiles
c---------------------------------------------------------------------
c$$$           if (i_genray_in_transformation.eq.1) then
c$$$c------------------------------------------------------
c$$$c             transformation of the density at nonuniform_profile_mesh 
c$$$c             from genray.in to genray.dat form 
c$$$c-------------------------------------------------------------------
c$$$              do k=1,ndensa
c$$$                do j=1,nbulka  
c$$$                   dens1(k,j)=dens1(k,j)*1.d19  !from 10**19/m**3 to 1/m**3
c$$$                enddo
c$$$              enddo
c$$$           endif

           call write_uniform_column_profile(i_unit,
     &     'dentab',nbulk,
     &     ndens,dens1,prof1_uniform)

           call write_uniform_column_profile(i_unit,
     &     'temtab',nbulk,
     &     ndens,temp1,prof1_uniform)

           call write_uniform_column_profile(i_unit,
     &     'tpoptab',nbulk,
     &     ndens,tpop1,prof1_uniform)

           call write_uniform_column_profile(i_unit,
     &     'vflowtab',nbulk,
     &     ndens,vflow1,prof1_uniform)
          
           call write_uniform_column_profile(i_unit,
     &     'zeftab',nbulk,
     &     ndens,prof2,zeff1)
          
          if( i_lsc_approach.eq.1) then

             call write_uniform_column_profile(i_unit,
     &      'EdcTSCtab',1,
     &      n_psi_TSC,prof2,EdcTSC_1D)

             call write_uniform_column_profile(i_unit,
     &      'JparTSCtab',1,
     &      n_psi_TSC,prof2,JparTSC_1D)

          endif
        
         endif
      
      endif! idens=1
c-----------------------------------------------------------
c     data for reading 3D distribution function from diskf      
c     written by CQL3D 
c---------------------------------------------------------------
      write(i_unit,read_diskf,iostat=kode)
      call check_read(kode,'read_diskf')
c------------------------------------------------------------
c     write data for emission calculations
c---------------------------------------------------------
      write(i_unit,emission,iostat=kode)
      call check_read(kode,'emission')     
c---------------------------------------------------------
c     data for for EC cone vertex coordinates calculation.
c     This case is used for the optimal OX mode conversion.
c-----------------------------------------------------------
      write(i_unit,ox,iostat=kode)
      call check_read(kode,'ox')
c-----------------------------------------------------------------
c     write data for for adj calculations
c-----------------------------------------------------------------      
      write(i_unit,adj_nml,iostat=kode)
      call check_read(kode,'adj_nml')

cSAP090203
c-------------------------------------------------------------------     
c     write data for density profile outside LCFS
c------------------------------------------------------------------
      write(i_unit,edge_prof_nml,iostat=kode)
      call check_read(kode,'edge_prof_nml')

c-----------------------------------------------------------------------
c     write data for LSC approach calculations
c--------------------------------------------------------------------- 
      write(i_unit,edge_prof_nml,iostat=kode)
      call check_read(kode,'edge_prof_nml')
c---------------------------------------------------------------------
      close(i_unit)      
      write(*,*)'in prepare_genray_input.f end reading genray.in file'
c     end of writing genray.in or genray.dat file

      return

      end !!!!!!!!!!!!! write_all_namelists (starts: l 1403) !!!!!!!!!



c     ****************CHECK_READ********************************
      subroutine check_read(iostat,name_of_namelist)
c     check the sign of iostat from operator read
c     iostat < 0   the end of input file was detected
c     iostat=0     reading has complited succefully
c     iostat>0     an error has occurred
c
c     input
      integer iostat
      character*(*) name_of_namelist

c      write(*,*)'in equilib.f in check_read name_of_namelist = '
c     .,name_of_namelist,'  iostat = ', iostat      

      if (iostat.gt.0) then
         write(*,1)name_of_namelist,name_of_namelist
 1       format('check_read has found the positive value of iostat',/,
     .   'The error has occurred in reading namelist = ',A,/,
     .   'Check input data in ',A,' in the file genray.in')
         stop
      endif

      if (iostat.lt.0) then
         write(*,*)
         write(*,*)'**************************************************'
         write(*,2) name_of_namelist,name_of_namelist,name_of_namelist
 2       format('check_read has found the negative value of iostat',/,
     .   'The end of input in the given namelist ',A,' was detected.',/,
     .   'It can be that namelist ',A,' is absent in genray.in file',/,
     .   'Check data in ',A,' in the file genray.in')
         write(*,*)'**************************************************'
cBH050322         pause     !Not necessary in some cases, since
                            !defaults may be OK.
      endif

      return
      end        !!!!!!!check_read!!!!!!!!!!!!!!!!!!!!


c     
c These functions are modified from onetwo for spline density profiles
c End of package at line ~2352.


      subroutine icsicu1 (x, y, nx, bpar, ct, ic, ier)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- modified version of IMSL subroutine named ICSICU
c
c ----------------------------------------------------------------------
c
c   computer            - dec10/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - interpolatory approximation by cubic splines
c                           with arbitrary second derivative end
c                           conditions.
c
c   usage               - call icsicu1 (x, y, nx, bpar, c, ic, ier)
c
c   arguments    x      - vector of length nx containing the abscissae
c                           of the nx data points (x(i),y(i)) i = 1,...,
c                           nx. (input) x must be ordered so that
c                           x(i) .lt. x(i+1).
c                y      - vector of length nx containing the ordinates
c                           (or function values) of the nx data points.
c                           (input)
c                nx     - number of elements in x and y. (input) nx
c                           must be .ge. 2.
c                bpar   - vector of length 4 containing the end
c                           condition parameters. (input)
c                           2.0*spp(1)+bpar(1)*spp(2) = bpar(2),
c                           bpar(3)*spp(nx-1)+2.0*spp(nx) = bpar(4),
c                           where spp(i) = second derivative of the
c                           cubic spline function s evaluated at x(i).
c                c      - spline coefficients. (output) c is an nx-1 by
c                           3 matrix. the value of the spline
c                           approximation at t is
c                           s(t) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
c                           where x(i) .le. t .lt. x(i+1) and
c                           d = t-x(i).
c                ic     - row dimension of matrix c exactly as
c                           specified in the dimension statement in
c                           the calling program. (input)
c                ier    - error parameter. (output)
c                         terminal error
c                           ier = 129, ic is less than nx-1
c                           ier = 130, nx is less than 2.
c                           ier = 131, input abscissa are not ordered
c                             so that x(1) .lt. x(2) ... .lt. x(nx).
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. IMSL routines - uertst1,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through IMSL routine uhelp
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - IMSL warrants only that IMSL testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c ----------------------------------------------------------------------
c
c     specifications for arguments
c
      integer            nx,ic,ier
c      real*8             x(nx),y(nx),bpar(4),c(ic,3)
      real*8             x(*),y(*),bpar(4),ct(3,*)
c
c     specifications for local variables
c
      integer            i,j,nxm1
      real*8             dx,dxj,dxjp1,dxp,dyj,dyjp1,half,one,pj,
     .                   six,sixi,two,yppa,yppb,zero
      equivalence        (dxj,yppb),(pj,sixi),(dxjp1,yppa)
      data               zero/0.0/,half/0.5/,one/1.0/,
     .                   two/2.0/,six/6.0/
c
      ier = 0
c
c     check error conditions
c
      nxm1 = nx-1
      if (ic .lt. nxm1)  go to 30
      if (nx .lt. 2   )  go to 35
      if (nx .eq. 2   )  go to 10
c
c     compute coefficients and right hand side of the tridiagonal
c     system defining the second derivatives of the spline interpolant for (x,y)
c
c     c(j,1) = lambda(j)
c     c(j,2) = mu(j)
c     c(j,3) = d(j)
c
      dxj = x(2)-x(1)
      if (dxj .le. zero)  go to 40
      dyj = y(2)-y(1)
      do 5 j=2,nxm1
         dxjp1 = x(j+1)-x(j)
         if (dxjp1 .le. zero)  go to 40
         dyjp1 = y(j+1)-y(j)
         dxp = dxj+dxjp1
         ct(1,j) = dxjp1/dxp
         ct(2,j) = one-ct(1,j)
         ct(3,j) = six*(dyjp1/dxjp1-dyj/dxj)/dxp
         dxj = dxjp1
         dyj = dyjp1
    5 continue
c
c     factor the tridiagonal matrix and solve for u
c
c     ct(2,j)  = u(j)
c     ct(1,j)  = q(j)
c     bpar(1) = lambda(1)
c     bpar(2) = d(1)
c     bpar(3) = mu(nx)
c     bpar(4) = d(nx)
c
   10 ct(1,1) = -bpar(1)*half
      ct(2,1) = bpar(2)*half
      if (nx .eq. 2)  go to 20
      do 15 j=2,nxm1
         pj = ct(2,j)*ct(1,j-1)+two
         ct(1,j) = -ct(1,j)/pj
         ct(2,j) = (ct(3,j)-ct(2,j)*ct(2,j-1))/pj
   15 continue
c
c     solve for cubic coefficients of spline interpolant
c     c(j,1), c(j,2), and c(j,3)
c
   20 yppb = (bpar(4)-bpar(3)*ct(2,nxm1))/(bpar(3)*ct(1,nxm1)+two)
      sixi = one/six
      do 25 i=1,nxm1
         j = nx-i
         yppa = ct(1,j)*yppb+ct(2,j)
         dx = x(j+1)-x(j)
         ct(3,j) = sixi*(yppb-yppa)/dx
         ct(2,j) = half*yppa
         ct(1,j) = (y(j+1)-y(j))/dx-(ct(2,j)+ct(3,j)*dx)*dx
         yppb = yppa
   25 continue
      go to 9005
   30 ier = 129
      go to 9000
   35 ier = 130
      go to 9000
   40 ier = 131
c
c 9000 call uertst1 (ier, 'icsicu1')
c 9000 write(*,*)'icsicu1 ier=',ier
 9000 continue  !Need to make this change to accomodate nowrite
                !option, which assumes write(*,*) without a label.
      write(*,*)'icsicu1 ier=',ier
 9005 return
c
      end

cSm060825
c      subroutine intrp_adp (r_in, f_in,n_in, r_out, f_out, n_out)
      subroutine intrp_adp (n_in,r_in,f_in,ct,n_out,r_out,f_out)

c
c -------------------------------------------------------------------
c     interpolate f_in onto f_out using splines or ?
c -------------------------------------------------------------------
c
cSm060825
c      USE param
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'   ! kj
c
c      dimension  r_in(*), f_in(*), r_out(*), f_out(*)
c      dimension  bpar(4), c(kj,3)
cSm060825
      integer n_in,n_out
c      dimension  r_in(n_in), f_in(n_in), r_out(n_out), f_out(n_out)
c      dimension  bpar(4), c(n_in-1,3)
      dimension  r_in(*), f_in(*), r_out(*), f_out(*)
      dimension  bpar(4), ct(3,*)

      kj=n_in-1  
  
c
      bpar(1) = 1.0       ! zero gradient at rho = 0
      bpar(2) = 6.0 *(f_in(2)-f_in(1))/((r_in(2)-r_in(1))**2)
      bpar(3) = 0.0       ! natural at rho = 1
      bpar(4) = 0.0
      ier     = 0
      call icsicu1 (r_in, f_in, n_in, bpar, ct, kj, ier)

      write(*,*)'----in subroutine intrp_adp-------------'
c      write(*,*)'n_out',n_out  
c      write(*,*)'n_in',n_in 
c      write(*,*)'kj',kj
c      write(*,*)'r_in',r_in
c      write(*,*)'f_in',f_in
c      write(*,*)'ct',ct
      write(*,*)'ier=',ier

      if (ier .ne. 0)  go to 10
cSm060825
     
c      do n=1,n_in 
c       write(*,*)'r_in(n),f_in(n)',r_in(n),f_in(n)
c      enddo

c      do n=1,kj 
c       do j=1,3
c         write(*,*)'n,j,ct(j,n)',n,j,ct(j,n)
c       enddo
c      enddo
c      write(*,*)'r_out',r_out
      call icsevu1 (r_in, f_in, n_in, ct, kj, r_out, f_out, n_out, ier)
c      write(*,*)'f_out',f_out
cSm060825
c      write(*,*)'n_out',n_out
c      do n=1,n_out           
c         write(*,*)'n,r_out(n),f_out(n)',n,r_out(n),f_out(n)
c      enddo

   10 if (ier .eq. 0)  return
cSm080825
c      call STOP ('subroutine INTRP_ADP: non-zero IER from IMSL', 247)
       write(*,*)'subroutine INTRP_ADP: non-zero IER from IMSL'
c
      end
      subroutine icsevu1 (x, y, nx, ct, ic, u, s, m, ier)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- modified version of IMSL subroutine named ICSEVU
c
c ----------------------------------------------------------------------
c
c   computer            - dec10/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - evaluation of a cubic spline
c
c   usage               - call icsevu1 (x, y, nx, c, ic, u, s, m, ier)
c
c   arguments    x      - vector of length nx containing the abscissae
c                           of the nx data points (x(i),y(i)) i = 1,...,
c                           nx (input). x must be ordered so that
c                           x(i) .lt. x(i+1).
c                y      - vector of length nx containing the ordinates
c                           (or function values) of the nx data points
c                           (input).
c                nx     - number of elements in x and y (input).
c                           nx must be .ge. 2.
c                c      - spline coefficients (input). c is an nx-1 by
c                           3 matrix.
c                ic     - row dimension of matrix c exactly as
c                           specified in the dimension statement
c                           in the calling program (input).
c                           ic must be .ge. nx-1
c                u      - vector of length m containing the abscissae
c                           of the m points at which the cubic spline
c                           is to be evaluated (input).
c                s      - vector of length m (output).
c                           the value of the spline approximation at
c                           u(i) is
c                           s(i) = ((c(j,3)*d+c(j,2))*d+c(j,1))*d+y(j)
c                           where x(j) .le. u(i) .lt. x(j+1) and
c                           d = u(i)-x(j).
c                m      - number of elements in u and s (input).
c                ier    - error parameter (output).
c                         warning error
c                           ier = 33, u(i) is less than x(1).
c                           ier = 34, u(i) is greater than x(nx).
c
c                           ********************************************
c                           output of warning errors has been suppressed
c                           ********************************************
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. IMSL routines - uertst1,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through IMSL routine uhelp
c
c   remarks  1.  the routine assumes that the abscissae of the nx
c                data points are ordered such that x(i) is less than
c                x(i+1) for i = 1,...,nx-1. no check of this condition
c                is made in the routine. unordered abscissae will cause
c                the algorithm to produce incorrect results.
c            2.  the routine generates two warning errors. one error
c                occurs if u(i) is less than x(1), for some i in the
c                the interval (1,m) inclusively. the other error occurs
c                if u(i) is greater than x(nx), for some i in the
c                interval (1,m) inclusively.
c            3.  the ordinate y(nx) is not used by the routine. for
c                u(k) .gt. x(nx-1), the value of the spline, s(k), is
c                given by
c                 s(k) = ((c(nx-1,3)*d+c(nx-1,2))*d+c(nx-1,1))*d+y(nx-1)
c                where d = u(k)-x(nx-1).
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - IMSL warrants only that IMSL testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c ----------------------------------------------------------------------
c
c     specifications for arguments
c
      integer            nx,ic,m,ier
c      real*8             x(nx),y(nx),c(ic,3),u(m),s(m)
      
      real*8             x(*),y(*),u(*),s(*)
      real*8             ct(3,*)
c
c     specifications for local variables
c
      integer            i,jer,ker,nxm1,k
      real*8             d,dd,zero
      data               i/1/, zero/0.0/

      write(*,*)'in subroutine icsevu1 ic=',ic
c      write(*,*)'nx',nx
c      write(*,*)'x',x
c      write(*,*)'y',y
c      write(*,*)'c',c
c      write(*,*)'u',u
c      write(*,*)'s',s
c      do i=1,ic
c        do j=1,3
c          write(*,*)'i,j,c(i,j)',i,j,c(i,j)
c        enddo
c      enddo
 
c
c     first executable statement
c
      jer = 0
      ker = 0
      if (m .le. 0)  go to 9005
      nxm1 = nx-1
      if (i .gt. nxm1)  i = 1
c
c     evaluate spline at m points
c
      do 40 k=1,m
c
c        find the proper interval
c
         d = u(k)-x(i)
         if (d) 5, 25, 15
    5    if (i .eq. 1)  go to 30
         i = i-1
         d = u(k)-x(i)
         if (d) 5, 25, 20
   10    i = i+1
         d = dd
   15    if (i .ge. nx)  go to 35
         dd = u(k)-x(i+1)
         if (dd .ge. zero)  go to 10
         if ( d .eq. zero)  go to 25
c
c        perform evaluation
c
c   20    s(k) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
   20    s(k) = ((ct(3,i)*d+ct(2,i))*d+ct(1,i))*d+y(i)
         go to 40
   25    s(k) = y(i)
         go to 40
c
c        u(k) < x(1)
c
   30    jer = 33
         go to 20
c
c        u(k) > x(nx)
c
   35    if (dd .gt. zero)  ker = 34
         d = u(k) - x(nxm1)
         i = nxm1
         go to 20
c
   40 continue
c
      ier = MAX0 (jer, ker)
c
cSm060825
c****  if (jer .gt. 0)  call uertst1 (jer, 'icsevu1')
c****  if (ker .gt. 0)  call uertst1 (ker, 'icsevu1')
****  if (jer .gt. 0)  write(*,*)'icsevu1 jer=',jer
****  if (ker .gt. 0)  write(*,*)'icsevu1 ker=',ker
c      write(*,*)'s',s
c
 9005 return
c
      end         !!!!!!!!!End modified spline package!!!!!!

  

      subroutine ddcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      real*8 dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end


     
      subroutine pack21(a,ibot,itop,jbot,jtop,b,iy,jx)
      implicit integer (i-n), real*8 (a-h,o-z)
c.......................................................................
c     It sometimes becomes necessary to take a
c     2-D array dimensioned ibot:itop by jbot:jtop
c     and repack it as though it were
c     dimensioned 1:iy by 1:jx, starting at a(1,1).
c     This routine does this, transfering relevant data
c     from array a to b.
c.......................................................................

      save
      dimension a(ibot:itop,jbot:jtop)
      dimension b(iy*jx)
c      write(*,*)'pack21 ibot,itop,jbot,jtop',ibot,itop,jbot,jtop
c      write(*,*)'pack21 iy,jx',iy,jx

      do 1 j=1,jx
        i1=(j-1)*iy+1
c        call scopy(iy,a(1,j),1,b(i1),1)

c         write(*,*)'pack21 a(1,j)'
c         do i=1,iy
c           write(*,*)'i,a(i,j)',i,a(i,j)
c         enddo
 
        call ddcopy(iy,a(1,j),1,b(i1) ,1)

c        write(*,*)'pack21 b(i1,j)'
c        do k=i1,i1+jx-1
c           write(*,*)'k,b(k)',k,b(k)
c        enddo

 1    continue
      return
      end


      subroutine ddcopy_integer(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c     dx and dy are integer
c
      implicit none 
cBH111031      integer*8 dx(*),dy(*)
      integer dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end

      subroutine pack21_integer(a,ibot,itop,jbot,jtop,b,iy,jx)
      implicit none
c.......................................................................
c     It sometimes becomes necessary to take a
c     2-D array dimensioned ibot:itop by jbot:jtop
c     and repack it as though it were
c     dimensioned 1:iy by 1:jx, starting at a(1,1).
c     This routine does this, transfering relevant data
c     from array a to b.
c 
c     a and b are integer
c.......................................................................
c-----input
      integer a,ibot,itop,jbot,jtop,iy,jx
      dimension a(ibot:itop,jbot:jtop)
c-----output
      integer b
      dimension b(iy*jx)
c-----local
      integer j,i1


      save
      
c      write(*,*)'pack21_integer ibot,itop,jbot,jtop',ibot,itop,jbot,jtop
c      write(*,*)'pack21_integer iy,jx',iy,jx

      do 1 j=1,jx
        i1=(j-1)*iy+1
c        call scopy(iy,a(1,j),1,b(i1),1)

c         write(*,*)'pack211 a(1,j)'
c         do i=1,iy
c           write(*,*)'i,a(i,j)',i,a(i,j)
c         enddo
 
        call ddcopy_integer(iy,a(1,j),1,b(i1) ,1)

c        write(*,*)'pack21_integer b(i1,j)'
c        do k=i1,i1+jx-1
c           write(*,*)'k,b(k)',k,b(k)
c        enddo

 1    continue
      return
      end

c
c

      subroutine pack21_dentab(a,ndensa,nbulka,b,ndens,nbulk)
c-------------------------------------------------------------
c     pack 2d array a(ndensa,nbulka) to 1D array b(ndens*nbulk)
c     to write profiles plasma profiles 
c     (dentab,temptab,tpoptab,vlowtab) to genray.dat file
c-------------------------------------------------------------
      implicit none
c-----input
      integer ndensa,nbulka,ndens,nbulk
      real*8 a(ndensa,nbulka)

c-----output
      real*8 b(ndensa*nbulka)

c-----locals
      integer i,j,k         

      k=0
      do j=1,ndens
         do i=1,nbulk
           k=k+1
           b(k)=a(j,i)
         enddo
      enddo

      return
      end

 
      subroutine read_nonuniform_line_profile(i_unit,
     &nametab,nbulk,ndens,
     &prof2_uniform,
     &prof2_nonuniform,radii2_nonuniform,nj_tab_nonuniform,kode)
!     [end at l ~2839].
c------------------------------------------------------------------------
c     Reads nonuniform profiles in line form.
c     nametab can be = 'dentab','temptab','tpoptab','vflowtab','zeftab'.
c     The profile will be in output file prof2_uniform(nbulka,ndensa)
c     It can read profiles:
c     prof_2_uniform= dens1,temp1,tpop1,vflow1,zeff1
c
c     Put profiles to nonuniform grid arrays:
c     prof2_nonuniform (ndensa,nbulka) - profiles
c     radii2_nonuniform(ndensa,nbulka) - small radius
c     nj_tab_nonuniform(nbulka)        - number of radial points
c
c     It puts the profile at uniform grid
c     rho(j)=(j-1)/(ndens-1) i=1,...,ndens 
c     to prof2_uniform(ndensa,nbulka)
c
c     output: kode  < 0   the end of input file was detected
c                   = 0     reading has complited succefully
c                   > 0     an error has occurred
c------------------------------------------------------------------------   
      implicit none
      include 'param.i'
   
c-----input
      character(*) nametab   ! name of namelist
      integer
     &i_unit,                ! is the number of opened input file 
     &nbulk,                 ! nbulk>=1 is a number of plasma components
                             !          It should be nbulk.le.nbulka

     &ndens ! The number of points used for spline plasma profiles
            ! for uniform grid
            ! and the number of input points for the radial plasma
            ! profiles for nonuniform grid for 
            ! nonuniform_profile_mesh='disabled'
c-----output
      real*8
     &prof2_uniform(ndensa,nbulka),        ! profile at uniform grid
     &prof2_nonuniform(ndensa,nbulka),     ! plasma profile at nonuniform grid
     &radii2_nonuniform(ndensa,nbulka)     ! radii2_nonuniform

      integer nj_tab_nonuniform(nbulka)    ! number of mesh points
       
      integer kode          !  < 0   the end of input file was detected
                            !  =0     reading has complited succefully
                            !  >0     an error has occurred
c-----locals
c      real*8
c     &zeff1(ndensa)                ! profile at uniform grid for zeff

      real*8 radii_1_in(ndensa),prof_1_in(ndensa)
      real*8 radii_1_out(ndensa),prof_1_out(ndensa)
      integer i,j
      real*8 ct(3,ndensa-1) !work array for spline
  

c-----------------------------------------------------------------
c      namelists for all table plasma profiles at non uniform
c      radial mesh written by rows:
c-----------------------------------------------------------------
      real*8 prof_2d(ndensa,nbulka),radii_2d(ndensa,nbulka)
      integer nj_tab(nbulka)
c      include 'name_non_uniform_mesh_profiles_line.i'
      namelist /dentab_nonuniform_line/    nj_tab,prof_2d,radii_2d
      namelist /temtab_nonuniform_line/    nj_tab,prof_2d,radii_2d
      namelist /tpoptab_nonuniform_line/   nj_tab,prof_2d,radii_2d
      namelist /vflowtab_nonuniform_line/  nj_tab,prof_2d,radii_2d
      namelist /zeftab_nonuniform_line /   nj_tab,prof_2d,radii_2d

      write(*,*)'in subroutine  read_nonuniform_line_profile'
    
      write(*,*)'nametab=',nametab
      write(*,*)'nbulk,ndens',nbulk,ndens
cSm070720
      rewind(unit=i_unit) 
c--------------------------------------------------------------------
c     read namelist 'nametab'
c--------------------------------------------------------------------
      call ibcast(nj_tab_nonuniform,0,nbulka)
      call bcast(radii2_nonuniform(1,1),0.d0,nbulka*ndensa)
      call bcast(prof2_nonuniform(1,1),0.d0,nbulka*ndensa)
      call bcast(prof2_uniform(1,1),0.d0,nbulka*ndensa)

      write(*,*)'in read_nonuniform_line_profile nametab',nametab

c-----------------------------------------------------------
c     read profile tables at nonuniform radial mesh 
c     by rows
c-----------------------------------------------------------
      call ibcast(nj_tab,0,nbulka)
      call bcast(prof_2d(1,1),0.d0,ndensa*nbulka)
      call bcast(radii_2d(1,1),0.d0,ndensa*nbulka)

      if (nametab.eq.'dentab_nonuniform_line') then
         read(i_unit,dentab_nonuniform_line,iostat=kode)
         call check_read(kode,'dentab_nonuniform_line')
         write(*,dentab_nonuniform_line)
      endif

      if (nametab.eq.'temtab_nonuniform_line') then
         read(i_unit,temtab_nonuniform_line,iostat=kode)
         call check_read(kode,'temtab_nonuniform_lined')
         write(*,temtab_nonuniform_line)
cSAP090311
          write(*,*)'after temtab_nonuniform_line kode',kode

      endif

      if (nametab.eq.'tpoptab_nonuniform_line') then
cSAP090311
         write(*,*)'in read_nonuniform_line_profile nametab'
         write(*,*)'before tpoptab_nonuniform_line'
         call bcast(prof2_nonuniform(1,1),1.d0,nbulka*ndensa)
         call bcast(prof2_uniform(1,1),1.d0,nbulka*ndensa)
         read(i_unit,tpoptab_nonuniform_line,iostat=kode)
         call check_read(kode,'tpoptab_nonuniform_line')
         write(*,*)'after tpoptab_nonuniform_line kode',kode

         write(*,tpoptab_nonuniform_line)

      endif

      if (nametab.eq.'vflowtab_nonuniform_line') then
         call bcast(prof2_nonuniform(1,1),0.d0,nbulka*ndensa)
         call bcast(prof2_uniform(1,1),0.d0,nbulka*ndensa)
         read(i_unit,vflowtab_nonuniform_line,iostat=kode) 
         call check_read(kode,'vflowtab_nonuniform_line')
         write(*,vflowtab_nonuniform_line)
      endif

      if (nametab.eq.'zeftab_nonuniform_line') then
          read(i_unit,zeftab_nonuniform_line,iostat=kode)
          call check_read(kode,'zeftab_nonuniform_line')
          write(*,zeftab_nonuniform_line)       
      endif

c--------------------------------------------------------      
      write(*,*)'in  read_nonuniform_line_profile after
     & read(1,',nametab,')'
          
c-------------------------------------------------------------------
c     check the number of input radial mesh points nj_tab(i)
c     for each species i=1,...,nbulk
c--------------------------------------------------------------------

      do i=1,nbulk
         if(nj_tab(i).gt.ndensa) then
           write(*,*)'nj_tab(i).gt.ndensa'
           write(*,*)'i=',i,'nj_tab(i)=',nj_tab(i),'ndensa=',ndensa
           write(*,*)'it should be nj_tab(i).le.ndensa'
           write(*,*)'Please increase ndensa in param.i and recompile'
           stop 'in read_nonuniform_line_profile'
         endif
         if(nj_tab(i).lt.0) then
           write(*,*)'nj_tab(i).lt.ndensa'
           write(*,*)'i=',i,'nj_tab(i)=',nj_tab(i)
          write(*,*)'it should be nj_tab(i).ge.0'
        write(*,*)'Please change nj_tab(i) in genray.in or genray.dat'
          stop 'in read_nonuniform_line_profile'
         endif
      enddo

      write(*,*)'nj_tab',nj_tab

c-----------------------------------------------------------------------
c     check that radial knots are positive and montonic
c----------------------------------------------------------------------
      do i=1,nbulk
           
         if(radii_2d(1,i).lt.0.d0) then
            write(*,*)'in genray.dat or genray.in namelist',nametab
            write(*,*)'radii_2d(1,i).lt.0.d0'
            write(*,*)'it should be radii_2d(1,i).gt.0.d0'
            write(*,*)'Please correct input radii_2d'
            stop 'in read_nonuniform_line_profile'
         endif  

         do j=2,nj_tab(i)  
            if(radii_2d(j-1,i).gt.radii_2d(j,i)) then
               write(*,*)'*****************************************'
               write(*,*)'in genray.dat or genray.in namelist',nametab
               write(*,*)'has nonmontonic radii knots at i,j',i,j
               write(*,*)'radii_2d(j-1,i).gt.prof_radii_2d(j,i)'
               write(*,*)'Please correct input radii_2d'
               stop 'in read_nonuniform_line_profile'
            endif  
         enddo
      enddo
c----------------------------------------------------------------------
c     put namelist's arrays to output arguments
c---------------------------------------------------------------------
      do i=1,nbulk
         nj_tab_nonuniform(i)=nj_tab(i)
         do j=1,nj_tab(i)
c            radii2_nonuniform(i,j)=radii_2d(j,i)
c            prof2_nonuniform(i,j)=prof_2d(j,i)
            radii2_nonuniform(j,i)=radii_2d(j,i)
            prof2_nonuniform(j,i)=prof_2d(j,i)
         enddo
      enddo
c-----------------------------------------------------------------------
c     put input profiles data y_out from non-uniform grid x,y
c     to uniform grid x_out using spline
c---------------------------------------------------------------------
      do j=1,ndens
         radii_1_out(j)=1.d0*(j-1)/(ndens-1)
      enddo     

      if (nametab.eq.'zeftab_nonuniform_line') then
         do j=1,nj_tab(1)
            radii_1_in(j)=radii_2d(j,1)
            prof_1_in(j)=prof_2d(j,1)
         enddo

         call put_to_uniform_mesh(nj_tab(1),ndensa,radii_1_in,
     &   prof_1_in,ct,ndens,radii_1_out,prof_1_out)

         do j=1,ndens         
            prof2_uniform(j,1)=prof_1_out(j)
         enddo  
      else
         do i=1,nbulk
            do j=1,nj_tab(i)
               radii_1_in(j)=radii_2d(j,i)
               prof_1_in(j)=prof_2d(j,i)

            enddo
            
c            write(*,*)'before put_to_uniform_mesh'
c            write(*,*)'i,nj_tab(i)',i,nj_tab(i)
c            write(*,*)'radii_1_in(j)',(radii_1_in(j),j=1,nj_tab(i))
c            write(*,*)'prof_1_in(j)',(prof_1_in(j),j=1,nj_tab(i))

            call put_to_uniform_mesh(nj_tab(i),ndensa,radii_1_in,
     &      prof_1_in,ct,ndens,radii_1_out,prof_1_out)

c            write(*,*)'after put_to_uniform_mesh'
c            write(*,*)'radii_1_out(j)',(radii_1_out(j),j=1,ndens)
c            write(*,*)'prof_1_out(j)',(prof_1_out(j),j=1,ndens)

            do j=1,ndens         
               prof2_uniform(j,i)=prof_1_out(j)
            enddo  
         enddo 
      endif 

c      do i=1,nbulk
c         do j=1,ndens         
c            write(*,*)'i,j,prof2_uniform(j,i)',i,j,prof2_uniform(j,i)
c         enddo
c      enddo 

      return
      end         !!!!!!!!!read_nonuniform_line_profile!!!!!!!!


      subroutine write_uniform_column_profile(i_unit,nametab,nbulk,
     & ndens,prof2_uniform,prof1_uniform)
!     [end at line ~2971].
c-----------------------------------------------------------------------
c     Writes the radial profile 'nametab' given at uniform space grid
c     rho(j)=(j-1)/(ndens-1) i=1,...,ndens, with namelist format
c     into file specified by i_unit.
c    
c     nametab can be = 'dentab','temptab','tpoptab','vflowtab','zeftab'.
c
c     The input profile will be in output file prof2_uniform(nbulka,ndensa)
c     The following profiles can be written: dens1,temp1,tpop1,vflow1,
c     or in prof1_unifom: zeff1, EdcTSC_1D(1:n_psi_TSC), 
c     JparTSC_1D(1:n_psi_TSC)    
c------------------------------------------------------------------------
      implicit none
      include 'param.i'
   
c-----input      
      character(*) nametab  ! name of namelist
      integer 
     &i_unit,    ! is the number of opened output file 
     &nbulk,     ! nbulk>=1 is a number of plasma components
                 !          It should be nbulk.le.nbulka

     &ndens ! The number of points used for spline plasma profiles
            ! at uniform grid
            ! and the number of input points for the radial plasma
            ! profiles for nonuniform grid at 
            ! nonuniform_profile_mesh='disabled'
            ! ndens=n_psi_TSC for  EdcTSC_1D(1:n_psi_TSC)
            ! JparTSC_1D(1:n_psi_TSC)
      real*8 prof2_uniform(ndensa,nbulka),  ! plasma profile at uniform grid
                                            ! profiles of density,temperature,
                                            ! tpop,vflow
     &prof1_uniform(ndensa)                 ! zeff profile
                                            ! EdcTSC_1D(1:n_psi_TSC)
                                            ! JparTSC_1D(1:n_psi_TSC)
c-----locals
c     Following are working arrays for the namelist input of density 
c     and temperature,... at nonuniform radial grids
      real*8 prof(nbulka*ndensa),
     &zeff1(ndensa)                ! profile at uniform grid for zeff

      real*8 prof2(nbulka,ndensa)
      integer i,j,kode
     
  
      namelist /dentab/ prof 
      namelist /temtab/ prof
      namelist /tpoptab/ prof
      namelist /vflowtab/ prof
      namelist /zeftab/ zeff1
      namelist /EdcTSCtab/ prof
      namelist /JparTSCtab/prof

      write(*,*)'in subroutine  write_uniform_column_profile'
      write(*,*)'nametab=',nametab
      write(*,*)'nbulk,ndens',nbulk,ndens


c--------------------------------------------------------------------
c     write namelist 'nametab'
c--------------------------------------------------------------------
c     Set default, in case no input data:
      call bcast(prof(1),1.d0,nbulka*ndensa)
c      call bcast(prof2(1,1),1.d0,nbulka*ndensa)
 
      write(*,*)'in write_uniform_column_profile, nametab= ',nametab

      if(nametab.eq.'dentab') then
        call pack21_dentab(prof2_uniform,ndensa,nbulka,prof,
     &                     ndens,nbulk)
        write(i_unit,dentab,iostat=kode)
c        write(*,dentab)
      endif

      if(nametab.eq.'temtab') then
        call pack21_dentab(prof2_uniform,ndensa,nbulka,prof,
     &                     ndens,nbulk)
        write(i_unit,temtab,iostat=kode)
c        write(*,temtab)
      endif

      if(nametab.eq.'tpoptab') then
        call pack21_dentab(prof2_uniform,ndensa,nbulka,prof,
     &                     ndens,nbulk)
        write(i_unit,tpoptab,iostat=kode)
c        write(*,tpoptab)
      endif

      if(nametab.eq.'vflowtab') then
        call pack21_dentab(prof2_uniform,ndensa,nbulka,prof,
     &     ndens,nbulk)
        write(i_unit,vflowtab,iostat=kode)
c        write(*,vflowtab)
      endif

      if(nametab.eq.'zeftab') then
        do j=1,ndens
           zeff1(j)=prof1_uniform(j)
        enddo
c        write(*,*)'zeff1',zeff1
        write(i_unit,zeftab,iostat=kode)
c        write(*,zeftab)        
      endif

      if(nametab.eq.'EdcTSCtab') then
        do j=1,ndens
           prof(j)=prof1_uniform(j)
        enddo
c        write(*,*)'EdcTSCtab prof1_uniform',prof1_uniform
        write(i_unit,EdcTSCtab,iostat=kode)
c        write(*,EdcTSCtab)        
      endif

      if(nametab.eq.'JparTSCtab') then
        do j=1,ndens
           prof(j)=prof1_uniform(j)
        enddo
c        write(*,*)'JparTSCtab'prof1_uniform',prof1_uniform
        write(i_unit,JparTSCtab,iostat=kode)
c        write(*,zeftab)        
      endif

      write(*,*)'in write_uniform_column_profile'
      write(*,*)'after write(i_unit,..,nametab,)'


    
      return
      end     !!!!!!!!!!!!!!write_uniform_column_profile!!!!!!!!
      
      subroutine read_uniform_column_profile(i_unit,
     & nametab,nbulk,
     & ndens,prof2_uniform,prof1_uniform,kode)
!     [end at line ~].
c-----------------------------------------------------------------------
c     Reads the radial profile 'nametab' in column form 
c     given at uniform grid: rho(j)=(j-1)/(ndens-1) i=1,...,ndens :
c    
c     nametab can be = 'dentab','temptab','tpoptab','vflowtab','zeftab'.
c
c     The profile will be in output file prof2_uniform(nbulka,ndensa)
c     It can read profiles:
c     prof2_uniform= dens1,temp1,tpop1,vflow1,
c     zeff1 =
c
c     profiles with names = 'EdcTSCtab', 'JparTSCtab'
c     are given at the uniform radial mesh at the center bin points
c     rho(j)=(j-0.5)/ndens i=1,...,ndens :
c     These cases are used to read profiles from TSC for lsc approach.
c     These profiles have n_psi_TSC points.
c     So, for reading of TSC profiles agruments are:
c     ndens = n_psi_TSC
c     nbulk=1
c     rof1_uniform_TSC = EdcTSC_1D or =JparTSC_1D
c
c     output: kode  < 0   the end of input file was detected
c                   = 0     reading has complited succefully
c                   > 0     an error has occurred
c------------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'lsc_approach_nml.i'
   
c-----input      
      character(*) nametab  ! name of namelist
      integer 
     &i_unit,    ! is the number of opened input file 
     &nbulk,     ! nbulk>=1 is a number of plasma components
                 !          It should be nbulk.le.nbulka

     &ndens ! The number of points used for spline plasma profiles
            ! at uniform grid
            ! and the number of input points for the uniform radial plasma
            ! profiles for nonunifoem grid at 
            ! nonuniform_profile_mesh='disabled'
c-----output
      real*8 prof2_uniform(ndensa,nbulka), ! plasma profile at uniform grid
                                           ! profiles of density,temperature,
                                           ! tpop,vflow
     &prof1_uniform(ndensa)                ! zeff1 profiles at uinform grid
                                           ! 1D profiles at uinform grid for
                                           ! ,EdcTSC_1D,JparTSC_1D

      integer kode          !  < 0   the end of input file was detected
                            !  =0     reading has complited succefully
                            !  >0     an error has occurred
c-----locals
c     Following are working arrays for the namelist input of density 
c     and temperature,... at nonuniform radii grids
      real*8 prof(nbulka*ndensa),
     &      zeff1(ndensa)      
      real*8 prof2(nbulka,ndensa)
      integer i,j
     
  
      namelist /dentab/ prof 
      namelist /temtab/ prof
      namelist /tpoptab/ prof
      namelist /vflowtab/ prof
      namelist /zeftab/ zeff1
      namelist /EdcTSCtab/ prof
      namelist /JparTSCtab/prof


      write(*,*)'in subroutine read_uniform_column_profile'
      write(*,*)'nametab= ',nametab
      write(*,*)'nbulk,ndens ',nbulk,ndens


c--------------------------------------------------------------------
c     read namelis 'nametab'
c--------------------------------------------------------------------
c     Set default, in case no input data:
      call bcast(prof(1),1.d0,nbulka*ndensa)
      call bcast(prof2(1,1),1.d0,nbulka*ndensa)
      call bcast(prof1_uniform(1),0.d0,ndensa)
      call bcast(zeff1(1),1.d0,ndensa)
 
cSAP080202
      call bcast(prof2_uniform(1,1),1.d0,ndensa*nbulka)
      write(*,*)'after bcast(prof2_uniform(1,1),1.d0,ndensa*nbulka)'
c      write(*,*)'prof2_uniform ',prof2_uniform

      rewind(unit=i_unit) 

      write(*,*)'in read_uniform_column_profile  1 nametab= ',nametab

      if(nametab.eq.'dentab') then
        read(i_unit,dentab,iostat=kode)  
        call check_read(kode,'dentab')
c        write(*,dentab)
      endif

      if(nametab.eq.'temtab') then
        read(i_unit,temtab,iostat=kode)
        call check_read(kode,'temtab')
c        write(*,temtab)
      endif

      if(nametab.eq.'tpoptab') then
        call bcast(prof(1),1.d0,nbulka*ndensa)
        call bcast(prof2(1,1),1.d0,nbulka*ndensa)
        call bcast(prof2_uniform(1,1),1.d0,ndensa*nbulka)
        read(i_unit,tpoptab,iostat=kode)
        call check_read(kode,'tpoptab')
c        write(*,tpoptab)
      endif

      if(nametab.eq.'vflowtab') then
        call bcast(prof(1),0.d0,nbulka*ndensa)
        call bcast(prof2(1,1),0.d0,nbulka*ndensa)
        call bcast(prof2_uniform(1,1),0.d0,ndensa*nbulka)
        read(i_unit,vflowtab,iostat=kode)
        call check_read(kode,'vflowtab')
c        write(*,*)' read_uniform_column_profile vlowtab'
c        write(*,vflowtab)
      endif

      if(nametab.eq.'zeftab') then
        call bcast(prof(1),1.d0,nbulka*ndensa) 
        call bcast(prof1_uniform(1),0.d0,ndensa)
        read(i_unit,zeftab,iostat=kode)
        call check_read(kode,'zeftab')
c        write(*,*)'read_uniform_column_profile zeftab'
c        write(*,zeftab)        
      endif

      if(nametab.eq.'EdcTSCtab') then
        call bcast(prof(1), EdcTSC_min,nbulka*ndensa) 
        call bcast(prof1_uniform(1),0.d0,ndensa)
        read(i_unit,EdcTSCtab,iostat=kode)
        call check_read(kode,'EdcTSCtab')
c        write(*,*)'read_uniform_column_profile EdcTSCtab'
c        write(*,EdcTSCtab)        
      endif

      if(nametab.eq.'JparTSCtab') then
        call bcast(prof(1),0.d0,nbulka*ndensa) 
        call bcast(prof1_uniform(1),0.d0,ndensa)
        read(i_unit,JparTSCtab,iostat=kode)
        call check_read(kode,'JparTSCtab')
c        write(*,*)'read_uniform_column_profile JparTSCtab'
c        write(*,JparTSCtab)        
      endif

      write(*,*)'in read_profile_tab_uniform_column
     & after read(i_unit,..,nametab,)'
c--------------------------------------------------------
c     Input radial profile is at uniform grid.
c     Put the input plasma profile from prof to prof2
c-----------------------------------------------------------------------      
      write(*,*)'nametab= ',nametab

      if (nametab.eq.'zeftab') then
         do j=1,ndensa
            prof1_uniform(j)=zeff1(j)
         enddo
      elseif 
     &  ((nametab.eq.'EdcTSCtab').or.(nametab.eq.'JparTSCtab')) then
          do j=1,ndens
             prof1_uniform(j)=prof(j)
          enddo
      else
         call pack12S(prof,ndensa*nbulka,
     &   nbulk,ndens,prof2,nbulka,ndensa)
         do i=1,nbulk
           do j=1,ndens
             prof2_uniform(j,i)=prof2(i,j)
           enddo
         enddo
      endif
     
c      write(*,*)'nametab=',nametab,' prof2',prof2
         

      do i=1,nbulk
         do j=1,ndens
            prof2_uniform(j,i)=prof2(i,j)
         enddo
      enddo

c      write(*,*)' nametab= ',nametab,' prof2_uniform= ',prof2_uniform

      write(*,*)'end of read_uniform_column_profile, nametab=',nametab

c      do i=1,nbulka
c         do j=1,ndensa
c            write(*,*)'j,i,prof2_uniform(j,i) ',j,i,prof2_uniform(j,i)
c         enddo
c      enddo

      return
      end        !!!!!!!!!read_uniform_column_profile!!!!!!!!!



      subroutine write_nonuniform_line_profile(i_unit,
     &nametab,nbulk,ndens,
     &prof2_nonuniform,radii2_nonuniform,nj_tab_nonuniform)
!     [end at line ~3322].
c------------------------------------------------------------------------
c     Writes nonuniform profiles in row form.
c     nametab can be = 'dentab','temptab','tpoptab','vflowtab','zeftab'.
c     The profile will be in output file prof2_uniform(nbulka,ndensa)
c     It can read profiles:
c     prof_2_uniform= dens1,temp1,tpop1,vflow1,zeff1
c
c     Put profiles to nonuniform grid arrays:
c     prof2_nonuniform (ndensa,nbulka) - profiles
c     radii2_nonuniform(ndensa,nbulka) - small radius
c     nj_tab_nonuniform(nbulka)        - number of radial points
c
c------------------------------------------------------------------------   
      implicit none
      include 'param.i'
   
c-----input
      character(*) nametab                 ! name of namelist
      integer
     &i_unit,    ! is the number of opened output file 
     & nbulk,    ! nbulk>=1 is a number of plasma components
                 !          It should be nbulk.le.nbulka

     &ndens ! The number of points used for spline plasma profiles
            ! at uniform grid
            ! and the number of input points for the radial plasma
            ! profiles for nonuniform grid at 
            ! nonuniform_profile_mesh='disabled'
      
      real*8
     &prof2_nonuniform(ndensa,nbulka),     ! plasma profile at nonuniform grid
     &radii2_nonuniform(ndensa,nbulka)     ! radii2_nonuniform

      integer nj_tab_nonuniform(nbulka)    ! number of mesh points
       
c-----locals
      integer i,j,kode
     
c-----------------------------------------------------------------
c      namelists for all table plasma profiles at nonuniform
c      radial mesh written by rows:
c-----------------------------------------------------------------
      real*8 prof_2d(ndensa,nbulka),radii_2d(ndensa,nbulka)
      integer nj_tab(nbulka)
c      include 'name_non_uniform_mesh_profiles_line.i'
      namelist /dentab_nonuniform_line/    nj_tab,prof_2d,radii_2d
      namelist /temtab_nonuniform_line/    nj_tab,prof_2d,radii_2d
      namelist /tpoptab_nonuniform_line/   nj_tab,prof_2d,radii_2d
      namelist /vflowtab_nonuniform_line/  nj_tab,prof_2d,radii_2d
      namelist /zeftab_nonuniform_line /   nj_tab,prof_2d,radii_2d

      write(*,*)'in subroutine write_nonuniform_line_profile'
      write(*,*)'input arguments:'
      write(*,*)'nametab=',nametab
      write(*,*)'nbulk,ndens',nbulk,ndens
      do i=1,nbulk
       write(*,*)'i=',i,'nj_tab_nonuniform(i)',nj_tab_nonuniform(i)
       do j=1,nj_tab_nonuniform(i)
         write(*,*)'j,i,radii2_nonuniform(j,i),prof2_nonuniform(j,i)',
     &              j,i,radii2_nonuniform(j,i),prof2_nonuniform(j,i)
       enddo
      enddo  
c-----------------------------------------------------
c     put input arguments to namelist's arrays
c---------------------------------------------------------------------
      do i=1,nbulk
         nj_tab(i)=nj_tab_nonuniform(i)
         do j=1,nj_tab(i)
           radii_2d(j,i)=radii2_nonuniform(j,i)
           prof_2d(j,i) =prof2_nonuniform(j,i)
         enddo
      enddo

c--------------------------------------------------------------------
c     write namelist 'nametab'
c--------------------------------------------------------------------
      write(*,*)'in write_nonuniform_line_profile'

c-----------------------------------------------------------
c     write profile tables at nonuniform radial mesh 
c     by rows
c-----------------------------------------------------------

      if (nametab.eq.'dentab_nonuniform_line') then
         write(i_unit,dentab_nonuniform_line,iostat=kode)
         write(*,dentab_nonuniform_line)
      endif

      if (nametab.eq.'temtab_nonuniform_line') then
         write(i_unit,temtab_nonuniform_line,iostat=kode)
         write(*,temtab_nonuniform_line)
      endif

      if (nametab.eq.'tpoptab_nonuniform_line') then
         write(i_unit,tpoptab_nonuniform_line,iostat=kode)
         write(*,tpoptab_nonuniform_line)
      endif

      if (nametab.eq.'vflowtab_nonuniform_line') then
         write(i_unit,vflowtab_nonuniform_line,iostat=kode)
         write(*,vflowtab_nonuniform_line)
      endif

      if (nametab.eq.'zeftab_nonuniform_line') then
          write(i_unit,zeftab_nonuniform_line,iostat=kode)
          write(*,zeftab_nonuniform_line)       
      endif

c--------------------------------------------------------      
      write(*,*)'in write_nonuniform_line_profile after
     &write(1,',nametab,')'
          
c-------------------------------------------------------------------
c     check the number of input radial mesh points nj_tab(i)
c     for each species i=1,...,nbulk
c--------------------------------------------------------------------

      do i=1,nbulk
         if(nj_tab(i).gt.ndensa) then
           write(*,*)'nj_tab(i).gt.ndensa'
           write(*,*)'i=',i,'nj_tab(i)=',nj_tab(i),'ndensa=',ndensa
           write(*,*)'it should be nj_tab(i).le.ndensa'
           write(*,*)'Please increase ndensa in param.i and recompile'
           stop 'in write_nonuniform_line_profile'
         endif
         if(nj_tab(i).lt.0) then
           write(*,*)'nj_tab(i).lt.ndensa'
           write(*,*)'i=',i,'nj_tab(i)=',nj_tab(i)
          write(*,*)'it should be nj_tab(i).ge.0'
          write(*,*)'Please change nj_tab(i)'
          stop 'in write_nonuniform_line_profile'
         endif
      enddo

      write(*,*)'nj_tab',nj_tab



      return
      end      !!!!!!!write_nonuniform_line_profile!!!!!!!!!!!!



      subroutine default_in
!     [end at line ~].
c--------------------------------------------------------------------------
c     It creates default input data in genray.dat (mixed units) file.
c     The results of the work  are in common block files: *nml.in
c-----------------------------------------------------------------------      
      implicit none

      include 'param.i'
c      include 'one.i'
c      include 'ions.i'
c      include 'three.i'
c      include 'five.i'
c      include 'cone.i'
c      include 'grill.i'
c      include 'onetwo.i'
c      include 'output.i'
c      include 'rkutta.i'
c      include 'six.i'
c      include 'scatnper.i'

      include 'one_nml.i'
      include 'ions_nml.i'
      include 'cone_nml.i'
      include 'dinit_nml.i'
      include 'emissa_nml.i'
      include 'grill_nml.i'
      include 'onetwo_nml.i'
      include 'output_nml.i'
      include 'rkutta.i'
      include 'six_nml.i'
      include 'scatnper_nml.i'
      include 'adj_nml.i'
cSAP090204
      include 'edge_prof_nml.i'
cSAP091225
      include 'lsc_approach_nml.i'
c-------------------------------------------------- 
c     input
c--------------------------------------------------
      integer i_genray_in_transformation

c      
c-----these two arrays are for namelist work only
c     namelist does not work at PC with the names te0 and teb
c      dimension ate0(nbulka),ateb(nbulka) !declared in one__nml.i
c-----tmprof is the working array for the namelist input of density 
c     and temperature
c     Following is working array for the namelist input of density 
c     and temperature,...
c--------------------------------------------------------
c      locals
c--------------------------------------------------------
      real*8 
     &rho,rmax,rmin
      integer 
     &i,k,j


c      namelist /genr/ r0x,b0,outdat,stat,mnemonic,rayop,dielectric_op,
c     +partner
c      namelist /tokamak/ indexrho,ipsi,ionetwo,ieffic,psifactr,
c     +deltripl,nloop,i_ripple,eqdskin,NR,
c     &n_wall,max_limiters,n_limimer,r_wall,z_wall,r_limiter,z_limiter,
c     &phi_limiter,h_add_wall 
c      namelist /wave/ frqncy,ioxm,ireflm,jwave,istart,delpwrmn,ibw,
cSAP090304
c     *no_reflection,
c     *i_vgr_ini,poldist_mx,ioxm_n_npar
c     &cnperp_plot_min,cnperp_plot_max,n_nperp_plot,
c     &cN_perp_root_max,n_points_root,
c     &i_look_roots,k_hot_root,
c     &i_rho_find_hot_nperp_root,
c     &rho_step_find_hot_nperp_roots,rho_min_find_hot_nperp_roots,
c     &shift_rho_geom
c      namelist /grill/ i_n_poloidal,n_theta_pol,ksi_nperp,
c     *i_rho_cutoff,,rho_step_find_LHFW_cutoff, 
c     *rho_initial_find_LHFW_cutoff, 
c     *ngrilld, !the old name for ngrill for old genray.in files
c     *ngrill,igrillpw,igrilltw,rhopsi0,thgrill,
c     * phigrill,height,nthin,anmin,anmax,nnkpar,powers,
c     &antormin,antormax,nnktor,anpolmin,anpolmax,nnkpol
c     &ilaunch,r0launch,z0launch,phi0launch,i_grill_pol_mesh
c     &i_grill_npar_ntor_npol_mesh:
c      namelist /rz_launch_grill_tab/prof_2d
c      namelist /scatnper/ iscat,scatd,rhoscat
c      namelist /dispers/ ib,id,iherm,iabsorp,iswitch,del_y,jy_d,
c     *idswitch,iabswitch,n_relt_harm,n_relt_intgr,iflux,
c     &i_im_nperp,i_geom_optic,ray_direction,errabs0,errrel0,navg,
c     &diff_err,relres,iabsorp_collisional,coll_mult,refl_loss,
c     &n_relt_harm,n_relt_harm1,i_salphal,ion_absorption,
c     &iabsorp_ql
c      namelist /numercl/ irkmeth,ndim1,isolv,idif,nrelt,
c     * prmt1,prmt2,prmt3,prmt4,prmt6,prmt9,icorrect,iout3d,
c     * maxsteps_rk,i_output,
c     & i_uh_switch,uh_switch,prmt6_uh_switch,
c     &toll_hamilt,
c     &i_resonance_curve_integration_method,epsi,
c     &eps_delta_pow
c      namelist /output/ iwcntr,iwopen,iwj,itools,i_plot_b,i_plot_d,
c     &n_plot_disp,r_plot_disp,id_plot_disp,
c     &s_poloid_plot_disp,point_plot_disp,     
c     &i_plot_disp_cold,
c     &n_plot_disp_cold,s_poloid_plot_disp_cold,r_plot_disp_cold,
c     &point_plot_disp_cold,     
c     &i_plot_wave_normal_cold,
c     &number_map_points_real_nperp,number_map_points_image_nperp,
c     &ratio_min_r_nperp,ratio_max_r_nperp,
c     &ratio_min_i_nperp,ratio_max_i_nperp,
c----frequencies plot along the straight line
c     &r_freq,z_freq,alpha_freq,beta_freq,dist_freq,max_plot_freq,
c     &nsteps_freq,n_ec_harmonics_freq,npar_freq
c      namelist /plasma/ nbulk,izeff,idens,temp_scale,den_scale,ndens,
c     &nonuniform_profile_mesh
c      namelist /species/ charge,dmas
c      namelist /varden/ var0,denn,denm,an,sigman
c      namelist /denprof/ dense0,denseb,rn1de,rn2de   
c      namelist /tpopprof/ tp0,tpb,rn1tp,rn2tp
c      namelist /vflprof/vfl0,vflb,rn1vfl,rn2vfl
c      namelist /tprof/ ate0,ateb,rn1te,rn2te
c      namelist /zprof/ zeff0,zeffb,rn1zeff,rn2zeff
c      namelist /dentab/ prof
c      namelist /temtab/ prof
c      namelist /tpoptab/ prof
c      namelist /vflowtab/ prof
c      namelist /zeftab/ zeff1
c      namelist /read_diskf/ i_diskf,
c    . netcdfnm,
c    . rtem0,
c     . rtail,r1t,r2t,ttail,rhot,r1h,r2h,thotpar,thotper,
c     . hotexp,hotmnpar,hotmxpar,hotmnper,hotmxper,
c     . rbeam,r1b,r2b,tbeam,ebeam,thbeam,
c     . jx,iym,lrz,ngen,
CENM 31Aug05 Added (optional) parameters at the end if dispers namelist
C    to be used in the relativistic dispersion relation in abhay_disp.f
c     & rvtovte1,rvtovte2,rtemp1,rtemp2,rtemp3
c      namelist /emission/i_emission,tol_emis,nharm1,nharm2,nfreq,
c     .freq00,ffreq01,wallr,i_rrind,i_r_2nd_harm,
c      i_emission_spectrum,jx_kin,max_kin_energy_kev
c      namelist /ox/ i_ox,
c       &theta_bot,theta_top,
c      i_ox_poloidal_max,eps_antenna,eps_xe


c     namelist /eccone/  NOT INCLUDED, but defaults set below.
c     NEED TO DIMENSION MRAY(),CR()   THIS needs regularization!, BH040412
c      character*8 raypatt !it is specified in cone.i
c      integer gzonemax
c      parameter (gzonemax=20) !set in param.i
c      integer gzone,nray_in,mray(gzonemax),cr(gzonemax) !in eccone.i
c     namelist /adj_nml/
c    & i_adj,
c    & t,ze,umax,dt,
c    & aerrmx,rerrmx,epsi_adj,
c    & npsi0,nthp0,nmax_chi,imax_chi,lmax_chi,tmax,alpha,rho_,
c    & n_harm_min,n_harm_max,
c    &n_relt_intgr_adj,i_resonance_curve_integration_method_adj,
c    & i_calculate_or_read_adj_function,
c    & i_chi_interpolation   
c     namelist/edge_prof_nml/ 
c     &i_edge_dens_anal,
c     &i_edge_dens_rz_mesh,
c     &n_pol_edge_dens,
cSAP090304
c     &sigmedgn,sigmedgt,
c     &theta_pol_edge_dens_ar_degree,
c     &sigmedgn_ar,
c     &dens_min_edge,temp_min_edge,
c     &theta_pol_edge_1_degree,    !for analytic formula of sigma_edge_n
c     &theta_pol_edge_2_degree,          
c     &sigma_theta_pol_edge_1_degree,
c     &sigma_theta_pol_edge_1_degree
c     &sigma_edgen_0,
c     &sigma_edgen_1,
c     &sigma_edgen_2,
c     &sigma_wall_n,
c     &sigma_lim_toroidal_degree,
c     &nxeqd_add,nyeqd_add
cSAP100406
c     namelist /lsc_approach_nml/  
c     &i_lsc_approach,
c     &vpar_dev_vte_max,  
c     &nv_lsc_ql_lh,
c     &n_psi_TSC,
c     &icount_iter_max, 
c     &eps_diff_d_q,
c     &n_power_steps,
c     &n_dim_radii_f_plots,
c     &n_radii_f_plots_ar
c&genr
      r0x=1.0d0 ![m]
      b0=1.0d0  ![Tl]
      outdat='zrn.dat'
      stat='new'
      mnemonic='genray'
      rayop='both'
      dielectric_op='disabled'
      partner='disabled'
c&end
      write(*,*)'in default_in after set data for /genr/'
      write(*,*)'r0x=',r0x
      write(*,*)'b0=',b0
      write(*,*)'outdat=',outdat
      write(*,*)'stat=',stat
      write(*,*)'mnemonic=',mnemonic
      write(*,*)'rayop=',rayop
      write(*,*)'dielectric_op=',dielectric_op
      write(*,*)'partner=',partner

!/genr/ namelist   	  (NSTX, FW,cold plasma,one ray)
!-------------------------------------------------------------------------
! mnemonic, is the run designator...to help keep track of runs.
!           It is used for naming the ray data output files:
!              mnemonic.txt and mnemonic.nc 
!              (Ray data o/p also depends on rayop nmlst variable.)
!            mnemonic is character*128, default="genray"
! rayop,     Specifies which of mnemonic.txt and mnemonic.nc files
!            are to be output:
!            "both", "text", "netcdf", or "none".  
!            rayop is character*8, default="both".
!            [Previous (related) out3d nml is no longer operative.]
! dielectric_op="enabled",adds output of the 9 complex dielectric tensor
!               elements to the .txt and .nc ray data files.
!               "disabled", omit such data from the data files.
!               dielectric_op is character*8, default="disabled"
!-------------------------------------------------------------------------
!Normalization constants:
! r0x (m) characteristic length, can be used as scale factor
! b0 (tl) characteristiic magnetic field, can be used as scale factor
!-------------------------------------------------------------------------
!Parameters for output files
!--------------------------------------------------------------------------
! outdat*20     name of output file
! stat*3        status of output file
!--------------------------------------------------------------------------
! partner  = 'disabled' to use input profiles from genray.dat or genray.in
!          = 'genray_profs_in.nc'  to use plasma profile data from the netCDF
!                       file: genray_profs_in.nc  created by transport code
!                       or otherwise;
!                  and to output power and current profiles on transport
!                      code radial grid to file: genray_profs_out.nc
!                      and genray_profs_out.txt (text file).
!          = 'genray_profs_in.txt'  to use plasma profile data from the text
!                      file: genray_profs_in.txt  created by transport code
!                      or otherwise;
!                  and to output power and current profiles on transport
!                      code radial grid to file: genray_profs_out.nc
!                      and genray_profs_out.txt (text file).
!------------------------------------------------------------------------
!/tokamak/
!-------------------------------------------------------------------------
!Tokamak
!--------------------------------------------------------------------------
! eqdskin=Name (character*512) of eqdsk equilibrium input file
!         "equilib.dat" (default)
!--------------------------------------------------------------------------
! Type of the radial coordinates
! indexrho  1 - sqrt(area), 2 - sqrt(torflux), 3 - sqrt(volume), 
!           4 - sqrt(psi-psimag), 5 - (psi-psimag)]
!           6 - (r_max(psi)-r_min(psi))/(r_max(psilim)-r_min(psilim))
! -------------------------------------------------
! ipsi=1 calculation of contours psi(z,r)=const
!     =0 -read in these contours from psi.bin file
! -------------------------------------------------
! ionetwo=1-calculation power and current radial
!           profiles, to the file onetwo.bin
!           0 - no calculations)
!--------------------------------------------------------------------------
! ieffic  choice of formula for the current drive efficiency
!        =1 asymptotic simple formula (homogeneous, nonrelativistic)
!        =2 asymptotic formula (East-Karney )
!        =3 asymptotic formula (curba subroutine)
!        =4 Lin-Liu (TorGA_curgap subroutine)
!        =5 using ADJ function for n_harmonic=0 case (LH wave)
!           it works for i_adj=1 case only
!        =6 using ADJ function for all harminics general case
!           it works for i_adj=1 case only
!--------------------------------------------------------------------------
! psifactr (it should be 0 < psifactr =<1,  psifcatr ~1)
!         is the parameter for the creation of the limiter points
!         using the closed flux surface:  psi(r,z)=psilim*psifactr 
!         psifactr is a parameter (it must be .le.1) to avoid
!         problems with the psi function near the separatrix.
!--------------------------------------------------------------------------
!deltripl is the relative amplitude of the ripple field at the
!         last flux surface (at rho=1)
!
!nloop    is number of toroidal field coils 
!
!i_ripple is the index to choose the ripple model
!         bripl_phi(z,r,phi)=(dF/dphi)/r=cos(N_loop*phi)*g(r,z)*N_loop/r
!         bripl_z(z,r,phi) =(dF/dz)    =sin(N_loop*phi)*(dg/dz)
!         bripl_r(z,r,phi) =(dF/dr)    =sin(N_loop*phi)*(dg/dr)
!         models for function g:	 
!         =1 the ripple model approximating the DIII-D field
!            g=beqd*reqd*deltripl*(r/rmax)**N_loop/N_loop 
!            beqd is the toroidal magnetic field at reqd (Tl)
!            reqd is the nominal major radius of the torus.
!            rmax is the max major radius at the last closed flux surface
!            r    is the major radius
!         =2 the ripple model using modified Bessel function I_0
!            g=beqd*reqd*deltripl*I_0(N_loop*rho(z,r))/(N_loop*I_(N_loop))
!--------------------------------------------------------------------------
! NR is the number of bin boundaries in the small radius direction
!    for the calculation of the power and current drive radial profiles.
!    Power and current is tabulated at (NR-1) bin centers.
!--------------------------------------------------------------------------
! wall and limiter coordinates
!             if n_wall= 0 then
!                no wall to be used, no reflection from the wall
!
!             if n_wall>0 then 
!                wall coordinates"  r_wall, z_wall,
!                and reflect rays from straight line segments between
!                the given points.  Count the number of coord pairs.
!                The coords must begin and end at the same physical
!                point.
!             It should be (n_wall .le. n_wall_a)
!             n_wall_a is a maximal value of n_wall, It is set in param.i file
!
! max_limiters  number of limiters
!               if max_limiters=0 no limiters will be used
!               It should be max_limiters.le.max_limiters_a
!
! n_limiter(1:max_limiters) is a number of limiter points
!             if n_limiter=0 the
!                 no limiter to be used, no reflection from the limiter
!             if n_limiter>0 then 
!                Reflect from the
!                chamber wall consisting of the wall coordinates and
!                those limiter coordinates to the plasma side of the
!                wall.
!             It should be (n_limiter .le. n_limiter_a)
!             n_limiter_a is a maximal value of n_limiter.
!             It is set in param.i file!
!
!  r_wall(n_wall),z_wall(n_wall) wall coordines [m]
!  r_limiter(n_limiter(i),max_limiters),
!  z_limiter(n_limiter(i),max_limiters) limiter coordines [m]
!
! phi_limiter(1:2,1:max_limiters) toroidal angles [degree] of the limiters
!                                  boundaries
! 0 =< phi_limiter(1,i) < phi_limiter(2,i) =< 360   
!
!      h_add_wall    the distance between wall points of the wall mesh
!                    with additional points. It is measured in m
!---------------------------------------------------------------------

c&tokamak
      eqdskin='equilib.dat'
      indexrho=2
      ipsi=1
      ionetwo=1
!      ieffic=2
      ieffic=6
!cSAP09116
!      psifactr=0.97d0
      psifactr=0.999d0
      deltripl=0.00d0
      nloop=24
      i_ripple=1
      NR=NRA
      n_wall=0
      do i=1,n_wall
        r_wall=0.d0
        z_wall=0.d0
      enddo

      max_limiters=0

      do i=1,max_limiters_a
        n_limiter(i)=0.d0
      enddo

      do k=1,max_limiters_a
         phi_limiter(1,k)=0.d0
         phi_limiter(2,k)=360.d0 
         do i=1,n_limiter_a
            r_limiter(i,k)=0.d0
            z_limiter(i,k)=0.d0
         enddo
      enddo
      h_add_wall=1.d-3
c&end

!/wave/
!-------------------------------------------------------------------------
!Waves
!-------------------------------------------------------------------------
! frqncy frequency f=w/2pi in GHz
!
! ioxm ( 1 - om, -1  - xm )   ! Wave mode in the equation N_perp=N_perp(gam).
!                             ! gam is the angle between the refractive vector
!                             ! and the magnetic field. It works if ioxm_n_npar=0. 
! For ioxm_n_npar=+1 or -1    ! dispersion equation: a*N**4+b*N*2+c=0 
! the wave mode will be       ! roots: N**2=(-b+ioxm*sqrt(b**4-4a*c))/2a 
! specified by  ioxm_n_npar   ! Here coefficients (A,B,C) are the functions
!                             ! of angle (gam)  
!                             ! a=A/delta**3, b=B/delta**3, c=C/delta**3
!                             ! delta=1-y_e for ib=1
!                             ! delta=1-y_i for ib=1 > 1 
!-----------------------------------------------------------------------
! 
! ioxm_n_npar - 
!           =0 (as default)ioxm_n_par will not be used,
!              the wave mode will be calculated using ioxm parameter          
! 
!           (=+1 or =-1) sign before square root in dispersion relation
!           delta**2*f~*N**4+delta*g~*N**2+w~=0
!           which gives N**2=N**2(N_parallel)
! 
!           f~=delta*eps_per
!     
!           g~=delta**2*N_par**2(eps_par-eps_per)+
!              delta**2(g**2-eps_per**2-eps_per*eps_par)
!
!           w~=delta**3*N_par**2(-eps_per*eps_par+eps_per**2-g**2)+
!              delta**3*eps_par(eps_per**2-g**2)
!
!           root=(-g+ioxm_n_npar*sqrt(g**2-4f*w))/2g
!
!           f=f~, g=g~/delta, w=w~/delta**2
!           delta=1-y_e for ib=1
!           delta=1-y_i for ib=1 > 1 
!
!------------------------------------------------------------------------- 
! ireflm  -max number of reflections =1 for EC
!-------------------------------------------------------------------------
! no_reflection !=1 switch off the artificial reflection from 
!               !the last closed flux surface.  Gives natural reflection
!               !from a density gradient outside the LCFS.
!               !=0 (default) switch on the artificial reflection from 
!               !the last closed flux surface.
!-------------------------------------------------------------------------
! jwave  (0 - LH wave, -1 AW, 1 - EC wave) wave harmonic used in calc.
!    of current drive efficiency (see ieffic).
! -------------------------------------------------
! istart  if start point outside the plasma=1 else=2
! if istart=1 use namelist &eccone below, =2 use &grill
! if istart=3 it use &grill and the additional calculations in dinit
! to launch the ECR ray inside the plasma in the O_X mode
! conversion point (rhoconv,theta), Theta is a poloidal angle (degree)
! for mode conversion point. It is given in dinit.f 
!--------------------------------------------------------------------------
! delpwrmn - Minimum power in each ray, as a fraction of
!            starting power in the ray, after which ray is stopped.
!--------------------------------------------------------------------------
! ibw=0 it is not the direct launch of the Bernstein waves
!    =1 the direct launch of electron Bernstein wave from dhot tensor
!       The last case works only for istart=2 and grill_lh conditions 
!--------------------------------------------------------------------------
! i_vgr_ini =+1 the wave is directed into the plasma (in the initial point)
!           =-1 the wave is directed out the plasma (in the initial point
!---------------------------------------------------------------------------
! poldist_mx is the maximal poloidal distance (m) along the ray
!            default=1.d+5
!------------------------------------------------------------------------
! i_look_roots=0    !do not plot D(N_perp) and do not calculate all 
!                   !hot roots
!             =1    !plot D(N_perp) and calculate all 
!                   !hot roots, but do not calculate ray   
!             =2    !calculate hot roots, use the root with number k_root
!                   !as the initial ray condition and calculate a ray
!----------------------------------------------------------------------       
! cnperp_plot_min,cnperp_plot_max !max and min Nperp to plot D(Nperp)      
! n_nperp_plot,                   !number of Nperp points to plot D(Nperp)
!-------------------------------------------------------------------------     
! cN_perp_root_max               !max value of n_perp to
!                                !find hot roots 
! n_points_root                  !number of  N_perp mesh points
!                                !to find hot plasma roots
!-------------------------------------------------------------------------
! k_hot_root   is the number of the hot plasma root
!              N_perp_root_ar(k_hot_root)
!              which will be used for ray initial condition
!              It works for i_look_roots=2 case only
!-------------------------------------------------------------------------
! i_rho_find_hot_nperp_roots=1  find the small radius rho_ini 
!             rho_ini > rho_min_find_hot_nperp_roots
!             at the vector rho^ where
!             hot plasma dispersion function D_hot(nper)=0
!             has one,two or three roots.
!             The vector rho^ is starting at the edge point 
!             (r_edge,z_edge,phi_edge),and directed to the 
!             magnetic axis O(xma,yma,phi_edge)
!             Write roots and polarization to  find_hot_roots.dat
!i_rho_find_hot_nperp_roots=0  do not find roots
!
! rho_step_find_hot_nperp_roots is the small radius step to find the hot
!                               plasma dispersion relation D_hot(N_perp)=0
!                               roots 
!
! rho_min_find_hot_nperp_roots  is the minimal rho
!-----------------------------------------------------------------------
!     For istart=1 case only
!     shift_rho_geom 
!                    It shifts the initial point of the ray from the plasma
!                    bounday (LCFS) inside the plasma.
!                    The ray initial coordinates at the LCFS are (r,z)
!                    and the geometrical radius of this initial point
!                    (the distance from the magnetic axis) is
!                    rho_geom=sqrt((r_xma)**2+(z-yma)**2)
!                    
!                    The initial point will be shifted inside the plasma
!                    rho_geom_shifted==rho_geom*shift_rho_geom
!                    So,the shifted initial coordinates are  
!                    r_shifted=xma + (r-xma)*shift_rho_geom
!                    z_shifted=yma + (z-xma)*shift_rho_geom
!
!                    By default shift_rho_geom=1.0d - 1.d-7
!-----------------------------------------------------------------------
c&wave
      frqncy=60.00d-3 !GHZ
      ioxm=-1
      ireflm=3
      no_reflection=0
cSAP090304
c      sigmedgn =0.02   moved to/ edge_prof_nml/ 
c      sigmedgt =0.02   moved to/ edge_prof_nml/ 
      jwave=1
      istart=2
      delpwrmn=1.d-2
      ibw=0
      i_vgr_ini=+1
      poldist_mx=1.d+5 ![m]
      ioxm_n_npar=0
      i_look_roots=0   
      cnperp_plot_min=0.d0
      cnperp_plot_max=5.d0      
      n_nperp_plot=50
      cN_perp_root_max=5.d0 
      n_points_root=50   
      k_hot_root=1 
      i_rho_find_hot_nperp_roots=0
      rho_step_find_hot_nperp_roots=1.d-2
      rho_min_find_hot_nperp_roots=0.9d0
      shift_rho_geom=1.d0
c&end 
!/scatnper/
!-------------------------------------------------------------------------
!N_perpendicular scattering
!-------------------------------------------------------------------------
! iscat it is the switch for the n_perp scattering
!       iscat=1 the scattering swithed on,
!            =0 the scattering swithed off
!-------------------------------------------------
! rhoscat(1:nscat_n) normalized small radii for the scattering location
!        The parameter nscat_n should be given in the param.i file
!
! The scattering of the polar angle deltheta will be
!      deltheta=dsqrt(2.d0*scatd)*ranorm(fseed)
! scatd(0) the mean square scattering angle (radians**2)
!          for the plasma boundary reflection points
!
! scatd(1:nscat_n) the mean square scattering angles (radians**2)
!          for the interior plasma boundary points
!-------------------------------------------------------------------
c&scatnper
      iscat=0
c-----here nscat_n=10
c-----scatd= 0.01, 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.01, 0.01
      do i=1,nscat_n
         scatd(i)=0.d0
      enddo
      scatd(1)=0.01d0
      scatd(nscat_n)=0.01d0
      scatd(nscat_n-1)=0.01d0   
c-----rhoscat= 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95, 0.97
      do i=1,nscat_n
         rhoscat(i)=0.1d0*i
      enddo
      rhoscat(nscat_n)=0.97d0
      rhoscat(nscat_n)=0.95d0
c&end

!/dispers/
!Dispersion relation
!-------------------------------------------------------------------------
! ib<=nbulk cyclotron resonance sort(=1 for ecr)
! the number in (1-y(ib)) for the multiplication of the
! dispersion relation to delete the singularity
! -------------------------------------------------
! id gives form of the dispersion relation
!            =1 AN**4+BN**2+C=0
!            =2 N**2=(-B+ioxm*Sqrt(B**2-4AC))/2A;
!	     =3 Appleton-Hartree;
!            =4 electron relativistic plasma from hermitian Mazzucato code
!            =5 electron relativistic plasma from total Mazzucato code
!            =6 hot non-relativistic plasma 
!            =7 electron relativistic plasma from Shkarofsky code
!            =8 Ono dispersion for fast waves
!            =9 hot non-relativistic plasma, full tensor [Bernstein-Friedland]
!            =10 Westerhof-Tokman dispersion with Mazzucato
!                relativistic electron dielectric tensor
!                Self-consistent absorption obtained with iabsorp=1.
!            =11 Eric Nelson-Melby relativistic tensor
!                Dispersion function = Re(Det eps^Hermitian, if iherm=1)
!		 calculated using Isaac Weiss' method
!            =12 Westerhof-Tokman dispersion with 
!                    Nelson-Melby dielectric tensor
!                Self-consistent absorption obtained with iabsorp=12.
!            =13 Westerhof-Tokman dispersion with
!                 hot non-relativistic plasma, full tensor
!            =14 Abhay Ram's dielectric tensor, using the Trubnikov integral
!	     =15 Westerhof-Tokman dispersion with
!	           Abhay Ram's dielectric tensor
!                Self-consistent absorption obtained with iabsorp=12.
!    **** NOTE: see notes around irkmeth if using id=11,12,14 or 15 *****
!              **** NOTE: irkmeth=1 should be used for id=14
!---------------------------------------------------------
! For use with Abhay Ram's dispersion relation (id=14 or id=15) parameters to
! control the integration routine for the Trubnikov integral.
!
! 'relres' offers 4 choices for the resolution to use for the relativistic
! dispersion relation using the Trubnikov integral:
! relres=1 low resolution, errabs0=1.d-4, errrel0=1.d-4, navg=3, diff_err=0.1
!       =2 medium res., errabs0=1.d-5, errrel0=1.d-5, navg=12, diff_err=1.d-3
!       =3 high res., errabs0=1.d-6, errrel0=1.d-6, navg=25, diff_err=1.d-6
!       =4 user-defined res., set the following parameters manually.
! default: relres=2
!
! The Trubnikov one-dimensional (complex) integral is performed by splitting
! up a region from 0 to 1.d8 into 10^6 pieces, and each piece is integrated
! using the SLATEC adaptive quadrature routine dqag. errabs0 and errrel0 are
! the absolute and relative error tolerances passed directly to dqaq.
! Then the adjacent pieces are compared (it is an oscillatory integrand)
! and using navg number of pieces, when the average difference between them
! are less then diff_err, the integration is presumed finished (Thus it may
! finish long before the upper limit of 1.d8).
!
! errabs0 - absolute error for dqag integration routine
! errrel0 - relative error for dqag integration routine
! navg - number of adjacent integration intervals to use in comparison
! diff_err - error tolerance using navg pieces, when the average difference
!         is less than diff_err, then the integration is done.
!
! To decide when one should use the low, medium, or high resolution 
! integration, here are some suggestions based on the behavior of the
! Trubnikov integrand: The integrand converges more slowly, and hence
! the resolutions should be set higher, for low electron temperature,
! low (i.e. near zero) magnitude of n_parallel, and for low (near or
! below the fundamental cyclotron frequency) frequency. 
! Examples: n_parallel = -0.05, Te=400 eV, omega/omega_ce=0.4 to 1.2, 
!    it was necessary to use errabs0=1.d-5,errrel0=1.d-5,navg=20,diff_err=1.d-5
!    to be completely converged.  By changing Te to 4000 eV, it was sufficient
!    to use 1.d-4,1.d-4,15 and 1.d-4. 
!   An easy case: n_parallel=0.3, Te=7 keV, omega/omega_ce=2.4 to 2.7, 
!    complete convergence already at errabs0=1.d-4,errrel0=1.d-4,navg=2,
!    diff_err=0.5
!   An intermediate case: n_parallel=0.1, omega/omega_ce=1.0, Te=300 eV
!    errabs0=1.d-5,errrel0=1.d-5,navg=12,diff_err=1.d-3 was OK.
! *** NOTE: Sometimes with too small n_parallel, the Trubnikov method does
!   not work well (id=14 or 15). Instead, use the Weiss method (id=11 or 12),
!  which works well for small n_parallel (but does not work for n_parallel>1).
!--------------------------------------------------------
! For Mazzucato plasma dispersion tensor:
! iherm =1 hermitian dielectric tensor, 2-full
!------------------------------------------------------------------------
!Absorption:
!iabsorp -choice of Imag(N_perp) (N_perp is the perpendicular refractive index)
!-------------------------------------------------------------------------
! iabsorp=1 for EC waves with Mazzucato solver (ki/kr<<1)
!        =2 for LH waves
!        =3 for FW waves, Chiu et al themal corr., NF 1989, with corrections.
!           At this time, this is only iabsorp value which provides power
!           profiles to individual ions (via powden_s/powtot_s).
!           Other models below (9, 91, 92, could be added).
!        =4 for all frequencies with Forest code (ki/kr<<1)
!        =5 for EC and EBW waves from Shkarofsky code
!        =6 for EC and BW anti-hermitian part relativistic tensor+
!                         hermitian_part (Forest code)
!        =7 for EC wave case.The complex electric field calculations
!           using Cold plasma tensor +antihermitian relativistic tensor
!           ---EC relativistic absorption 
!           dielectric tensor=hermitian part(cold plasma)+
!                             anti-hermitian part(full relativistic) 
!        =8 for id=10,12,13 and Westerhof-Tokman dispersion 
!           function  [Projection method from real part of
!           eigenvalue.]  Integration is with respect to distance.
!        =9 The absorption is calculated for hot dispersion
!           using the formula from Stix book p.74 (17,18,21)
!           Im(k_perp)= 0.5*Power_abs/(P^+T^)
!           It uses hot plasma dielectric tensor.
!           It calculates hot plasma dielectric tensor reps() and 
!           electric field polarization (cex,cey,cez) using this 
!           hot plasma tensor.
!        =91 FW absorption is calculated by
!            subroutine absorpfw_pinsker_1
!            using formula:               	                            
!            1) for ion absorption from                      
!            R.i.Pinsker,M.Porkolab                          
!            e-mail: pinsker@fusion.gat.com                  
!            05/02/11                                         
!            2) for electron absorption from                 
!            S.C.Chiu,V.S.Chan,R.W.Harvey,M.Porkolab	 
!            Theory of fast wave current drive for tokamak   
!            plasma,Nuclear fusion,1989,vol.29,No.12, p.2175-2186
!   
!            electric field polarization will be for the cold plasma dielectric tensor
!
!        =92 FW absorption is calculated by
!            subroutine absorpfw_pinsker_2 using formula:   
!            for electron and ion absorption from                      
!            R.i.Pinsker,M.Porkolab                          
!            e-mail: pinsker@fusion.gat.com                  
!            05/02/11                                         
!          
!            electric field polarization will be for the cold plasma dielectric tensor
!
!        =10 The absorption is calculated for relativistic tensor
!           (A.Ram id=14 or Nelson-Melby id=11)
!           using the formula from Stix book p.74 (17,18,21)
!           Im(k_perp)= 0.5*Power_abs/(P^+T^)
!           It uses relativistic dielectric tensor.
!           It calculates relativistic dielectric tensor reps() and 
!           electric field polarization (cex,cey,cez) using this 
!           tensor.
!        =11 The absorption Im(N_perpendicular) is calculated 
!           for the relativistic dispersion function Determinant=0.
!           The dispersion function uses the full Abhay Ram dielectric
!           tensor (with the Trubnikov integral) (for id=14) or
!	    Nelson-Melby with Weiss method (for id=11).
!           The projection method is used to find Im(N_perpendicular).
!           That is, Im(N_perp)=-Im(D)/(dD/dRe(N_perp)), evaluated 
!           at give Real(N_perp) with Im(N_perp)=0.
!           The real part Real(n_perp) is obtained from the trajectory
!           solution for the given id value. 
!           Then complex n_perp is used with full dielectric tensor
!           to get the polarizations.  Good for arbitrary id value,
!           like other iabsorp calcs.
!	 =12 Find Im(N_perpendicular) by finding the exact solution to
!            Det(Complex n_perp)=0. This returns Im(N_perp), like the
!            projection method described above, but is more accurate,
!            especially when Im(N_perp)/Re(N_perp) is not negligible.
!            Uses Muller algorithm to solve for the complex root in
!            the complex plane. See mullerfun2.f
!            For id.eq.11.or.id.eq.12, uses Eric Nelson-Melby tensor.
!            For id.eq.14.or.id.eq.15, uses Abhay Ram tensor.
!------------------------------------------------------------------------
!iabsorp_ql =0  do not use QL flux for absorption calculations
!           =1  to use QL flux for electron absorption calculations
!               In this case QL flux will be calculated for harmonics
!               numbers nharm in the following interval:
!               n_harm_adj_min =< nharm   =< n_harm_adj_max 
!               n_harm_adj_min   number of minimal and maximal harmonics
!               n_harm_adj_max   for power and CD calculations  
!               Electric field polarization will be calculated according to
!               the value of the index: iabsorp 
!               Energy flux "fluxn" will be calculated according to 
!               the value of the index: iflux 
!------------------------------------------------------------------------
!                 To switch off ion absorption
! ion_absorption ='enabled' to add ion absorption (by deafault)
!                  It works at iabsorp=3,9,91,92
!                ='disabled' do not add ion absorption
!------------------------------------------------------------------------
! iabsorp_collisional =0 no additional collisional absorption
!                     =1 collisional absorption  using formula
!                        Im(N)=dabs(nu_ei/(gr_perp))*clight/omega)
! coll_mult =1.d0(default), multiplies above coll absorp expression
!------------------------------------------------------------------------
! The change of the dispersion relation and absorption
! near the gyro-frequency points
!-------------------------------------------------
! iswitch=1   To use the change of the dispersion relation and
!             absorption
!        =0   Do not use the change of the dispersion relation and
!             absorption 
!     del_y   If the difference |1-nY(jy)|<del_y 
!             (jy=1-nbulk ,n=...-2,-1,0, 1,2,3,...)
!             then switch on the new 
!             given type of the dispersion and absorption.
!   jy_d      is the type of plasma species 1<=jy<=nbulk
!   idswitch  is the type of the dispersion function near the 
!             gyro-frequency points
!             It can be equal 1,2,3,4,5,6
!   iabswitch is the type of the absorption near the gyro-frequency point  
!----------------------------------------------------------------------- 
!   n_relt_harm1 is the lowest, i.e., minimum harmonic used in the
!               anti-hermitian dielectric tensor calculations.
!               It can be positive of negative.
!               Default value is +9999, in which case this input
!               is ignored.
!   n_relt_harm (.ge.1) gives the number of EC harmonics used 
!               in anti-hermitian dielectric tensor calculations
!               If n_relt_harm1=9999, then harmonics from 
!                  -n_relt_harm to +n_relt_harm are used.
!               If n_relt_harm1.ne.9999, then harmonics from
!                  n_relt_harm1 to n_relt_harm1+n_relt_harm are used.
!   It is necessary that the harmonics used in this calculation
!     be within the range of parameters [n_relt_harm1a,n_relt_harm2a]
!     set in the param.i file.s
!     These conditions are checked in the code.
!-------------------------------------------------------------------                
!   n_relt_intgr is the number of points for integration along the
!     resonance curve  (default=50).  Note, this variable is used
!     below with namelist i_resonance_curve_integration_method (in
!     numercl).
!---------------------------------------------------------------------
!  flux=B~.B+E~.d(omega*eps_herm)/(domega).E
!   iflux=1 the flux will be calculated using the the group velocity from
!           the chosen dispersion relation (with given id) and the electric
!           field calculated for the chosen iabsorp
!   iflux=2 the flux will be calculated using V_gr for the electron cold plasma
!           dispersion and polarization (using subroutine  grpde2)  
!-------------------------------------------------W----------------------
! i_im_nperp choice of the method to find Im_N_perp 
!    for hot plasma(iabsorp=4):
! i_im_nperp=1 Im_N_perp=abs(ImD_full/(dD_hermitian/dReN_perp)) 
!              (This method has been found to give poor accuracy
!               for FW in a DIII-D FW situation, see CompX
!               report CompX-2005-1.)
! i_im_nperp=2 (Re_N_perp,Im_N_perp) is the complex root 
!              (of the complex dispersion relation)
!              calculated by Newton iterations with the numerical
!              derivatives (the chord method)
!------------------------------------------------------------------
! i_geom_optic sets  the form of the ray equations
!              =1  integration in time (default):
!                  ray-tracing equations right hand side=
!		   dr^/dt=-(dD/dN^)/(dD/domega)
!                  dN^/dt=(dD/dr^)/(dD/domega)
!                  In this case rside1 gives v_group
!              =2  integration is space,
!                  ray-tracing equations right hand side=
!		   dr^/dl=- ray_direction * (dD/dN^)p
!                  dN^/dl=  ray_direction * (dD/dr^)p
!                  p=1.d0/dsqrt(deru(1)**2+deru(2)**2+(r*deru(3))**2)
!                  deru(1)=dD/dN_z, deru(2)=dD/dN_r,deru(3)=dD/dCM,
!                  N_phi=cm/r
!----------------------------------------------------------------------
! ray_direction =+1 as default it
!               -1 !Only for i_geom_optic=2
! It is a multiplier in right hand side of ray-tracing equations
! It is used for i_geom_optic=2 case
!----------------------------------------------------------------------
! i_salphal(nbulka)  sets which damping will be in salphal_nc
!                   for iabsorp=3 or for iabsorp=9 cases.
!                   For other 'iabsorp' cases 'saplhal' contains 
!                   the electron damping coefficients
!
!	     Default:i_salphal(1)=1,i_salphal(2:nbulk)=0 electron damping only  
!            A particular species contribution to salphal_nc is added if
!            i_salphal(species_number)=1. That is, damping coefficients
!            for all species with i_salphal(k).ne.0 are summed into saplhal_nc 
!			  
!----------------------------------------------------------------------------
! refl_loss fraction of power lost at each reflection
!----------------------------------------------------------------------------
c&dispers
      ib=2
      id=2
      relres=2
      errabs0=1.d-5
      errrel0=1.d-5
      navg=5
      diff_err=0.01
      iherm=1
      iabsorp=2
      iabsorp_ql =0  
      iabsorp_collisional=0
      coll_mult=1.d0
      iswitch=0
      del_y=1.d-2
      jy_d=2
      idswitch=2
      iabswitch=2
      n_relt_harm1=9999
      n_relt_harm=1
      n_relt_intgr=50
      iflux=1
      i_im_nperp=1
      i_geom_optic=1      
      ray_direction=1.d0
     
      do i=1,nbulka
        i_salphal(i)=0
      enddo
      
      i_salphal(1)=1
      refl_loss=0.0d0      
      ion_absorption ='enabled'
c&end

!/numercl/
!------------------------------------------------------------------------
!Numerical method
!-------------------------------------------------------------------------
! irkmeth (0-constant 1-variable step 5 order,2- variable step in RK 4 order)
!  irkmeth=0: Poloidal distance of output is at intervals .ge.prmt6.
!             Checks time step for passing outside plasma and reflects.
!  irkmeth=1: Only poloidal distance for control of output point (prmt6,
!             i_output has no effect). Output at distance.ge.prmt6,
!             i.e, the first code step beyond prmt6 distance.
!             No control for being outside the plasma and reducing
!             the step.  Correction method specified by icorrect is
!             operative.
!             Time or length for integration according to i_geom_optic.
!     *** NOTE: irkmeth=1 may not work well with fully relativistic dispersion
!        relations (id=11,12,14,15), unless prmt4 is quite small (e.g. 2.0d-6) ***
!  irkmeth=2: Most developed method of ray equation integration.
!             Time or space step in the of the equations
!             (according to setting of i_geom_optic) is controlled
!             so that output is at intervals prmt6 (meters).
!             As ray approaches the plasma edge, it is reflected
!             at the last closed flux surface.
!     *** NOTE: irkmeth=2 works best for fully relativistic dispersion
!        relations (id=11,12,14,15). An example of what prmt parameters work
!        well for relativistic EBW and irkmeth=2:
!                prmt1=0.000d+00
!                prmt2=9.999d+05
!                prmt3=1.000d-04
!                prmt4=5.000d-04
!                prmt6=1.0d-03
!---------------------------------------------------------------------------              
! ndim1 (number of the ray tracing equations)
! isolv=1 correction,=2 expl.solution
! idif=1 analytic differentiation, =2 numerical
! -------------------------------------------------
! nrelt   Maximum number of ray elements per ray.
!         Must be .le. nrelta (a parameter)
!--------------------------------------------------------------------------
! -------------------------------------------------
! Runge-Kutta method parameters
! -------------------------------------------------
! prmt1=tau initial= prmt(1)
! prmt2=tau final=prmt(2)
! prmt3=initial tau step=prmt(3) non-dimensional
! prmt4=required accuracy=prmt(4)
! prmt6=hprint=prmt(6) [meter] poloidal or total distance for results output
! prmt9 accuracy of Hamiltonian in the Runge-Kutta subroutine at irkmeth=2
!       it will reduce Runge-Kutta time step if dabs(ham).ge.eps_ham)
!--------------------------------------------------------------------------
!icorrect= switch for Hamiltonian correction in subroutine outpt
!          [See manual].
!icorrect=0 switch off the correction
!        =1 switch on the correction
!-------------------------------------------------------------------------
! iout3d [obsolete] ='enable'  to write output 3d.dat wile   [DEFUNCT]
!           ='disable' not write 3d.dat file         [DEFUNCT]
!--------------------------------------------------------------------------
! maxsteps_rk the maximal number of the time steps of the Runge-Kutta
!             solver (in default =10000)
!--------------------------------------------------------------------------
! i_output is used for irkmeth=2 only
! i_output=1 output is at the equal poloidal distance prmt6
!         =2 output is at the equal total distance prmt6
!--------------------------------------------------------------------------
!
! The following has been used for OXB in cases where the UH layer
!   is very close to the plasma boundary.   Then, in the vicinity of
!   the UH layer, switch the step size along the ray to a shorter
!   value.
! i_uh_switch=1    if uh=dsrt(xe+ye**2) < uh_switch then change 
!                  the output step prmt(6) to prmt6_uh_switch 
!            =0    do not change the output step prmt(6)
!
! prmt6_uh_switch  [meter] is the output step for i_uh_switch=1 case
!
! uh_switch       if uh<uh_switch then change the output step for
!                 i_uh_switch=1 case
!--------------------------------------------------------------------------
! Measure error in the dispersion relation.
! If  toll_hamilt <D/(N|gradD|) then stop ray calculation
!-------------------------------------------------------------------------
!    
! i_power_switch_resonance   =1  to use  prmt6_power_switch_resonance
!                            =0  do not change the output step prmt(6)
!
! prmt6_power_switch_resonance   is the output step for 
!                                i_power_switch_resonance=1 case
!                                 in resonace area
!
! n_power_switch_resonance   is the number of different used EC resonances 
!
! y_power_switch_resonance(n_resonance_a) are used ratios omega_ce/omega
!
! del_y_power_switch_resonance determines the resonance area
! The condition for resonance area 
! abs(Ye-y_power_switch_resonance(k))< del_y_power_switch_resonance
! k=1,...,n_power_switch_resonance 
! n_power_switch_resonance_a is max of  n_power_switch_resonance
!                            It is set in param.i
!-------------------------------------------------------------------
!       i_resonance_curve_integration_method=1 !rectangle integration
!                                              !over angle,
!                                              !for ellipse case only
!       i_resonance_curve_integration_method=2 !rectangle formula,
!                                              !p_perp integration
!                                              !elliptic or hyperbolic case
!       i_resonance_curve_integration_method=3 !trapezoidal formula,
!                                              !p_perp integration
!                                              !elliptic or hyperbolic case
!       i_resonance_curve_integration_method=4 !adaptive Simpson integration 
!                                              !p_perp integration
!                                              !elliptic or hyperbolic case
!
!       i_resonance_curve_integration_method is used in subroutine intgr_rl
!       to choose the numerical integration method for 
!       anti-hermitian relativistic dielectric tensor calculation.
!       This applies for iabsorp=6,7 and for emission calculations.
!       n_relt_intgr (from namelist &disper) is number of number of points 
!         in integration for i_resonance_curve_integration_method=1,2,3.
!-----------------------------------------------------------------------
!       epsi  absolute accuracy used in adaptive Simpson 
!-----------------------------------------------------------------------
!       eps_delta_pow   ! If the power variation along a ray at one output step
!                       ! bigger than eps_delta_pow:
!                       ! (delpwr(is-1)-delpwr(is))/delpwr(is-1) > eps_delta_pow
!                       ! and if the ray power is bigger than initial ray power
!                       ! multiplied by eps_delta_pow:
!                       ! delpwr(is-1)> delpwr(1)*eps_delta_pow 
!                       ! then the code will reduce the output step prmt(6).
!                       !
!                       ! Then the code will return the output step prmt(6)
!                       ! to the initial value prmt6 given in the input file.
!                       !
!                       ! eps_delta_pow is used in prep3d.f
!--------------------------------------------------------------------------
c&numercl
      irkmeth=2
      ndim1=6
      isolv=1
      idif=1
      nrelt=5000
      prmt1=0.000d+00
      prmt2=9.999d+05
      prmt3=2.000d-02
      prmt4=1.000d-02 
      prmt6=2.000d-02 ![m]
      prmt9=1.d15
      icorrect=1
      iout3d='enable'  !  [DEFUNCT]
      maxsteps_rk=10000
      i_output=1
      i_uh_switch=0
      uh_switch=1.5d0 
      prmt6_uh_switch=1.d-5
      toll_hamilt=1.d-3
c
      i_power_switch_resonance   =0 
      prmt6_power_switch_resonance=prmt6*1.d-1
      n_power_switch_resonance=1
      y_power_switch_resonance(1)=0.5d0
      del_y_power_switch_resonance=1.d-2 

      i_resonance_curve_integration_method=4
      epsi=1.d-5
c&end
      prmt(1)=prmt1  
      prmt(2)=prmt2  
      prmt(3)=prmt3   
      prmt(4)=prmt4   
      prmt6=prmt6/r0x  ! normalization of the output step
      prmt(6)=prmt6   
      prmt(9)=prmt9
c      ihlf=prmt9
c      prmt(9)=ihlf
      eps_delta_pow=1.d+10 
!/output/
!----------------------------------------------------------------------
! iwcntr =1 genray.f will calculate the contours of the
!            gyrofrequency omega_c=n at the poloidal cros-section (r,z) plane
!        =0 genray.f will not do it
! iwopen =1 mk_grapc will calculate open contours omega_c_iwj=n  (using contrb1)
!         2 mk_grapc will calculate close contours omega_c_iwj=n (using contrb2)
! iwj     1 <= iwj <= nbulk the number of plasma component
!          mk_grapc will calculate contours omega_c_iwj=n, iwj is a kind of the plasma
!         component must be.le.nbulk, iwj=1 for the electron gyrofrequency
!         iwj.ge.2  for the ion (iwj kind) gyrofrequency
! itools =0 do not use mkgrtool
!        =1 to use mkgrtool
!----------------------------------------------------------------------
!
! i_plot_b =1 create figures for the magnetic field,density and temperature 
!             profiles in plot.ps file using subroutine map_b based on PGplot
!             Also, plot characteristic frequencies to *.bin files.
! i_plot_b =0 do not write the b,n,T figures to plot.ps file 
!
!-----------------------------------------------------------------------
!  to plot dispersion function contours D_cold(ReN_perp,ImN_perp) at given points
!  for hot plasma absorption: iabsorp.eq.4 and iabsorp.eq.6 
!-------------------------------------------------------------------------
! i_plot_d =1 create the dispersion function contours d(ReNperp,ImN_perp)
!             in plot.ps file using PGplot
!             The code will plot these contours in all output points
!             in prep3d.f specified by output step prmt6
! i_plot_d =0 do not create the dispersion function contours in prep3d
!             
!---------------------------------------------------------------------------
! for plotting dispersion function contours D(ImN_perp,ReN_perp) at specified
! poloidal lengths or major radii in plot.ps file using PGplot
!---------------------------------------------------------------------------
!
! n_plot_disp=0 do not plot contours D(ReN_perp,Im_N_perp)
!          0< n_plot_disp=<n_plot_dispa is the number of major radius
!          points where contours will be ploted
! id_plot_disp  determines the dispersion function D type
!          used for contours plots
! r_plot_disp(n_plot_disp) major radiusl [m] where contours will be plotted
!
! s_poloid_plot_disp(n_plot_disp) poloidal distance [m] where contours
!                                 will be plotted
!
! point_plot_disp ='poloidl_dist' to create D contours at given
!                   s_poloid_plot_disp() as default
!                 ='major_radius'  to create D contours at given
!                   r_plot_disp()
!
! number_map_points_real_nperp  is the number of map points
!                               in Real(N_perp) direction
!
! number_map_points_image_nperp  is the number of map points
!                                in Image(N_perp) direction
!
! ratio_min_r_nperp,ratio_max_r_nperp  set the ratio of
!          minimal and maximal map boundaries in Real N_perp direction
!          (min_r_nperp < Real(N_perp) < max_r_nperp) to the value of
!          Real(N_perp_ray_=cnper along the ray:
!          min_r_nperp= Real(N_perp_ray)*ratio_min_r_nperp
!          max_r_nperp= Real(N_perp_ray)*ratio_max_r_nperp
!          These parameters should be: 
!           0 =< ratio_min_r_nperp < 1
!           1 <  ratio_max_r_nperp 
!
! ratio_min_i_nperp,ratio_max_i_nperp  set the ratio of
!          minimal and maximal map boundaries in Image N_perp direction
!          (min_i_nperp < Image(N_perp) < max_i_nperp) to the value of
!          Image(N_perp_ray)=cnprim along the ray:
!          min_i_nperp= Image(N_perp_ray)*ratio_min_i_nperp
!          max_i_nperp= Image(N_perp_ray)*ratio_max_i_nperp
!          These parameters should be: 
!           0 =< ratio_min_i_nperp < 1
!           1 <  ratio_max_i_nperp 
!
!          If Image(N_perp_ray) < 1 then the code will set
!          following map boundaries: min_i_nperp=0 and max_i_nperp=1. 
!
! n_contour_plot_disp is the number of contours for D(ReN_perp,Im_N_perp)
!          It should be =< n_contour_plut_disp_a
!-----------------------------------------------------------------------
! to plot cold plasma dispersion function D_cold(N_perp) at given points
! to plot.ps file            
!-------------------------------------------------------------------------
! i_plot_disp_cold  It used only in grill_lh to plot D in initial point
!                   =0 do not plot D_cold(N_perp)
!                   =1 plot D(N_perp)
!-----------------------------------------------------------------------
! n_plot_disp_cold=0 do not plot D(ReN_perp)
!          0< n_plot_disp_cold =< n_plot_disp_colda is the number of major radius
!          points where D_cold(Re0N_perp) will be ploted
!
! r_plot_disp_cold(n_plot_disp) major radius[m] where D_cold(N_perp)
!                               will be plotted
!
! s_poloid_plot_disp_cold(n_plot_disp) poloidal distance [m] where D_cold(N_perp)
!                                      will be plotted
!
! point_plot_disp_cold ='poloidl_dist' to create D(ReN_perp) plots at given
!                   s_poloid_plot_disp_cold() as deafault
!                 ='major_radius'  to create D(ReN_perp) plots at given
!                   r_plot_disp_cold()
!----------------------------------------------------------------------
! i_plot_wave_normal_cold! for plotting cold wavenormal to plot.ps file
!                        !  =1 to plot wavenormal at the initial point
!                        !     and output points, specified for cold plasma
!                        !     dispersion function D_cold plotting, at
!                        !     r_plot_disp_cold or s_poloid_plot_disp_cold
!                        !  =0 do not plot (as default)
!------------------------------------------------------------------------
!  For characteristic frequencies plotted along the straight line:
!  It works for i_plot_b.eq.1
!  Frequencies are electrons:  plasma, gyroharmonics, UH, f_R=0, f_L=0
!                  ions:       ion plasma, gyroharmonic, LH
!
!  r_freq,z_freq, cordinates of the line edge point [m]
!  alpha_freq  is the toroidal angle of the line [degree] 0 <alpha_freq<360
!              =0 r coordinate of the line is directed along 
!                 the major radius
!  beta_freq   is the angle between the line and the verticle axis Z
!              0 < beta_freq <180
!              =0 the line is directed along Z axis
!  dist_freq   is the line length   [m]
!
!  nsteps_freq  is the number of points used for plot.
!               It should be  nsteps_freq .le. 1000
!  n_ec_harmonics_freq  is the number of plotted 
!                       ec harmonics
!
!  max_plot_freq is the  maximal frequency at the plot [GHZ]
!-------------------------------------------------------------
!  npar_freq   N_parallel to plot X mode cutoof, 
!              It works for i_plot_b.eq.1
!
!  Plot with xdraw freqelec
!            xdraw freqion
!             Also, plot characteristic frequencies to *.bin files.

!--------------------------------------------------------------
c&output
      iwcntr=0
      iwopen=1
      iwj=2
      itools=0
      i_plot_b=1
      i_plot_d=0
      n_plot_disp=0
      point_plot_disp ='poloidl_dist'
      id_plot_disp=4
      if (n_plot_disp.gt.0) then
         do k=1,n_plot_disp
            s_poloid_plot_disp(k)=prmt6*k                    
         enddo
      endif

      rmax=2.d0
      rmin=1.d0
      if (n_plot_disp.gt.0) then
         do k=1,n_plot_disp
            r_plot_disp(k)=1.01*rmin+
     &                    (0.99*rmax-1.01*rmin)/n_plot_disp*k                    
         enddo
      endif

      point_plot_disp_cold ='poloidl_dist'
      i_plot_disp_cold=0
      n_plot_disp_cold=0
      if (n_plot_disp_cold.gt.0) then
         do k=1,n_plot_disp
            s_poloid_plot_disp_cold(k)=prmt6*k                    
         enddo
      endif
      if (n_plot_disp_cold.gt.0) then
         do k=1,n_plot_disp_cold
            r_plot_disp_cold(k)=1.01*rmin+
     &                    (0.99*rmax-1.01*rmin)/n_plot_disp_cold*k                    
         enddo
      endif
      i_plot_wave_normal_cold=0

      number_map_points_real_nperp=10
      number_map_points_image_nperp=10  
      ratio_min_r_nperp=0.5d0
      ratio_max_r_nperp=1.5d0
      ratio_min_i_nperp=0.d0
      ratio_max_i_nperp=2.5d0  
      n_contour_plot_disp=n_contour_plot_disp_a
      r_freq=1.49d0
      z_freq=0.d00 
      dist_freq=1.28d0 
      alpha_freq=180.d0
      beta_freq=90.d0
      nsteps_freq=780
      n_ec_harmonics_freq=6 
      max_plot_freq=200.d0 !GHZ 
      npar_freq=0.d0
c&end
!/plasma/
!-------------------------------------------------------------------------
!Plasma parameters
!-------------------------------------------------------------------------
! nbulk>=1 is a number of plasma components
!        It should be nbulk.le.nbulka
!        nbulka is a maximal number of plasma components
!        nbulka is a parameter which is set in param.i file
!----------------------------------------------------
! izeff =0 zeff will be calculated using the given ions;
!          electron density will be calculated using ions;
!       =1 zeff, electron density and ion densities with(i), i=2,nbulk-2
!          are given,
!          ion densities(i) i=nbulk and i= nbulk-1 will be calculated 
!          using Zeff, electron density and ion's densities(i), i=2,nbulk-2.
!          In this case it should be nbulk.ge.3
!       =2 zeff, electron and ion (if nbulk>1) densities are given,
!          and zeff is not recalculated from the plasma components;
!       =3 Use eqdsk pres (pressure). Let temperature T_E=T_i
!          pres=dens1(k,1)*temp1(k,1)+
!          Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
!          In this case we will calculate Zeff(rho),
!          dens_electron(rho) and T_e(rho)=T_i(rho)
!       =4 Use eqdsk pres (pressure), the given temperature
!          profiles T_i(rho) (i=1,nbulk) and the given Z_eff(rho).
!          nbulk should be .ge. 3
!          pres=dens1(k,1)*temp1(k,1)+
!          Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
!          In this case we will calculate dense(1)(rho),
!          dense(nbulk)(rho) and dense(nbulk-1)(rho)
! -----------------------------------------------------
! idens (0 - analytic, 1 - spline) representation of
! the density, temperature and zeff radial profiles
! -----------------------------------------------------
!   temp_scale(nbulka),den_scale(nbulka) are the parameters to multiply
!   the given temperature and density profiles
! -----------------------------------------------------
! ndens is the number of points for the input radial density and 
!       temperature profiles
!------------------------------------------------------
! nonuniform_profile_mesh= 'enabled' use nonuniform small radius mesh for input
!                                spline profiles (works for idens=1 only)
!                = 'disabled'    do not use nonuniform mesh (default)
!--------------------------------------------------------
c&plasma
      nbulk=nbulka
      izeff=2
      idens=1
      do i=1,nbulka          
         temp_scale(i)=1.d0
         den_scale(i)=1.d0
      enddo
      ndens=ndensa
      nonuniform_profile_mesh='disabled'
c&end

!/species/
! plasma component charges charge(i)=mod(charge(i)/charge_electron)
! -----------------------------------------------------
! charge(1) =1 electrons
! charge(i) i=1,nbulk   charge(i+1) must be ge.charge(i)
! charge(i) i=1,nbulk
! -----------------------------------------------------
! plasma components mass dmas(i)=Mass(i)/Mass_electron
! -----------------------------------------------------
! dmas(1) 1 electrons
! dmas(i)   i=1,nbulk
! -----------------------------------------------------
c&species 
      do i=1,nbulka          
         charge(i)=1.d0
         dmas(i)=3674.d0
      enddo
c      charge(1)=1.d0
c      charge(2)=1.d0
      dmas(1)=1.d0
c      dmas(2)=3674.d0
c&end

!/varden/
! the density variation
! -----------------------------------------------------
!   var0 is an amplitude of the density variation (del_n_0)
!   (see Manual, 3.37...)
!   denm is the number of the poloidal mode in the density variation(l_theta)
!   denn is the number of the toroidal mode in the density variation(l_phi)
!   an   is the radial localization of the variation (rho_0)
!   sigman is the parameter that characterizes the radial thickness
!          of the density fluctuation    
! -----------------------------------------------------
c&varden
      var0=0.d0
      denm=1.d0
      denn=15.d0
      an=0.5d0
      sigman=0.1d0 
c&end


!/denprof/
! -----------------------------------------------------
!Analytic radial profiles (idens=0).  Splines (idens=1).
!dense(i)=(dense0(i)-denseb(i))*(1-rho**rn1de(i))**rn2de(i)+denseb(i)
!-----------------------------------------------------
!        if(izeff.eq.0) then
!           zeff will be calculated using the given ions;
!	    nbulk1=nbulk
!	 if(zeff.eq.1 ) then          
!            zeff, the electron density and ion densities with(i) with i=2,nbulk-2
!            are given,
!            Ions components with i=nbulk and i= nbulk-1 will be
!            calculated will be calculated using Zeff, the electron density and 
!            ion's densities(i) with i=2,nbulk-2.
!            In this case it should be nbulk.ge.3
!
!            if (nbulk.eq.1) nbulk1=1
!            if (nbulk.eq.2) then
!	         nbulk1=2
!	     endif
!            if (nbulk.gt.2) nbulk1=nbulk-2
!	 endif
! -----------------------------------------------------
! dense0(i)   central density in 10**19 m**(-3) i=1,nbulk1
! -----------------------------------------------------
! denseb(i)  edge density in 10**19 m**(-3) i=1,nbulk1
! -----------------------------------------------------
! rn1de(i) i=1,nbulk1
! -----------------------------------------------------
! rn2de(i) i=1,nbulk1
! -----------------------------------------------------
c&denprof 
      do i=1,nbulka 
         dense0(1)=6.200d+0
         dense0(2)=6.200d+0
         rn1de(i)=2.d+0
         rn2de(i)=1.00d+0
      enddo
!      dense0(1)=6.200d+0
!      dense0(2)=6.200d+0
!      denseb(1)=0.2d-0
!      denseb(2)=0.2d-0
!      rn1de(1)=2.d+0
!      rn1de(2)=2.d+0
!      rn2de(1)=1.00d+0
!      rn2de(2)=1.00d+0
c&end
!/tpoprof/
! Ratio tpop=T_perp/T_parallel
! tpop(i)=(tp0(i)-tpb(i))*(1-rho**rn1tp(i))**rn2tp(i)+tpb(i)
! -----------------------------------------------------
! tp0(i) =           central T_perp/T_parallel i=1,nbulk
! -----------------------------------------------------
! tpb(i) =           boundary T_perp/T_parallel i=1,nbulk
! -----------------------------------------------------
! rn1tp(i) i=1,nbulk
! -----------------------------------------------------
! rn2tp(i)  i=1,nbulk
! -----------------------------------------------------

c&tpopprof
      do i=1,nbulka
         tp0(i)=1.0d0
         tpb(i)=1.0d0 
         rn1tp(i)=2.0d0
         rn2tp(i)=1.0d0
      enddo
!     tp0(1)=1.0d0
!     tp0(2)=1.d0
!     tpb(1)=1.0d0
!     tpb(2)=1.0d0
!     rn1tp(1)=2.0d0
!     rn1tp(2)=2.0d0
!     rn2tp(1)=1.0d0
!     rn2tp(2)=1.0d0
c&end

!/vflprof/
! drift velocity parallel B (m/sec) 
! vflow(i)=(vfl0(i)-vflb(i))*(1-rho**rn1vfl(i))**rn2vfl(i)+vflb(i)
! -----------------------------------------------------
! vfl0(i)     central vflow in m/sec  i=1,nbulk
! -----------------------------------------------------
! vflb(i)     boundary vflow in m/sec i=1,nbulk
! -----------------------------------------------------
! rn1vfl(i) i=1,nbulk
! -----------------------------------------------------
! rn2vf(i)  i=1,nbulk
! -----------------------------------------------------

c&vflprof 
      do i=1,nbulka 
         vfl0(i)=0.0d+0
         vflb(i)=0.0d+0
         rn1vfl(i)=2.d0
         rn2vfl(i)=1.0d0
      enddo
!      vfl0(1)=0.0d+0
!      vfl0(2)=0.0d0
!      vflb(1)=0.0d+0
!      vflb(2)=0.0d0
!      rn1vfl(1)=2.d0
!      rn1vfl(2)=2.0d0
!      rn2vfl(1)=1.0d0
!      rn2vfl(2)=2.0d0
c&end


!/zprof/
! -----------------------------------------------------
! zeff=(zeff0-zeffb)*(1-rho**rn1zeff)**rn2zeff+zeffb
! -----------------------------------------------------
! zeff0   central Z_eff
! zeffb   boundary Z_eff
! rn1zeff zeff=(zeff0-zeffb)*
! rn2zeff      (1-rho**rn1zeff)**rn2zeff+zeffb
!-----------------------------------------------------
c&zprof
      zeff0=2.0d0
      zeffb=2.0d0
      rn1zeff=2.d0
      rn2zeff=1.d0
c&end

!/tprof/
! Average temperature tempe=(T_parallel+2*T_perp)/3
! tempe(i)=(te0(i)-teb(i))*(1-rho**rn1te(i))**rn2te(i)+teb(i)
! -----------------------------------------------------
! te0(i) =at0(i)    central temperature in kev	i=1,nbulk
! -----------------------------------------------------
! teb(i) =ateb(i)    boundary temperature in kev i=1,nbulk
! -----------------------------------------------------
! rn1te(i) i=1,nbulk
! -----------------------------------------------------
! rn2te(i)  i=1,nbulk
! -----------------------------------------------------

c&tprof
      do i=1,nbulka 
         ate0(i)=3.0d0
         ateb(i)=5.0d-2
         rn1te(i)=2.d0
         rn2te(i)=1.0d0
      enddo

!      ate0(1)=3.0d0
!      ate0(2)=3.0d0
!      ateb(1)=5.0d-2
!      ateb(2)=5.0d-2
!      rn1te(1)=2.d0
!      rn1te(2)=2.d0
!      rn2te(1)=1.0d0
!      rn2te(2)=1.0d0
c&end



!/grill/
!------------------d3d---LH------one ray--------------------------------
!------------------LH/EBW-Starting-inside-plasma-----------------------
!  Grill conditions  for istart=2 (start point inside the plasma)
!----------------------------------------------------------------------
! i_n_poloidal =1         The input parameter is N_parallel(from grill).
!  (by default =1)        N_phi,N_theta are calculated from given N_parallel 
!                         N_rho=N_perpendicular(N_parallel) is determined 
!                         from the dispersion relation. It is directed
!                         along +,- gradient(psi) 
!
! i_n_poloidal =2         The input parameters: N_parallel(from grill)
!                         and  n_theta_pol. By default N_theta=0. 
!                         N_perpendicular(N_parallel) is determined 
!                         from the dispersion relation. 
!                         N_phi is calculated from N_parallel and N_theta
!                         N_rho is calculated form N_perpendicular, N_parallel
!                         and N_theta. 
!                         It is directed along +,- gradient(psi)
!
! i_n_poloidal=3          The given variables: N_parallel and the angle
!                         0<ksi_nperp<180 between the vector N_perpendicular 
!                         and gradient(psi). By default ksi_nperp=0.
!                         N_perpendicular(N_parallel) is determined 
!                         from the dispersion relation.
!                         N_phi,N_theta and N_rho are calculated from
!                         N_parallel,N_perpendicular and ksi_nperp.
!
! i_n_poloidal=4          The given variables:N_toroidal and
!                         N_poloidal. It case uses i_vgr_ini set in /waves/
!                         to choose the direction of the small radial N_rho
!                         component. To launch the ray inside the plasma
!                         i_vgr_ini=1 or to the plasma edge i_vgr_ini=-1 
!---------------------------------------------------------------------
! n_theta_pol            The poloidal refractive index component
!                         It is used for i_n_poloidal =2     
!                         By_default n_theta_pol=0.
!----------------------------------------------------------------------
! ksi_nperp               (degrees) the angle 0<ksi_nperp<180
!                         between the vector N_perpendicular 
!                         and gradient(psi). By default ksi_nperp=0.
!---------------------------------------------------------------------
!  Calculation of the small radius value near the plasma edge
!  where LH or FW have cutoff:  
!  i_rho_cutoff=0 (default) no calculations
!              =1 use these calculations
!------------------------------------------------------------------
!  rho_step_find_LHFW_cutoff  is the non dimensional small radius step
!                            used in subroutine  rho_ini_LHFW
!                            It is used at i_rho_cutoff=1
!-------------------------------------------------------------------
!  rho_initial_find_LHFW_cutoff  is the initial small radius point. 
!                            As default rho=1- rho_step_find_LHFW_cutoff
!                            It is used at i_rho_cutoff=1
!                            (default=1.-1.d-3)
!--------------------------------------------------------------------
!  ngrill  is a number of the poloidal grill angles
!          It is required that ngrill.le.ngrilla, 
!          where ngrilla is parameter in param.i
!----------------------------------------------------------------------
!  igrillpw options specifying N_parallel power spectra
!           =1 power=powers/nnkpar,
!           =2 power=sin**2x/x**2,
!           =3 power=exp-((npar-anmin)/anmax)**2    [default=1]
!----------------------------------------------------------------------
!  igrilltw specifies the form poloidal variation of power,
!                    =1 uniform over height, =2 cos**2 variation.
!----------------------------------------------------------------------
!  rhopsi0(1:ngrill) initial small radius for wave front
!                    (0<rhopsi0<1)
!  rhopsi0(i)=...    i=1,ngrill
!----------------------------------------------------------------------
!  thgrill(1:ngrill) poloidal  angle of grill, measured counter
!                    clockwise from horizontal through the
!                    magnetic axis (degrees).
!  thgrill(i)=...    i=1,ngrill (degree)         [default=0.d0]
!---------------------------------------------------------------------
!  phigrill(1;ngrill) is a toroidal grill angle of grill
!                              (degrees)
!  phigrill(i)=... i=1,ngrill (degree)         [default=0.d0]
!----------------------------------------------------------------------
! height(1:ngrill) is a poloidal length (m) of grill
!                 (giving poloidal power distribution of each grill).
! height(i)=...   i=1,ngrill                  [default=0.2d0]
!----------------------------------------------------------------------
! nthin(1:ngrill) is a number of rays near the each poloidal
!                 center, simulating a grill
! nthin(i)=...    i=1,ngrill       [default: nthin(1)=1]
!----------------------------------------------------------------------
!  anmin(1:ngrill)  position of the left bound
!                   of power spectrum P(n_parallel) (Can be neg).
!  anmin(i)=...     i=1,ngrill
!----------------------------------------------------------------------
!  anmax(1:ngrill)  position of the right bounds
!                   of power spectrum P(n_parallel) (Can be neg).
!  anmax(1)=...     i=1,ngrill
!---------------------------------------------------------------------
!  nnkpar(1:ngrill)  number of points  of power spectrum
!                    P(n_parallel)
!  nnkpar(i)=...     i=1,ngrill
!----------------------------------------------------------------------
!  powers(1:ngrill)  power in one grill (MWatts)
!  (total power of grill(in MWatts) will be powtott=sum{powers}
!  powers(i)=...     i=1,ngrill
!-------------------------------------------------------------------
!  below are for i_n_poloidal=4 case, set (N_toroidal, N_poloidal)
!----------------------------------------------------------------------
!  antormin(1:ngrill)  position of the left bound
!                   of power spectrum P(n_toroidal) (Can be neg).
!  antormin(i)=...     i=1,ngrill
!----------------------------------------------------------------------
!  antormax(1:ngrill)  position of the right bounds
!                   of power spectrum P(n_toroidal) (Can be neg).
!  antormax(1)=...     i=1,ngrill
!---------------------------------------------------------------------
!  nnktor(1:ngrill)  number of points  of power spectrum
!                    P(n_toroidal)
!  nnktor(i)=...     i=1,ngrill
!----------------------------------------------------------------------
!  anpolmin(1:ngrill)  position of the left bound
!                   of power spectrum P(n_poloidal) (Can be neg).
!  anpolmin(i)=...     i=1,ngrill
!----------------------------------------------------------------------
!  anpolmax(1:ngrill)  position of the right bounds
!                   of power spectrum P(n_poloidal) (Can be neg).
!  anpolmax(1)=...     i=1,ngrill
!---------------------------------------------------------------------
!  nnkpol(1:ngrill)  number of points  of power spectrum
!                    P(n_poloidal)
!  nnkpol(i)=...     i=1,ngrill
!---------------------------------------------------------------------
!  ilaunch=0, no effect (the default)
!
!  ilaunch=1, to launch a single ray at r0launch,phi0launch,z0launch
!             in the plasma (meters and degs)
!
!  This option is added for comparison with other codes.
!  r0launch is the major radius of the launch point [m]
!  z0launch is the vertical position of the launch point [m]
!  phi0launch is the toroidal angle of the launch point [degree] 
!
!  ilaunch=2, to launch rays at rlaunch(j,i),phigrill(i),zlaunch(j,i)
!             in the plasma (meters and degs)
!           
!             i=1,ngrill
!               j=1,nthin(i)
!
!             total power at each grill (i) is  powers(i)
!             with weight_power_launch(j,i) non dimensional             
!             sum{j=1,ntin(j)}(weight_power_launch(j,i))=1 Normalization
!             power at each ray witn index (j,i) is 
!                   powers(i)*weight_power_launch(j,i)!
!---------------------------------------------------------------------
!       i_grill_pol_mesh: option specifying the poloidal mesh wtheta(j)
!                         near the central grill angle thgrill(i)
!                         =1 equispaced mesh
!                            wtheta(j)-wtheta(j-1)=zdth=Const(default)
!                         =2 poloidal mesh will be chosen to get the equal
!			     power fpwth(j) for all rays near the central
!                            grill angle fpwth(j)=1/nthini
!
!------------------------------------------------------------------------
!       i_grill_npar_ntor_npol_mesh: option specifying the refractive
!                         index meshes.
!
!                         For  i_n_poloidal=1,2,3 it specifies the
!                         n_parallel mesh anzin(n) for the power
!                         spectrum pwcpl(n) n=1,...,nnkpari
!                         =1 equispaced mesh
!                            anzin(n)-anzin(n-1)=hnpar=Const (default)
!                         =2 n_parallel mesh will be chosen to get equal
!			     power pwcpl(n) for all rays in the given power
!                            spectrum  pwcpl(n)=1.d0/nnkpari
!                            pwcpl(n)=power_spectrum(anzin(n))*
!                                     delta_npar_bin(n)= 1.d0/nnkpari
!
!                            For  i_n_poloidal=4 it specifies two meshes:
!                            a) n_toroidal mesh anztorin(ntor) and
!                            b) n_poloidal mesh anzpolin(npolmesh)
!                            for the power spectrum
!                            pwcpl_tp(1:nnktori,1:nnkpoli)=pwcpl_t*pwcpl_t
!                         =1 equispaced meshs (default)
!                            anztorin(ntor)- anztorin(ntor-1)=hntor=Const
!                            anzpolin(npol)- anzpolin(npol-1)=hnpol=Const
!                         =2 the meshes anztorin(1:nntori) anzpolin(1:nnkpoli)
!                            will be chosen to get the equal
!			     power pwcpl_tp(ntor,npol) for all rays in
!                            the given power spectrum
!                            pwcpl_tp(ntor,npol)=1.d0/(nnktori*nnkpoli)
!
!------------------------------------------------------------------------
c&grill
      i_n_poloidal=1
      n_theta_pol=0.d0
      ksi_nperp=0.0d+0
      i_rho_cutoff=0
      rho_step_find_LHFW_cutoff=1.0d-3
      rho_initial_find_LHFW_cutoff=1.d0-rho_step_find_LHFW_cutoff
      ngrill=1
      igrillpw=1
      igrilltw=2
      rhopsi0(1)=0.97d+00
      thgrill(1)= 0.0d+0
      phigrill(1)=0.0d+0
      height(1)=0.20d+0
      nthin(1)=1
      anmin(1)=4.000d+0
      anmax(1)=6.000d+0
      nnkpar(1)=1
      powers(1)=1.0d+0     !MWATT
      antormin(1)=0.1d0
      antormax(1)=1.d0
      nnktor(1)=1 
      anpolmin(1)=0.1d0
      anpolmax(1)=1.d0 
      nnkpol(1)=1    
      ilaunch=0 
      r0launch=1.6d+0
      z0launch=0.d+0
      phi0launch=0.d+0
      i_grill_pol_mesh=1
      i_grill_npar_ntor_npol_mesh=1
c&end
!     namelist /rz_launch_grill_tab/  It is used at ilaunch=2 only
!-----------------------------------------------------------
!     direct grill launch points from namelist /rz_launch_grill_tab/
!     It is used at ilaunch=2 case only
!
!     rlaunch(j,i),zlaunch(j,i) [meters]
!     and weight_power_launch(j,i),[non dimentional]
!     Here
!     i=1,ngrill
!       j=1,nthin(i)     
!     
!     weight_power_launch(j,i) should be normalized using
!     following formula:             
!     sum{j=1,ntin(j)}(weight_power_launch(j,i))=1  
!
!     power at each ray witn index (j,i) is powers(i)*weight_power_launch(j,i)
!---------------------------------------------------------
!    Bellow are described used input data for rz launch
!    which were set in namelist /grill/.
!
!     ngrill is a number of grills was set in namelist /grill/  
!
!     phigrill(i=1:ngrill) is a toroidal angle of grill, measured counter
!                          clockwise from horizontal through the
!                          magnetic axis (degrees).
!                          They were set in namelist /grill/.
!  
!     nthin(i=1:ngrill)    numbers of poloidal points at each grill
!                          They were set in namelist /grill/.
!  
!     total powers at each grill (i) is  powers(i=1,ngrill) .
!                          It were set in namelist /grill/.
!------------------------------------------------------------------  
!  The example for ngrill=2, ntin(1)=3, ntin(2)=4
!--------------------------------------------------------------
!&rz_launch_grill_tab
!for grill with i=1
!
!      prof_2d(1,1)=
!      rlaunch(j=1,i=1), zlaunch(j=1,i=1), weight_power_launch(j=1,i=1),
!      rlaunch(j=2,i=1), zlaunch(j=2,i=1), weight_power_launch(j=2,i=1),
!      rlaunch(j=3,i=1), zlaunch(j=3,i=1), weight_power_launch(j=3,i=1)
!
!for grill with i=2
!
!      prof_2d(1,2)=
!      rlaunch(j=1,i=2), zlaunch(j=1,i=2), weight_power_launch(j=1,i=2),
!      rlaunch(j=2,i=2), zlaunch(j=2,i=2), weight_power_launch(j=2,i=2),
!      rlaunch(j=3,i=2), zlaunch(j=3,i=2), weight_power_launch(j=3,i=2),
!      rlaunch(j=4,i=2), zlaunch(j=4,i=2), weight_power_launch(j=4,i=2),
!&end
!-------------------------------------------------------------------------
!&rz_launch_grill_tab
!for grill with i=1
!     
! prof_2d(1,1)=
!     1.6d0, 0.d0,  0.3333d0,
!     1.7d0, 0.d0,  0.3333d0,
!     1.8d0, 0.d0,  0.3333d0,
!for grill with i=2
! prof_2d(1,2)=
!     1.60d0,  0.d0, 0.25d0,
!     1.65d0,  0.d0, 0.25d0,
!     1.70d0,  0.d0, 0.25d0,
!     1.75d0,  0.d0, 0.25d0,     
!&end
!---------------------------------------------------------------------------
      if (ilaunch.eq.2) then
         do i=1,ngrill
            do j=1,nthin(i)
               rlaunch(j,i)=1.6d0
               zlaunch(j,i)=0.0d0
               weight_power_launch(j,i)=1.d0/dfloat(nthin(i))
            enddo
         enddo
      endif 
!----------------------------------------------------------------------
!/eccone/
!-----------------------CANONICAL 2004 ITER TEST data-------------
!     Use equilib.dat= g521022.01000, or equivalent.
!
!     The namelist section specifies ECR cones  for istart=1 
!           (ray cones start outside the plasma).
!
!     Multiple source locations and launch conditions are implemented.
!     ncone=number of source cones. [default=1] [Max is parameter nconea]. 
!
!     For multiple sources (ncone.gt.1), is is necessary to set
!     ncone values for each of the the namelist variables given below:
!     powtot,zst,rst,phist,alfast,betast,alpha1,alpha2(only for raypatt
!     ="genray").  
!     The specifications of number of rays per cone do not not vary 
!     from cone to cone (i.e., na1,na2,gzone,mray(*), cr(*)) do not 
!     vary with cone number.
!    
!     powtot= total power from antenna(MW)
!
!     Two systems for specification of the ECR cone are provided,
!       chosen by raypatt:
!
!     raypatt='genray',  specify ray pattern per following
!                        genray method:	
!     zst (m)   initial z of the cone vertex
!     rst(m)    initial r of the cone vertex
!     phist(degree) initial toroidal angle phi of cone vertex,
!                   measure from x-z plane.
!     alfast(degree) toroidal angle measure from R-vector through
!                    source
!     betast(degree) poloidal angle measured from z=constant plane,
!                      positive above plane, negative below.
!     alpha1(degree) angle cone width
!     alpha2(degree) starting angle along cone
!     na1 number of cones (0 for central ray only)
!     na2 number of rays at cone(for na1.ge.0)
!
!     raypatt='toray',  specify ray pattern per the following
!                       toray method:  [Defn of betast changed,
!                       and there are additional namelist inputs.]
!     zst (m)   initial z of the cone vertex
!     rst(m)    initial r of the cone vertex
!     phist(degree) initial toroidal angle phi of cone vertex,
!                   measure from x-z plane.
!     alfast(degree) toroidal angle measure about z-axis
!                    from R-vector through the source.
!     betast(degree) poloidal polar angle measured from positive
!                    z-axis. [DIFFERENT FROM RAYPATT='GENRAY'!]
!     alpha1(degree) angle cone width, half-power width of the beam.
!     gzone:    if 0 then 48 ray case, as specified by mray() below
!               if 1 then there can only be 1 ray, the central ray.
!               if .gt.1 then describes number of elements in mray
!     mray(*):  if gzone .gt.0, use gaussian formulation with this number
!        of rays in corresponding annular zone, otherwise use the
!        usual 1,5,12,12,18  (48 ray) arrangement.
!        mray(1) is effectively 1.
!     cr(*):    azimuthal phase of ray pattern for each zone, in radians;
!       same size array as mray for gaussian formulation.
!      Standard setting for gzone=0 is 0.0,0.1,0.05,-0.05,0.05.
!      
! ----input data for disk to disk launching, raypatt='diskdisk':
!     power distribution at the first disk has the gaussian variation
!     w.r.t. disk radius rho:
!     power(rho)=Const*exp(-rho/sigma_launching_disk)**2
!
!     The power destribution at the launching disk will be determine by the
!     parameter : 0. < part_gauss_power =< 1.
!     part_gauss_power=Integral(0,rho_launchin_disk){rho*d(rho)*
!     2/(sigma_launching_disk)**2*exp(-rho/sigma_launching_disk)**2}
!
!     input parameters for diskdisk case are
!     sigma_launching_disk [m] It works at 1.<part_gauss_power
!     d_disk is distance between the disks perpendicular to disks [m]
!     part_gauss_power  It is from 0. to 1.
!                              if 0.<part_gauss_power<1.    
!                              sigma_launching_disk will be calculated using:
!                                sigma_launching_disk=rho_launching_disk/
!                                dsqrt(dlog(1.d0/(1-part_gauss_power))))
!                             
!                              If part_gauss_power.ge.1 then
!                              sigma_launching_disk will be taken from         
!                              genray.dat and the code will recalculate
!                              part_gauss_power using given
!                              sigma_launching_disk
!     rho_focus_disk  [m] the second disk radius
!     n_mesh_disk_radial_bin is the number of radial bins at the first disk
!     n_mesh_disk_angle_bin(n_mesh_disk_radial_bin) are the number 
!                 of angle bins at each radius bin
!     initial_azimuth_angle_degree((n_mesh_disk_radial_bin) are initial
!                  angles on the first disk around the central ray
!
!     The central ray will be directed from the center of the first
!     disk to the center of the second disk
!     The other rays will be directed from the first disk
!     to edge of the second disk
!
! ----the input data for diskbeam launching, raypatt='diskbeam'
!     Rays will be launched from the launching disk parallel to
!     the central ray
!     power distribution at the launhing disk has the gaussian form
!     on disk radius: rho
!     power(rho)=Const*exp(-rho/sigma_launching_disk)**2
!
!     The power at the launching disk will be determine by the
!     parameter : 0. < part_gauss_power =< 1.
!     part_gauss_power=Integral(0,rho_launchin_disk){rho*d(rho)*
!     (2/sigma_launching_disk)**2*exp(-rho/sigma_launching_disk)**2}
!
!                               
!     input parameters for diskdisk case are
!     sigma_launching_disk [m] It works for 1.<part_gauss_power
!     part_gauss_power  It is from 0. to 1.
!                              if <part_gauss_power<1    
!                              sigma_launching_disk will be calculated using:
!                                sigma_launching_disk=rho_launching_disk/
!                                dsqrt(dlog(1.d0/(1-part_gauss_power)))
!                             
!                              If part_gauss_power.ge.1 then
!                              sigma_launching_disk will be taked from 
!                              genray.dat.  Then the code will recalculate
!                              part_gauss_power using given
!                              sigma_launching_disk
!     rho_launching_disk  [m] radius of the launching disk
!     rho_focus_disk      [m] radius of the second disk   
!     n_mesh_disk_radial_bin is the number of radial bins at the first disk
!     n_mesh_disk_angle_bin(n_mesh_disk_radial_bin) are the number 
!                 of angle bins at each radius bin
!                 It should be n_mesh_disk_angle_bin(1)=1.
!
!     initial_azimuth_angle_degree(n_mesh_disk_radial_bin) are initial
!                  angles on the first disk around the central ray 
!                  directed clockwise from the vector R_0 
!                
!     The central ray will be directed from the center of the launching
!     disk. The central ray direction is set by angles:
!     alfast(degree) and betast(degree) 
!-----------------------------------------------------------------
c &eccone
      ncone=1
      powtot(1)=1.0d0          !MWatt

      raypatt='genray'
      zst(1)=+4.11d+0
      rst(1)=6.4848+00
      phist(1)=+0.000d+0
      betast(1)=-56.075d0       !Equals -(polar_angle-90.)
      alfast(1)=+137.84d0
      alpha1(1)=1.177d+00
      alpha2(1)=+1.500d+1
      na1=3
      na2=10
      na2=10

! raypatt='toray'
! zst(1)=+4.11d+0
! rst(1)=6.4848+00
! phist(1)=+0.000d+0
! betast(1)=146.075d0          !Measured from z-axis
! alfast(1)=+137.84d0
! alpha1(1)=1.177d+00
! gzone=5
! mray(1)=1
! mray(2)=5
! mray(3)=12
! mray(4)=12
! mray(5)=18
! mray(6)=24
! mray(7)=24
! cr(1)=0.0d0
! cr(2)=0.1d0
! cr(3)=0.05d0
! cr(4)=-0.05d0
! cr(5)=0.05d0
! cr(6)=-0.025d0
! cr(7)=0.025d0

! nray_in=48 Need not be specified for gzone.ne.0
! cr=defaults.

!     raypatt='diskdisk'
! zst(1)=+4.11d+0
! rst(1)=6.4848+00
! phist(1)=+0.000d+0
! betast(1)=-56.075d0  !Equals -(polar_angle-90.)
! alfast(1)=+137.84d0
! alpha1(1)=1.177d+00
      d_disk=0.50d0 !m
      sigma_launching_disk=0.025d0 !m
      part_gauss_power=1.1d0 !in this case sigma_launching_disk
                             !will be taken from input genray.dat
                             !then the code will recalculate
                             !part_gauss_power using given
                             !sigma_launching_disk
           
      rho_launching_disk=0.1d0 !m
      rho_focus_disk=0.015d0   !m
      n_mesh_disk_radial_bin=1
      n_mesh_disk_angle_bin(1)=1
      n_mesh_disk_angle_bin(2)=1
      n_mesh_disk_angle_bin(3)=1
      n_mesh_disk_angle_bin(4)=1
      n_mesh_disk_angle_bin(5)=1
      initial_azimuth_angle_degree(1)=0.d0 !degree
      initial_azimuth_angle_degree(2)=0.d0 !degree
      initial_azimuth_angle_degree(3)=0.d0 !degree
      
      initial_azimuth_angle_degree(4)=0.d0 !degree
      initial_azimuth_angle_degree(5)=0.d0 !degree

!     raypatt='diskbeam'
! zst(1)=+4.11d+0
! rst(1)=6.4848+00
! phist(1)=+0.000d+0
! betast(1)=-56.075d0  !Equals -(polar_angle-90.)
! alfast(1)=+137.84d0
! alpha1(1)=1.177d+00
!      sigma_launching_disk=0.025d0 ![m] wokrs at 0<part_gauss_power<1
!      part_gauss_power=1.1d0 !in this case sigma_launching_disk
                             !will be taken from input genray.dat
                             !then the code will recalculate
                             !part_gauss_power using given
                             !sigma_launching_disk   
!      rho_launching_disk=0.1d0 !m
!      
!      n_mesh_disk_radial_bin=1
!      n_mesh_disk_angle_bin(1)=1
!      n_mesh_disk_angle_bin(2)=1
!      n_mesh_disk_angle_bin(3)=1
!      n_mesh_disk_angle_bin(4)=1
!      n_mesh_disk_angle_bin(5)=1
!      initial_azimuth_angle_degree(1)=0.d0 !degree
!      initial_azimuth_angle_degree(2)=0.d0 !degree
!      initial_azimuth_angle_degree(3)=0.d0 !degree
!      initial_azimuth_angle_degree(4)=0.d0 !degree
!      initial_azimuth_angle_degree(5)=0.d0 !degree

c &end
!/dentab/at uniform grid rho(i)=(i-1)/(ndens-1) i=1,...,ndens
!         for nonuniform_profile_mesh='disabled'
!--------------------------------------------------------------------------
! density profiles (table data, case: idens=1)	dens1(ndens,nbulk)
!--------------------------------------------------------------------------
! ndensa (a parameter in param.i) is max number of points in the 
!   small radius direction.
! nbulka (a parameter in param.i) is a max number of plasma components.
! Input of profiles is set up so spline profiles can be input in tables
!   of size specified through namelist variables ndens and nbulk.
! ndens (variable)    is  number of points in small radius direction
!                     (set in namelist /plasma/). Must be .le. ndensa.   
! nbulk (variable)    is number of plasma components, must have: 
!                     nbulk.le.nbulka, and
!                     first component is for electrons
! nbulk1 is number of density components which must be specified.
! nbulk1 is calculated in dinit_mr subroutine,
! The fragment of dinit_mr is given here to understand the
! nbulk1 value
!
! The number of columns in dentab should be equal to nbulk.
! If nbulk1 < nbulk then we should put the density profiles
!   for the first nbulk1 plasma components in the table.
! The profiles for last (nbulk-nbulk1) plasma components can be arbitrary.
!
!--------------------------------------------------------------------------
!c     calculation of nbulk1
!      if(((izeff.eq.0).or.(izeff.eq.2)).or.(izeff.eq.3)) then
!c        izeff=0, zeff will be calculated using the given ions;
!c                 electron density will be calculated using ion's densities;
!c        izeff=2, zeff will not coincide with the plasma components
!c             =3  it uses eqdsk pres (pressure) and ion densities_i 
!c                 for i=2,... nbulk
!c                 Let temperature T_E=T_i
!c                 pres=dens1(k,1)*temp1(k,1)+
!c                      Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
!c                 In this case we will calculate Zeff(rho),
!c                 dens_electron(rho) and T_e(rho)=T_i(rho)
!c                  This case works for nbulk >1 only.
!c             =4  it uses eqdsk pres (pressure), zeff,ions densities
!c                 for i=2,... nbulk-2 (nbulk>2) and T_e(rho),T_i(rho)
!c                 pres=dens1(k,1)*temp1(k,1)+
!c                      Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
!c                 In this case we will calculate dens_electron(rho) and
!c                 ion densities for i=nbulk and i=nbulk-1)
!         nbulk1=nbulk
!      else
!c        (izeff=1 or izeff=4) ion densities(i) with i= nbulk and i=(nbulk-1) will
!c                 be calculated  using given
!c                 Zeff, the electron density and ion's densities(i), i=2,nbulk-2;
!         if (nbulk.le.2) nbulk1=nbulk
!         if (nbulk.eq.2) then
!	    write(*,*)'nbulk=2, Zeff must be equal charge(2)'
!           write(*,*)'Please check it or use the option izeff=0'
!	    stop
!	 endif
!         if (nbulk.gt.2) nbulk1=nbulk-2
!      endif !izeff
!
!      The case nbulk=1 is used often for ECR and EBW cases.
!      In these cases only the electron component is essential.
!
!      For nbulk=1 and izeff=2 case only the electron density
!      is used in dispersion relation.In this case  Z_effective
!      is used for current drive efficiency calculations.   
!------------------------------------------------------------------------
! dens1(ndens,nbulk) (10!!3/cm!!3)
!------------------------------------------------------------------------
! If  ((izeff.eq.0).or.(izeff.eq.3)) then the electron density
! will be calculated from the charge neutrality.
! In that case we can set the arbitrary values for the electron density
! dens1(k,1), k=1:ndens and should set nbulk1-1 ion densities:
! dens1(k,i), k=1:ndens, i=2:nbulk1.
! A constant radial step is assumed, 
! The first line (k=1, i=1:nbulk1) dens1(1,i) is for rho=0
! The last line (k=ndens, i=1:nbulk1) dens1(ndens,i) is for rho=1
! The example for  izeff=0, ndens=5, nbulk=3, nbulka=4
!  nbulk is a number of plasma species
!  nbulka is a maximal number of plasma species.
!  nbulka is set in param.i file.
!  It should be nbulka.ge.nbulk
!       
! column:    1       2      nbulk
!         electron   ion    ion
!         
! prof=      0.,     1.2,   1.3,       
!            0.,     2.2,   2.3,       
!            0.,     3.2,   3.3,       
!            0.,     4.2,   4.3,       
!            0.,     5.2,   5.3, 
!
!For izeff=1 case we should set the profiles of the electron density and 
!ion densities(i) with i=1,nbulk-2      
!The columns of ion densities(i) with i=nbulk-1 and i=nbulk should be
!fill in by arbitrary numbers (they can be zeros). 
!
!The example for izeff=1, ndens=5, nbulk=4, nbulka=5
!  nbulk is a number of plasma species
!  nbulka is a maximal number of plasma species.
!  nbulka is set in param.i file.
!  It should be nbulka.ge.nbulk
!
!
!colomn:      1       nbulk-2  nbulk-1  nbulk
!         electron    ion      ion      ion   
!         
! prof=      1.1,     1.0,     0.0,     0.0,       
!            0.9,     0.85,    0.0,     0.0,      
!            0.6,     0.55,    0.0,     0.0,     
!            0.4,     0.32,    0.0,     0.0,     
!            0.2,     0.15,    0.0,     0.0,
!
! ------------------------------------------------------------------------
! Here array prof(ncomp,ndens) was used for convenience in namelist.
! dens1(k,i)=prof(i,k) k=1:ndens, i=1:ncomp
! ------------------------------------------------------------------------
! If (izeff.ne.0) and (izeff.ne.3) then we should set the electron density
! dens1(k,1)  and ion densities dens1(k,i) i=2:nbulk1
!------------------------------------------------------------------------
c&dentab
c prof=   6.200000000E+00,  6.200000000E+00,  0.,  
c         6.140000000E+00,  6.140000000E+00,  0.,  
c         5.960000000E+00,  5.960000000E+00,  0., 
c         5.660000000E+00,  5.660000000E+00,  0.,  
c         5.240000000E+00,  5.240000000E+00,  0., 
c         4.700000000E+00,  4.700000000E+00,  0., 
c         4.040000000E+00,  4.040000000E+00,  0., 
c         3.260000000E+00,  3.260000000E+00,  0., 
c         2.360000000E+00,  2.360000000E+00,  0.,  
c         1.340000000E+00,  1.340000000E+00,  0.,  
c         2.000000000E-01,  2.000000000E-01,  0.,      
c&end
      do i=1,nbulka
         do k=1,ndens
            rho=1.d0*(k-1)/dfloat(ndens-1)
            dens1(k,i)=(dense0(i)-denseb(i))*(1-rho**rn1de(i))**rn2de(i)
     .                 +denseb(i)
         enddo
      enddo

c---------------------------------------------------------------------
c     Here dens1 will be in [10**19/m**3] for the genray.dat input file
c---------------------------------------------------------------      

!/temtab/
!--------------------------------------------------------------------------
! temperature profiles (table data, case: idens=1)	temp1(ndens,nbulk)
! Average temperature temp1=(T_parallel+2*T_perp)/3
!--------------------------------------------------------------------------
! It this namelist we must set electron temp1(ndens,1) and all ion
! species temp1(ndens,i) temperature (keV) {i=2:nbulk}
! Constant radial step is assumed.
! The first line (k=1, i=1:nbulk) temp1(1,i) is for rho=0
! The last line (k=ndens, i=1:nbulk) temp1(ndens,i) is for rho=1
! The example for   ndens=5, nbulk=3, ncomp=4
!                 
!         electron-1 ion-2  ion-nbulk 
! prof=      0.,     1.2,   1.3,       
!            0.,     2.2,   2.3,       
!            0.,     3.2,   3.3,       
!            0.,     4.2,   4.3,       
!            0.,     5.2,   5.3,
!
! In all cases temtab should has nbulk columns with the temperature
! profiles for all nbulk plasma components.    
! ------------------------------------------------------------------------
! Here array prof(ncomp,ndens) was used for convenience in namelist.
! temp1(k,i)=prof(i,k) k=1:ndens, i=1:nbulk
! ------------------------------------------------------------------------
c&temtab
c prof=    3.000000000E+00,  3.000000000E+00, 0.,  
c          2.970500000E+00,  2.970500000E+00, 0., 
c          2.882000000E+00,  2.882000000E+00, 0., 
c          2.734500000E+00,  2.734500000E+00, 0., 
c          2.528000000E+00,  2.528000000E+00, 0., 
c          2.262500000E+00,  2.262500000E+00, 0., 
c          1.938000000E+00,  1.938000000E+00, 0., 
c          1.554500000E+00,  1.554500000E+00, 0., 
c          1.112000000E+00,  1.112000000E+00, 0., 
c          6.105000000E-01,  6.105000000E-01, 0.,  
c          5.000000000E-02,  5.000000000E-02, 0., 
c&end

      do i=1,nbulka
         do k=1,ndens
            rho=1.d0*(k-1)/dfloat(ndens-1)
            temp1(k,i)=(ate0(i)-ateb(i))*(1-rho**rn1te(i))**rn2te(i)
     .                 +ateb(i)
         enddo
      enddo

!/tpoptab/
!--------------------------------------------------------------------------
! Tpop=T_perp/T_parallel profiles (table data, case: idens=1) 
!      tpop1(ndens,ncomp)
!--------------------------------------------------------------------------
! It this namelist we must set electron tpop1(ndens,1) and all ion
! species tpop1(ndens,i)  {i=2:nbulk}
! Constant radial step is assumed.
! The first line (k=1, i=1:nbulk) tpop1(1,i) is for rho=0
! The last line (k=ndens, i=1:nbulk) tpop1(ndens,i) is for rho=1
! The example for   ndens=5, nbulk=3, ncomp=4
!               
!         electron-1 ion-2  ion-nbulk 
! prof=      1.,     1.2,   1.3,       
!            1.,     1.2,   2.3,       
!            1.,     1.2,   3.3,       
!            1.,     4.2,   4.3,       
!            1.,     5.2,   5.3,       
!
! In all cases tpoptab should has nbulk columns with Tpop
! profiles for all nbulk plasma components.
! ------------------------------------------------------------------------
! Here array prof(ncomp,ndens) was used for convenience in namelist.
! tpop1(k,i)=prof(i,k) k=1:ndens, i=1:nbulk
! ------------------------------------------------------------------------
c&tpoptab
c prof=     1.0,  1.0, 1.0,   
c          1.0,  1.0,  1.0,   
c          1.0,  1.0,  1.0,   
c          1.0,  1.0,  1.0,   
c          1.0,  1.0,  1.0,   
c          1.0,  1.0,  1.0,   
c          1.0,  1.0,  1.0,   
c          1.0,  1.0,  1.0,   
c          1.0,  1.0,  1.0,   
c          1.0,  1.0,  1.0,   
c          1.0,  1.0,  1.0,   
c &end
      do i=1,nbulka
         do k=1,ndens
            rho=1.d0*(k-1)/dfloat(ndens-1)
            tpop1(k,i)=(tp0(i)-tpb(i))*(1-rho**rn1tp(i))**rn2tp(i)
     .                 +tpb(i)
         enddo
      enddo
c &end



!/zeftab/
!--------------------------------------------------------------------------
! Zeff profiles (table data, case: idens=1)	zeff1(ndens)
!--------------------------------------------------------------------------
! It this namelist we must set zeff1(ndens)
! Constant radial step is assumed.
! The first value  zeff1(1) is for rho=0
! The last value zeff2(ndens) is for rho=1
! The example for  ndens=5
! zeff1= 1., 1., 1., 1. 1.
!--------------------------------------------------------------------------
c&zeftab
c zeff1= 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
c&end
      do k=1,ndens
	 rho=1.d0*(k-1)/dfloat(ndens-1)
         zeff1(k)=(zeff0-zeffb)*(1-rho**rn1zeff)**rn2zeff+zeffb
      enddo
! &end

!/dentab_nonuniform_line/
!--------------------------------------------------------------------------
! density profiles at arbitrary nonuniform radial mesh
! given in line form (idens=1): dens1(ndens,nbulk)
! i.e., density profile as rows of values, for each species, 1:nbulk.
!  namelist /dentab_nonuniform_line/nj_tab,prof_2d,radii_2d
!
! integer nj_tab(nbulka): i=1,nbulk1 nj_tab(i) is the number of grid points
!                         for 'i' specie density 
!                         It should be nj_tab(i).le.ndensa
!
! real*8 prof_2d(ndensa,nbulka) : density profiles
!
! real*8 radii_2d(ndensa,nbulka): small radius meshes used for density profiles
!--------------------------------------------------------------------------
! ndensa (a parameter in param.i) is max number of points in the 
!        small radius direction.
! nbulka (a parameter in param.i) is a maximal number of plasma components.
! Input of profiles is set up so spline profiles can be input in tables
! of size specified through namelist variables ndens and nbulk.
!
! nbulk1           is a number of densities.
!                  nbulk1 depends on nbulk and izeff model (nbulk1.le. nbulk)
!                  See description of /dentab/ namelist
! nbulk (variable)    is number of plasma components must be: nbulk.le.ncompa
!                     (first component is for electrons)
! nbulk1 is number of densities components which should be given
! nbulk1 was calculated in dinit_mr subroutine ,
! The fragment of dinit_mr is given here to understand the
! nbulk1 value
!---------------------------------------------------------------------------
!Example for nbulk1=3
! /dentab_nonuniform_line/
!nj_tab(1)=3
!nj_tab(2)=5       
!nj_tab(3)=7
!
!prof_2d(1,1)=1.d0,0.5d0,0.001d0,
!prof_2d(1,2)=1.d0,0.7d0,0.5d0,0.25.d0,0.001d0,
!prof_2d(1,3)=1.d0,0.9.d0,0.7d0,0.5d0,0.25.d0,0.12.d0,0.001d0,
!
!radii_2d(1,1)=0.d0,0.5d0,1.d0,
!radii_2d(1,2)=0.d0,0.2d0,0.5d0,0.7.d0,1.d0,
!radii_2d(1,2)=0.d0,0.2d0,0.3.d0,0.5d0,0.7d0,0.9.d0,1.d0
!---------------------------------------------------------------------------
c &dentab_nonuniform_line
      do i=1,nbulka
         nj_tab_dens1(i)=ndens
         do k=1,ndens 
            rho=1.d0*(k-1)/dfloat(ndens-1)
            dens1(k,i)=(dense0(i)-denseb(i))*(1-rho**rn1de(i))**rn2de(i)
     .                 +denseb(i)
         enddo
      enddo

c     prof_2d(k,i)=
c     radii_2d(k,i)=
c&end

!/temtab_nonuniform_line/
!--------------------------------------------------------------------------
! temperature profiles at arbitrary nonuniform radial mesh
! given in line form (idens=1): temp1(ndens,nbulk)
!  namelist /temtab_nonuniform_line/nj_tab,prof_2d,radii_2d
!
! integer nj_tab(nbulka): i=1,nbulk nj_tab(i) is the number of grid points
!                         for 'i' specie temperature
!                         It should be nj_tab(i).le.ndensa
!
! real*8 prof_2d(ndensa,nbulka) : temperature profiles
!
! real*8 radii_2d(ndensa,nbulka): small radius meshes used for temperature profiles
!--------------------------------------------------------------------------
! ndensa (a parameter in param.i) is max number of points in the 
!        small radius direction.
! nbulka (a parameter in param.i) is a maximal number of plasma components.
! Input of profiles is set up so spline profiles can be input in tables
! of size specified through namelist variables ndens and nbulk.
!
!---------------------------------------------------------------------------
!Example for nbulk=3
!& temtab_nonuniform_line
!nj_tab(1)=3
!nj_tab(2)=5       
!nj_tab(3)=7
!
!prof_2d(1,1)=1.d0,0.5d0,0.001d0,
!prof_2d(1,2)=1.d0,0.7d0,0.5d0,0.25.d0,0.001d0,
!prof_2d(1,3)=1.d0,0.9.d0,0.7d0,0.5d0,0.25.d0,0.12.d0,0.001d0,
!
!radii_2d(1,1)=0.d0,0.5d0,1.d0,
!radii_2d(1,2)=0.d0,0.2d0,0.5d0,0.7.d0,1.d0,
!radii_2d(1,2)=0.d0,0.2d0,0.3.d0,0.5d0,0.7d0,0.9.d0,1.d0
!&end
!---------------------------------------------------------------------------
c &temtab_nonuniform_line
      do i=1,nbulka
         nj_tab_temp1(i)=ndens
         do k=1,ndens 
            rho=1.d0*(k-1)/dfloat(ndens-1)
            temp1(k,i)=(ate0(i)-ateb(i))*(1-rho**rn1te(i))**rn2te(i)
     .                 +ateb(i)
         enddo
      enddo

c     prof_2d(k,i)=
c     radii_2d(k,i)=
c&end

!/tpoptab_nonuniform_line/
!--------------------------------------------------------------------------
! Tpop=T_perp/T_parallel profiles
!-----------------------------------------------------------------------
! tpop profiles at arbitrary nonuniform radial mesh
! given in line form (idens=1): tpop1(ndens,nbulk)
!  namelist /tpoptab_nonuniform_line/nj_tab,prof_2d,radii_2d
!
! integer nj_tab(nbulka): i=1,nbulk nj_tab(i) is the number of grid points
!                         for 'i' specie temperature
!                         It should be nj_tab(i).le.ndensa
!
! real*8 prof_2d(ndensa,nbulka) : tpop profiles
!
! real*8 radii_2d(ndensa,nbulka): small radius meshes used for tpop profiles
!--------------------------------------------------------------------------
! It this namelist we must set electron tpop1(ndens,1) and all ion
! species tpop1(ndens,i)  {i=2:nbulk}
! Constant radial step is assumed.
! The first line (k=1, i=1:nbulk) tpop1(1,i) is for rho=0
! The last line (k=ndens, i=1:nbulk) tpop1(ndens,i) is for rho=1
! The example for   ndens=5, nbulk=3, nbulka=4
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
c&tpoptab_nonuniform_line
      do i=1,nbulka
         nj_tab_tpop1(i)=ndens
         do k=1,ndens 
            rho=1.d0*(k-1)/dfloat(ndens-1)
            tpop1(k,i)=(tp0(i)-tpb(i))*(1-rho**rn1tp(i))**rn2tp(i)
     .                 +tpb(i)
            write(*,*)'default_in k,i,tpop1(k,i)',k,i,tpop1(k,i)
         enddo
      enddo
      write(*,*)'nj_tab_tpop1',nj_tab_tpop1
c     prof_2d(k,i)=
c     radii_2d(k,i)=
c&end

!/vflowtab_nonuniform_line/
!--------------------------------------------------------------------------
! vflow  profiles at arbitrary nonuniform radial mesh
! given in line form (idens=1): vflow11(ndens,nbulk)
!  namelist /vflowtab_nonuniform_line/nj_tab,prof_2d,radii_2d
!
! integer nj_tab(nbulka): i=1,nbulk nj_tab(i) is the number of grid points
!                         for 'i' specie vflow
!                         It should be nj_tab(i).le.ndensa
!
! real*8 prof_2d(ndensa,nbulka) : vflow profiles
!
! real*8 radii_2d(ndensa,nbulka): small radius meshes used for vflow profiles
!--------------------------------------------------------------------------
! Constant radial step is assumed.
! The first line (k=1, i=1:nbulk) vflow1(1,i) is for rho=0
! The last line (k=ndens, i=1:nbulk) vflow1(ndens,i) is for rho=1
!---------------------------------------------------------------------------
c &vflowtab_nonuniform_line
      do i=1,nbulka    
         nj_tab_vflow1(i)=ndens  
         do k=1,ndens
            rho=1.d0*(k-1)/dfloat(ndens-1)
            vflow1(k,i)=(vfl0(i)-vflb(i))*(1-rho**rn1vfl(i))**rn2vfl(i)
     .                 +vflb(i)
         enddo
      enddo

c     prof_2d(k,i)=
c     radii_2d(k,i)=

c&end!

!/zeftab_nonuniform_line/
!--------------------------------------------------------------------------
! zeff profile at arbitrary nonuniform radial mesh
! given in line form (idens=1): zeff1(ndens,1)
!  namelist /temtab_nonuniform_line/nj_tab,prof_2d,radii_2d
!
! integer nj_tab(1): i=1,nbulk nj_tab(1) is the number of grid points
!                         for zeff 
!                         It should be nj_tab(1).le.ndensa
!
! real*8 prof_2d(ndensa,1) : zeff profile
!
! real*8 radii_2d(ndensa,1): small radius meshes used for zeff profile
!-----------------------------------------------------------------------
! It this namelist we must set zeff1(ndens)
! Constant radial step is assumed.
! The first value  zeff1(1) is for rho=0
! The last value zeff1(ndens) is for rho=1
!---------------------------------------------------------------------------
!Example
! /zeftab_nonuniform_line/
!nj_tab(1)=3
!
!prof_2d(1,1)=1.d0,01.d0,1.d0,
!
!radii_2d(1,1)=0.d0,0.5d0,1.d0,
!---------------------------------------------------------------------------
c&zeftab_nonuniform_line
      nj_tab(1)=ndens
      do k=1,ndens
	 rho=1.d0*(k-1)/dfloat(ndens-1)
         zeff1(k)=(zeff0-zeffb)*(1-rho**rn1zeff)**rn2zeff+zeffb
      enddo
c     prof_2d(i,1)=
c     radii_2d(i,1)=
c&end


c/read_diskf/ 
!   i_diskf=0 usage analytic Maxwellian electron distribution 
!   i_diskf=1 read in the file diskf
!   i_diskf=2 read in the file netcdfnm.nc 
!   i_diskf=3 analytic calculation of the non-Maxwellian distribution:
!             f(p,theta,rho)=ne(rho)*[(1-rtail-rhot-rbeam)*f_max(T(rho))
!                                     +rtail*f_tail+rhot*f_hot+rbeam*f_beam]
!   i_diskf=4 analytic +3D spline calculation of continuous non-Maxwellian
!              distribution with three temperatures in  three energy ranges.
!             Uses 3D spline calculation of distributions on a grid.
!   i_diskf=5 Fully analytic velocity space  +1D  spline vs radois for
!             density calculation of continuous non-Maxwellian distribution
!             with three temperatures in  three energy ranges.
!             Generally, much faster than i_diskf=4 approach.
!------------------------------------------------------
!   the data for analytic non-Maxwellian electron distribution
!    
!   jx   - the number of used normalized momentum/ mesh points
!   lrz  - the number of used radial mesh points
!   iym  - the number of used pitch-angle mesh points 
!          (here the same at each radius)
!   ngen - the number of plasma species (here we use only electron specie with 
!          the number of specie k=1) 
!   jxa,iya,lrza,ngena ! the max values for jx,iym,lrz,ngen 
!   rtem0 - ratio tem0/electron_temperature(rho=0)
!           tem0 is the max energy for the momentum normalization (KeV) 
!-----tail parameters  (i_diskf=3)
!     f_tail=H(rho,rt1,rt2)*f_rel_Maxw(ttail),
!     H(x,x1,x2) is the box function. H=1 for x1<x<x2 otherwise H=0  
!   r1t,r2t    small normalized radii for the tail localization  
!   rtail      the relation the tail density to the total density
!   ttail      tail temperature (KeV)!tail temperature (KeV)
!-----hot parameters
!     f_hot=H(rho,rh1,rh2)*H(epar,hotmnpar,hotmxpar)*H(eper,hotmnper,hotmxper)*
!     *(p_per/mc)**hotexp*exp{-mu(thoppar)(p_par/m_ec)**2-mu(thopper)(p_per/m_ec)**2}.
!     Here mu(T)=m_e*c**2/T  
!   r1h,r2h             - small normalized radii for the hot localization
!   rhot                - the relation of hot hot density to the total density
!   thotpar,thotper     - parallel and perpendicular hot temperatures (KeV)
!   hotmnpar,hotmxpar   - the boundaries of the parallel energy box(KeV)
!                         hotmnpar < epar < hotmxpar  
!   hotmnper,hotmxper   - the boundaries of the perpendicular energy box(KeV)
!                         hotmnper < eper < hotmxper  
!   hotexp              - the degree of the perpendicular momentum: (p_per/mc)**hotexp
!-----beam parameters  (i_diskf=3)
!     f_beam=H(rho,rb1,rb2)*exp{-0.5*mu(tbeam)*
!              [(p_par-p_beam_par)**2+(p_per-p_beam_per)**2]/(m_e*c)**2}
!     Here
!          (p_beam /m_e*c)**2=ebeam**2/(m_e**2*c**4)-1
!           p_beam_par=p_beam*cos(thbeam)
!           p_beam_per=p_beam*sin(thbeam)
!   r1b,r2b      - small normalized radii for the beam localization
!   rbeam        - the relation of the beam density to the total density
!   ebeam        - beam energy (KeV)
!   thbeam       - beam pitch angle (0=<degree=<180) 
!   tbeam        - beam temperature (KeV)
!-----Three temperature case  (i_diskf=4)
!   rvtovte1,rvtovte2 = ratio of momentum-per-mass (electrons) to on-axis
!                       thermal velocity vte0= sqrt(Te/me), defining the
!                       three velocity ranges for the temperatures.
!                       defaults=1.e6,1.e6 [i.e., effectively infinity]
!   rtemp1, rtemp2, rtemp3 = ratios of temperatures in each of the
!                       three velocity (energy) bins to the radially
!                       local temperature.
!   In summary:  The three momentum-per-mass bins are [0.,rvtovte1*vte0],
!                [rvtovte1*vte0,rvtovte2*vte0], and [rvtovte2*vte0,infinity].
!                These bins are constant as a function of radius.
!                The temperatures in each bin are given by rtemp[1-3]
!                and vary as a function of radius as the bulk temperature.
!--------------------------------------------------------------------
c &read_diskf
      i_diskf=0
      netcdfnm='netcdfnm.nc'
      jx=100
      lrz=10
      iym=100  
      ngen=1  
      rtem0=10.d0
      r1t=0.d0
      r2t=1.d0      
      rtail=0.d0      
      ttail=1.d0
      r1h=0.d0
      r2h=1.d0             
      rhot=0.d0               
      thotpar=1.d0
      thotper=1.d0     
      hotmnpar=1.d0
      hotmxpar=2.d0   
      hotmnper=1.d0
      hotmxper=2.d0    
      hotexp=1.d0    
      r1b=0.d0
      r2b=1.d0     
      rbeam=0.d0        
      ebeam=1.d0       
      thbeam=30.d0       
      tbeam=1.d0      
c--------for i_diskf=4 and i_diskf=5
      rvtovte1=1.d6
      rvtovte2 =1.d6
      rtemp1=1.d0
      rtemp2=1.d0
      rtemp3=1.d0
c&end

!--------------------------------------------------
!/emission/ 
!the data for emission calculations
!   i_emission=0 no emission
!             =1 emission calculations
!
!   0<tol_emis=<1 tolerance parameter to add the new mesh point s_n
!                if in_0(n)>tol_emis*i_0 
!   nharm1=< nharm =<nharm2 gives the used EC harmonics (Not work now)
!   nfreq=<nfreqa (see param.i) is the number of frequencies
!         nfreq=1 gives detailed plots of emission from a single ray     
!         nfreq.gt.1 gives spectra covering the specified frequency
!         range
!   if nfreq=1 code will use the frequency determined by 
!              frqncy(GHz) given in namelist wave
!   freq00=<freq01 are ratios of the minimal and maximal emission
!         frequencies to electron gyro-frequency f_ce(rho=0) at the
!         magnetic axis.  Search screen dump for "magnetic axis bmag"
!         to value.  Gyrofrequency (GHz) is 28.0*<magnetic axis bmag>.
!    wallr is the wall reflection coefficient, {0=< wallr =<1}
!
!    i_rrind chooses the subroutine to calculate N_ray (default =0)
!      i_rrind=0  N_ray=N
!      i_rrind =1 from cold electron plasma using rrind   
!      i-rrind =2 from hot non-relativistic dispersion relation  
!--------------------------------------------------
!    i_r_2nd_harm=1 to calculate the major radius of the EC 2nd harmonic
!                =0 do not calculate (default =0) 
!              (to use drawemfr.in file i_r_2nd_harm=1)
!              (to use drawemf1.in file i_r_2nd_harm=0)
!---------------------------------------------------
!    i_emission_spectrum =0  do not calculate emission spectrums
!                         1  calculate emission spectrums
!    jx_kin               The number of  kinetic energy grid points
!                         for emission spectrum.
!    max_kin_energy_kev   maximal kinetic energy [KeV] of the grid 
!                         for emission spectra. 
!---------------------------------------------------
c&emission
      i_emission=0
      tol_emis=1.d0
      nharm1=1
      nharm2=1
      nfreq=1
      freq00=1.d0
      freq01=1.d0
      wallr=1.d0
      i_rrind=1
      i_r_2nd_harm=0
      i_emission_spectrum =0  
      jx_kin=50  
      max_kin_energy_kev=10.d0   
c&end

!/ox/
!------------------------------------------------------------------------
!namelist related to calculation of EC cone vertex coordinates
!for the optimal OX mode  conversion (details on this calculation
!will be given in the Genray manual)
!------------------------------------------------------------------------
! i_ox =0 /default/ do not use these calculations.
!      =1 calculations of the optimal EC cone vertex for OX conversion. 
!         The optimal direction will give optimal N_parallel
!         in OX conversion point at Xe=1.
!         In this case should have istart=1.
!         Optimal launch angles are output to ECcone_optimal.dat.
!      =2 launch the ray using EC cone vertex coordinates calculations
!         and using OX transmission procedure
!------------------------------------------------------------------------
! theta_bot(icone) < theta_top(icone) 
! are the poloidal angle boundaries (degree) at
! Xe=1 surface. They are used in genray.f to find the optimal ray
! direction at the given EC cone vertex icone=1,...,ncone
!----------------------------------------------------------
! i_ox_poloidal_max is the maximal number of the poloidal angles used
!                   in the bisection method. This method calculates the
!                   the poloidal angle theta_pol of the point M at the
!                   flux surface Xe(rho=1) : M(poloidal_angle,rho).  
!                   The ray launched from M to the plasma edge will
!                   go to the EC cone vertex.
!---------------------------------------------------------
! eps_antenna       is the accuracy with which the ray launched from the
!                   flux surface Xe(rho)=1 at the given poloidal angle
!                   reaches the given cone vertex (i_ox=1).
!                   This accuracy is calculated using the distance of 
!                   the vertex radial coordinate from the specified
!                   antenna position, 
!                   sqrt((rst_ox-rst)**2+(zst_ox-zst)**2)< eps_antenna
!----------------------------------------------------
!eps_xe   The parameter which sets the vicinity
!         of the O-mode cutoff surface.
!         If xe < (1-eps_xe) then this subroutine will 
!         creates the ray jump in small radius direction
!         and find the X mode.
!         eps_xe=1.d-2 is seted as default in dinit.f
!----------------------------------------------------
c&ox
      i_ox=0
      theta_bot(1)=0.0d0
      theta_top(1)=180.d0
      i_ox_poloidal_max=20
      eps_antenna=1.d-4
      eps_xe=1.d-2
c&end


!---------------------------------------------------------
!/adj_nml/
!---------------------------------------------------------
! INTEGER
! i_adj   =1 to use adj calculations
!        =0 do not use adj_calculastions 
! npsi0  is a number of radial points  
! nthp0  is a number of poloidal points for the integration
!        along b filed line
! nmax_chi   is the number of grid points in u_0 direction (for adj mesh)
! imax_chi   is the number of grid points in pitch angle 
! lmax_chi   is the maximal number of Legandre harmonic
!        lmax=1 corresponds to collision with a Maxwellian.
!        Usually lmax=21
! tmax   is the maximal number of time steps
!
! n_harm_adj_min   number of minimal and maximal harmonics
! n_harm_adj_max   for power and CD calculations  
!                  It should be:  n_harm_adj_min.le.n_harm_adj_max
! n_relt_intgr_adj        number of points for Power and CD integration
!                         in p_perp direction along resonce curve
! REAL*8
! t      is the temperature in units of m*u_n**2
!        If you want velocities normalized to sqrt(T/m), set t=1.
! ze     is the electron charge state. It is usually =1, 
!        but for Lorentz limit, ion charge Zi ==> infinity,
!         set  ze=0 and zi=1
! umax   is the maximal value of u_0. 
!        The grid is spaced between 0 and umax.
!        The momentum in the code is normalized at unorm=sqrt(T/m),
!        u^=u/unorm.
!        T is a temperature at the flux surface where the adj equation
!        is solving
!        So, the maximal momentum at the given flux surface is
!        u_maximal=(p/m)_maximal=unorm*umax
! dt     is time step in adj solver, usually =1
! alpha  is a weighting factor in the time stepping algorithm.
!        =0 for explicit algorithm
! rho_   governs the treatment  of the boundary at u=umax
!        for rho_=1 the second derivative from the adjoint function chi
!        d^2(chi)/dx^**2(at u=umax)=0
!        for rho_=0 the first derivative =0
! aerrmx  is some error for subadj 
! rerrmx  is some error for subadj 
! epsi_adj  the accuracy used in integration by adaptive Simpson
!           Works at i_resonance_curve_integration_method_adj=4 case only
!-------------------------------------------------------------------------------
!                  choose the integration method for power and CD calulations 
!       i_resonance_curve_integration_method_adj=1 !rectangle integration
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
!----------------------------------------------------------------
!  i_calculate_or_read_adj_function      =1 calculate chi function
!                                        =0 read chi function from
!                                           the file: adjout
!---------------------------------------------------------------
! i_chi_interpolation            to choose interpolation method
!                                for derivatives from chi function
!                                =1 using spline [default value]
!                                   [has given problem near t-p bndry]
!                                =2 linear interpolation
!                                   of numerical derivatives
!--------------------------------------------------------------
!&adj_nml
        i_adj=0
cSAP090827
c        npsi0=npsi0_a
c        nthp0=nthp0_a
        npsi0=10
        nthp0=1000
        nmax_chi=100 
        imax_chi=200
        lmax_chi=21
        tmax=300
        t=1.d0
        ze=1.d0
        umax=10d0
        dt=1.d0  
        alpha=0.55d0
        rho_= 1  
        aerrmx=1.d-5  
        rerrmx=2.d-6
        n_harm_adj_min=0
        n_harm_adj_max=1
        n_relt_intgr_adj=100
        i_resonance_curve_integration_method_adj=4
        epsi_adj=1.d-5
        i_calculate_or_read_adj_function =1 
        i_chi_interpolation=1   
!&end
!-------------------------------------------------------------
!/edge_prof_nml/  to set density profile outside LCFS
!     integer
!     &i_edge_dens_anal,  ! =0 to use sigmedgn=constant
!                         ! =1 the analytic formula for sigmedgn(theta_pol)
!                         ! =2 table data for  sigmedgn(theta_pol)
!                            [default = 0]
!---------------------------------------------------------------------
!     &n_pol_edge_dens    ! number of poloidal points of
!                         ! the poloidal mesh [See below].
!-----------------------------------------------------------------------
!   i_edge_dens_rz_mesh is used at n_wall.gt.0. 
!                       It means that wall is used
!
!                       For n_wall=0 no wall to be used,
!                       no reflection from the wall
!                       and rz_mesh will not be used
!--------------------------------------------------------------
!    &i_edge_dens_rz_mesh, ! =0 No write (density_rz,temperautre_rz, 
!                          !    flux_rz, rho_rz)
!                          !    and no other calculations related
!                          !    i_edge_dens_rz_mesh
!                          !
!                          ! =1 Write to the file dens_temp_rz_out.txt:
!                          !    density_r_z (i,j,k,0)  and temperature_r_z (i,j,k)
!                          !    (1:nxeqd_add,1:nyeqd_add,1:nbulk)
!                          !    on rz mesh for test or modification.
!                          !    Format (as given in equilib.f) is:
!                          !    write(i_unit,2) nxeqd_add,nyeqd_add,nbulk !mesh dimensions
!                          !    write(i_unit,3)(((dens_r_z_in_MKS(i,j,k),
!                          !   &                i=1,nxeqd_add),j=1,nyeqd_add),k=1,nbulk) !MKS
!                          !    write(i_unit,3)(((temp_prof(i,j,k),
!                          !   &                i=1,nxeqd_add),j=1,nyeqd_add),k=1,nbulk) !KeV
!                          ! 2  format(52x,3i5) 
!                          ! 3  format(5e16.9)
!                          ! 4  format(52x,2i5)
!                          !  
!                          !    Write to the output file: psi_rho_rz_out.txt:
!                          !    nxeqd_add,nyeqd_add, !dimensions of rz mesh
!                          !    xdimeqd,ydimeqd,     !horizontal and vertical
!                          !                         !  full-widths of the mesh [meters]
!                          !    redeqd,              !major radius of the inner edge 
!                          !                         !  of the rectangular grid.
!                          !    ymideqd,             !vertical shift of the
!                          !                         !  rectangular grid.
!                          !    flux_rz, rho_rz (on the rz mesh).
!                          !    Format is as in equilib.f:
!                          !    write(i_unit,4) nxeqd_add,nyeqd_add  !mesh dimensions
!                          !    write(i_unit,3) xdimeqd,ydimeqd,reqd,redeqd,ymideqd
!                          !    write(i_unit,3)(rr_add(i),i=1,nxeqd_add)
!                          !    write(i_unit,3)(zz_add(j),j=1,nyeqd_add)
!                          !    write(i_unit,3)((rho_r_z(i,j),i=1,nxeqd_add),j=1,nyeqd_add)
!                          !    write(i_unit,3)((psi_r_z(i,j),i=1,nxeqd_add),j=1,nyeqd_add)
!                          !
!                          !    Here density_rz profile are determined from
!                          !    i_edge_dens_anal model.
!                          !
!                          ! =2 Use  the density and temperature data from the file
!                          !      dens_temp_rz_in.txt.  NOTE:  _in.dat, 
!                          !                                   NOT _out.dat, as above.
!                          !      That is, read density_r_z and
!                          !      temperature_r_z as in above dens_temp_rz_out.txt format.
!                          !      i_edge_dens_anal value not used in this case.
!                          !      Data file psi_rho_rz_out.txt is not needed.
!                          !
!                          ! The idea of writing dens_temp_rz_out.txt with =1 option,
!                          ! is that it provides a template for modification to the file
!                          ! dens_temp_rz_in.txt for the =2 option.  The file
!                          ! psi_rho_rz_out.txt can provide help is making the density/
!                          ! temperature modification, but is not required for the
!                          ! i_edge_dens_rz_mesh=2 option.
!                          !     
!                          ! [default = 0]
!-------------------------------------------------------------------------
! For the temperature outside LCFS (rho>1) at no_reflection=1
! the code will use following formula at all i_edge_dens_anal values:
!
! temperature(i,rho)=temperature(i,rho=1)*exp(-(rho-1)/sigmedgt)
!
! Here
! i is a number of plasma specie i=1,..,nbulk,
! sigmedgt=const is normalized to the plasma radius.
! By default
! sigmedgt=0.02
!      
! For densities outside LCFS (rho>1) at no_reflection=1 the code will use
! different formula according to i_edge_dens_anal value.
!.......................................................................
! At i_edge_dens_anal=1 the code will use following formula
!
! density(i,rho)=density(i,rho=1)*exp(-(rho-1)/sigmedgn)
! Here
! i is a number of plasma specie i=1,..,nbulk,
! sigmedgn=const is normalized to the plasma radius.
! By default
! sigmedgn=0.02
!........................................................................
! At i_edge_dens_anal=1  the code will use following formula
!
! density(i,rho,theta_pol)=density(i,rho=1)*exp(-(rho-1)/sigma_edge_n(theta_pol))
!
! The function  sigma_edg_n(theta_pol) is normalized to the plasma radius.
!
! sigma_edge_n(theta_pol_radian)=sigma_edgen_0+
!                             +del1*exp1(theta_pol_radian)+del2*exp2(theta_pol_radian)
!
! exp1(theta_pol_radian)=exp(-((theta_pol_radian-theta_pol_edge_1_radian)/
!                        sigma_theta_pol_edge_1_radian)**2)
! exp2(theta_pol_radian)=exp(-((theta_pol_radian-theta_pol_edge_2_radian)/
!                        sigma_theta_pol_edge_2_radian)**2)
! del1=sigma_edgen_1-sigma_edgen_0
! del2=sigma_edgen_2-sigma_edgen_0
!
! Input angles in the namelist should be set in degree.
!
! So, sigma_edge_n is equal to
! 1) =sigma_edgen_0 at poloidal angles which are far from given poloidal angles
!     theta_pol_edge_1 and theta_pol_edge_2
! 2) =sigma_edgen_1  poloidal angles near theta_pol_edge_1
! 3) =sigma_edgen_2  poloidal angles near theta_pol_edge_2
!.........................................................................
! At i_edge_dens_anal=2  the code will use following formula
!
! density(i,rho,theta_pol)=density(i,rho=1)*exp(-(rho-1)/sigma_edge_n(theta_pol))
!
! The function sigma_edg_n(theta_pol) is normalized to the plasma radius.
! sigma_edg_n(theta_pol) is a spilne approximation at the given tables:
!
! theta_pol_edge_dens_ar_degree(1:n_pol_edge_dens) poloidal mesh 
!
! sigmedgn_ar(1:n_pol_edge_dens) normalized to small radius and
!                               set in poloidal mesh points 
!------------------------------------------------------------------------
! Minimal values of density and temperature outside LCFS
!
! dens_min_edge                   !minimum edge density
!                                 !10**13/cm**3 for all plasma species
!                                 !Applies for all i_edge_dens_anal values
! temp_min_edge,                  !minimal edge temperature
!                                 ![KeV] for all plasma species
!
! Outside LCFS the code uses dens_min_edge to sets
! 1) the minimal electron density equal to "dens_min_edge"
! 2) for ions the minimal ion density will be different for each ion
!    specie i=2,,....,nbulk to create the plasma charge equal to Z_eff 
!    at LCFS (rho=1) 
!    dens_min(i)=dens_min_edge*ratio(i)
!    Here ratio(i)=ion density(i,rho=1)/electron density(rho=1) 
!    is the ratio of the ion density(i) at rho=1
!    to the electron density at rho=1
!     
!-------------------------------------------------------------------------
!     Data to set 
!     &n_pol_edge_dens    ! number of poloidal points (integer) of
!                         ! the poloidal mesh for i_edge_dens_anal=2!
!
!     theta_pol_edge_dens_ar_degree, !poloidal angle mesh [0=<degree=<360]
!     It is assumed that:
!      theta_pol_edge_dens_ar_degree(1)=0
!      theta_pol_edge_dens_ar_degree(i) < theta_pol_edge_dens_ar_degree(i+1)
!      theta_pol_edge_dens_ar_degree(n_pol_edge_dens )= 360 degree
!
!    For example
!     theta_pol =90 at the top of the poloidal cross-section
!     theta_pol=270 at the bottom of the poloidal cross-section
!
!                                     !to set sigmedgn_ar
!     sigmedgn_ar                     !exponential density fall-off distance
!                                     !outside LCFS starting at rho=1 density
!                                     !Normalized to plasma radius 
!
!     theta_pol_edge_1_degree         !for analytical formula of sigma_edge_n
!     theta_pol_edge_2_degree            
!     sigma_theta_pol_edge_1_degree,
!     sigma_theta_pol_edge_2_degree
!     sigma_edgen_0  !Normalized to plasma radius
!     sigma_edgen_1  !Normalized to plasma radius
!     sigma_edgen_2  !Normalized to plasma radius
!
!     The density profile outside LCFS has the form
!
!     dens_rho_theta=densedge*dexp(-(rho_small-1.d0)/sigma_edge_n)
!

!     [densedge is the density for the given plasma specie at LCFS (rho=1).
!       Not namelist, but input with the plasma profiles within the LCFS.]
!
!     sigma_edge_n depends on the poloidal angle theta_pol.
!      
!     If i_edge_dens_anal=1 then the analytic formula for sigma_edge_n is used:
!
!     sigma_edge_n(theta_pol)=
!     (sigma_edgen_1-sigma_edgen_0)*     
!     exp(-((theta_pol_radian-theta_pol_edge_1_radian)/sigma_theta_pol_edge_1_radian)**2)+
!     (sigma_edgen_2-sigma_edgen_0)* 
!     exp(-((theta_pol_radian-theta_pol_edge_2_radian)/sigma_theta_pol_edge_2_radian)**2)
!
!     If i_edge_dens_anal=2 then the table will be used to set sigma_edge_n(i)
!     at the poloidal mesh theta_pol_edge_dens_ar_degree(i) ,i=1,...,n_pol_edge_dens
!
!--------------------------------------------------------------------------------
!
!     The code can create exponential density fall outside LCFS near chamber
!     wall and limiters.
!     This density fall can be used for natural ray reflection
!     from the chamber wall and limiters.
!
!     The chamber wall has the same poloidal crossection for all toroidal angles.
!     So, the chamber wall is determined by its poloidal boundary.
!
!     The code can use several limiters. The number of limiters max_limiters
!     is set in the namelist /tokamak/. Each limiter m=1:max_limiter has the poloidal boundary
!     and it lokalised in the toroidal direction by toroidal angles :
!     phi_limiter(1,m) < phi < phi_limiter(2,m)
!     which are set in the namelist /tokamak/.
!
!     The poloidal limiter boundary consists from the line L_limiter which
!     is close to LCFS and two horizontal lines. 
!     The top horizontal line connects the top point of L_limiter line with
!      the chamber wall. 
!     The bottom horizontal line connects the bottom point of L_limiter line
!      with the chamber wall. 
!       
!     So, each limiter has the poloidal boundary at RZ plane and vertical
!     toroidal plane bounaris perpendicular to the toroidal direction.
!                               
!................density wall-limiters fall at the poloidal plane..............
!
!     The poloidal character length of the chamber wall and limiters density 
!     fall at the poloidal plane is set by  the input variable:
!     sigma_wall_n  -[meters] exponential density fall off polidal distance
!                   near chamber wall outside LCFS  
!
!     To get the ray reflection close in to the wall or limiters sigma_wall_n
!     should be small.
!     By default sigma_wall_n=3.d-3 [meter]. 
!
!     The code uses spline at RZ poloidal mesh to create density with density
!     fall near chamber wall and each limiters.
!     This RZ mesh rr_add(nxeqd_add),zz_add(nyeqd_add) should have
!     small steps in R and Z directions to resolve
!     the density fall at the given character length sigma_wall_n
!     It means that RZ mesh step should me smaller than sigma_wall_n.
!
!     The numbers of RZ mesh points are set by
!       
!     nxeqd_add       =< nxeqd_add_a  number of points at RZ mesh at  R-direction 
!     nyeqd_add       =< nyeqd_add_a  number of points at RZ mesh at  Z-direction 
!                        
!     The step of rr_add mesh is xdimeqd/(nxeqd_add-1)
!     The step of zz_add mesh is ydimeqd/(nyeqd_add-1)
!     Here xdimeqd, ydimeqd  are the horizontal(R-direction) and vertical
!     (Z-direction) full-widths [meters] of the  rectangle given it the input
!     equilib.dat eqdsk file
!
!     So, for the resolution of the wall-limiters density falls at the
!     poloidal plane it needs
!     xdimeqd/(nxeqd_add-1) =< sigma_wall_n
!     ydimeqd/(nyeqd_add-1) =< sigma_wall_n
!
!      sigma_wall_n    [meter] exponential density fall off poloidal distance
!                     near the chamber wall-limiters outside LCFS       
!     nxeqd_add       number of points at RZ mesh at  R-direction 
!     nyeqd_add       number of points at RZ mesh at  Z-direction 
!                     They are used for creation density_r_z   
!
!-----The calculation of the poloidal wall-limiters density fall
!
!     For each RZ point of the poloidal mesh P(rr_add(i),zz_add(j)) the code
!     calculates the minimal poloidal distance distance_to_wall(i,j,m)
!     between mesh point P(i,j) and
!     1) all chamber wall points  
!        M(r_wall_add(k,m),z_wall_add(k,m)) for m=0
!     2) all limiter points points
!         M(r_wall_add(k,m),z_wall_add(k,m)) for each
!        limiter m=1,...,max_limiters.
!
!
!     Then the code at each point P(i,j) outside LCFS calculates factor 
!     for the chamber wall m=0 
!
!     factor=1.d0-dexp(-(distance_to_wall(i,j,m=0)/sigma_wall_n)**2) 
!
!     which is equal to zero at the wall and it is equal to unit
!     far from the wall.
!     
!     Then for each limiter m=1,max_limiters and each point P(i,j) outside
!     LCFS the code calculates factor_limiter 
!     factor_limiter=1.d0-
!     &            dexp(-(distance_to_wall(i,j,m)/sigma_wall_n)**2)
!
!    Then for each limiter m=1,max_limiters 
!    if (factor_limiter.lt.factor) factor=factor_limiter
!
!    The density in each P(i,j) point outside LCFS is
!    density_r_z(i,j,k,m)= dens_loc(zz_add(i),rr_add(j), k)*factor
!    Here 
!    k=1,nbulk is a number of plasma specie           
!    dens_loc(zz_add(i),rr_add(j), k) is a density value calculated
!    at the point P(i,j) according to the given i_edge_dens_anal
!    without chamber-limiter density fall
!
!    After that the code calculates spline coefficients for 2D approximation
!    of the density at RZ mesh (zz_add,rr_add) using the found density with
!    density wall-limiter fall dens_loc(zz_add(i),rr_add(j), k).
!
!    These spline coefficients will be used to find density with
!     wall-limiters fall.
!
!    If the toroidal angle phi of the point Q is inside the given limiter
!    vertical boundaries with number 1=< m =<max_limiters
!    phi_limiter(1,m)< phi <phi_limiter(2,m) 
!    In this case the density in this point will be calculated using 
!    spline coefficients for found limiter m:   
!    density_r_z(i,j,k,m),
!    density_r_z_rr(i,j,k,m),density_r_z_zz(i,j,k,m),
!    density_r_z_rrzz(i,j,k,m).
!
!    If the toroidal angle phi of the point Q is outside of any 
!    given limiters then the density in this point will be calculated using 
!    spline coefficients for m=0:   
!    density_r_z(i,j,k,m=0),
!    density_r_z_rr(i,j,k,m=0),density_r_z_zz(i,j,k,m=0),
!    density_r_z_rrzz(i,j,k,m=0),
!    In this case the density fall is near the chamber wall only.
!
!.....density limiters fall near the vertical plane limiter boundaries
!      perpendicular to the toroidal direction........................... 
!
!     The toroidal character angle of limiters density fall at the vertical
!     limiter boundaries is:
!     sigma_lim_toroidal_degree [degree]  [degree] is the density fall
!                                         off toroidal angle
!                                         at the vertical limiter
!                                          boundaries
!
!     The code checks if the point P is inside one of limiter with number m
!     then the code founds the vertical limiter boundary nearest to the given
!     point P: nearest_lim_boundary=1 or =2. 
!     Then the code calculates factor 
!
!     factor=dexp(-(phi-phi_limiter(nearest_lim_boundary,m)/
!                   sigma_lim_toroidal_radian)**2)
!
!     Then the code calculates density at the point P according to
!     the given i_edge_dens_anal without chamber-limiter
!     density fall: dens=density_r_z_i_m(z,r,k,m)
!
!     After that the density with limiter fall is
!     dens_lim=dens*factor
!------------------------------------------------------------------------
!&edge_prof_nml
      sigmedgt =0.02d0  
      i_edge_dens_anal=0  
      i_edge_dens_rz_mesh=0
!-----_edge_dens_anal=0 case-------------------
      sigmedgn =0.02d0          
!      sigmedgn =0.1d0
!-----i_edge_dens_anal=1 case-------------------       `      
      sigma_edgen_0=0.02d0              !Normalized to plasma radius
      sigma_edgen_1=sigma_edgen_0       !Normalized to plasma radius
      sigma_edgen_2=sigma_edgen_0       !Normalized to plasma radius
      theta_pol_edge_1_degree =90.d0
      theta_pol_edge_2_degree =270.d0   
      sigma_theta_pol_edge_1_degree=90.d0
      sigma_theta_pol_edge_2_degree=90.d0 
   
!-----i_edge_dens_anal=2  case-------------------
      n_pol_edge_dens=n_pol_edge_dens_a
      do i=1,n_pol_edge_dens
         theta_pol_edge_dens_ar_degree(i)=360.d0/(n_pol_edge_dens-1)
     &                                    *(i-1)
         sigmedgn_ar(i)=0.02d0 
      enddo    
!-----minimal density and temperature outside LCFS------------------
      dens_min_edge=1.d-6
      temp_min_edge=1.d-3
!-----density fall near the chamber wall and limiters---------------- 
      sigma_wall_n=1.d-3  ![m]
      nxeqd_add=1000
      nyeqd_add=1000 
      sigma_lim_toroidal_degree=1.d-3  
!&end
!--------------------------------------------------------------------------
!     data for namelist /lsc_approach_nml/  
!     to calculate LH absorption and CD  for nonmaxwellian
!     electron distribution function
!     and given DC electric field radial profile
!-------------------------------------------------------------------------
!     i_lsc_approach     !=0 do not use lss_approach
!                        !=1 to use lsc_approach
!-------------------------------------------------------------------------
!     nv_lsc_ql_lh      !number of 1the parallel velocity uniform mesh points 
!                       !used for QL diffusion coefficient and
!                       !electron distribution function
!                       !v_par_mesh_lsc(-nv_lsc_ql_lh:nv_lsc_ql_lh)
!---------------------------------------------------------------------
!     v_par_lsc _dev_vte_max !maximal of the parallel velocity mesh
!                            !normalized to thermal electron velocity  
!                            !v_par_mesh_lsc(nv_lsc_ql_lh)=
!                            !                   v_par_lsc_dev_vte_max
!---------------------------------------------------------------------
!    n_psi_TSC               !number of central bin radial mesh points for TSC
!---------------------------------------------------------------------- 
!    EdcTSC_min              !minimal value of  abs(E_DC) [V/m] 
!                            !if  abs(EdcTSC_1D(i)).lt.EdcTSC_min then
!                            !EdcTSC_1D(i)=EdcTSC_min    
!-----------------------------------------------------------------------
!    icount_iter_max         ! max number of iteratins to get QL diffusion
!                            ! coefficients. 
!                            ! It is used in subroutine LSC_power_iterations
!------------------------------------------------------------------------
!    eps_diff_d_ql           ! accuracy of iteratins to get QL diffusion
!                            ! coefficients:
!                            ! diff_d_q < eps_diff_d_ql.   
!                            ! diff_d_q=||d_ql_lsc_lh_ar_old- d_ql_lsc_lh_ar||
!                            ! It is used in subroutine LSC_power_iterations
!------------------------------------------------------------------------
!   n_power_steps            !number of iterations in initial power
!                            !initial power at each ray will be
!                            !(i_power_step*1.d0/n_power_steps)*powinilh(iray)
!------------------------------------------------------------------------
!  n_dim_radii_f_plots       ! number of radii points from TSC mesh
!                            ! rho_bin_center(1:n_psi_TSC)
!                            ! at which obtained distribution function 
!                            ! ln(f_e(v_parallel)) will be plotted 
!                            ! into the plot.ps file
!                            ! It should be  n_dim_radii_f_plots.le.n_psi_TSC 
!                            ! 
!                            ! If n_dim_radii_f_plots=0 then
!                            ! the codfe will not plot f_e
!----------------------------------------------------------------------------
! n_radii_f_plots_ar(1:n_dim_radii_f_plots) !numbers of radii points
!                                           ! from TSC mesh
!                                           ! rho_bin_center(1:n_psi_TSC)
!                                           ! at which obtained distribution
!                                           ! function ln(f_e(v_parallel))
!                                           ! will be plotted 
!                                           ! into the plot.ps file
!----------------------------------------------------------------------------
!& lsc_approach_nml
      i_lsc_approach=0  !do not use lsc_approach!
      nv_lsc_ql_lh=100  
      v_par_lsc_dev_vte_max=10
      n_psi_TSC=NR-1  
      EdcTSC_min=1.d-5    ![V/m]
      icount_iter_max=20 
      eps_diff_d_ql=1.d-6
      n_power_steps=1  
c      n_dim_radii_f_plots=1               !for one outpur point
c      n_radii_f_plots_ar(1)=n_psi_TSC/2
      n_dim_radii_f_plots=2                !for two output points
      n_radii_f_plots_ar(1)=n_psi_TSC/2    !for two output points
      n_radii_f_plots_ar(2)=2*n_psi_TSC/3  !for two output points
!&end
!--------------------------------------------------------------------------
! /EdcTSCtab/
!    DC electrical field (along B field) profle at [v/m] 
!    at uniform grid with points given at the centers of radial bins
!    rho(i)=(i-0.5)/n_psi_TSC  i=1,..., n_psi_TSC
!
!    Example fo case  n_psi_TSC=3
!
!&EdcTSCtab
! prof1_uniform_TSC=   1.d-2   1.d-2   1.d-2      
!&end
      do i=1,n_psi_TSC
         EdcTSC_1D(i)=EdcTSC_min
      enddo
!--------------------------------------------------------------------------
! /JparTSCtab/ [A/m**2]
!    Parallel total current (along B field) profile [A/m**2]
!    at uniform grid with points given at the centers of radial bins
!    rho(i)=(i-0.5)/n_psi_TSC  i=1,..., n_psi_TSC
c
c    Example fo case  n_psi_TSC=3
c&JparTSCtab
c prof1_uniform_TSC=   1.d-2   1.d-2   1.d-2      
c&end
      do i=1,n_psi_TSC
         JparTSC_1D(i)=0.d0      
      enddo
!---------------------------------------------------------------------------

      return
      end          !!!!!!!!!!!!default_in!!!!!!!!!!!!!!!!

      

      subroutine pack12S(a,nmax,k,p,b,kmax,pmax)
c     put 1D array a(1...nmax)to 2d array b(1..k,1...p),
c       where b is dimensioned b(kmax,pmax)     
      implicit none
c-----input
      integer nmax  !maximal value of the dimension a
      real*8  a(nmax)
      integer k,p       !dimensions of 2d array of b
      integer kmax,pmax !dimensions of 2d array b, equal
                        !maximal values of the dimensions k,p
c-----output
      real*8 b(kmax,pmax)
c-----locals
      integer ik,ip,j

      do ip=1,p
         do ik=1,k
             j=(ip-1)*k+ik
             b(ik,ip)=a(j)
         enddo
      enddo

      return
      end 
          

      subroutine put_to_uniform_mesh(nmax,nmaxa,x,y,ct,
     &n_out,x_out,y_out)
c-----------------------------------------------------------------------
c     put the input array data y_out from non-uniform grid x,y
c     to uniform grid x_out using spline from onetwo code
c---------------------------------------------------------------------
      implicit none

c-----input
      integer nmax,   ! the number of input non-uniform grid points 
     &nmaxa          ! the maximal number of input non-uniform grid
                      ! and output uniform grid 
   
      real*8 x(nmaxa),! input non-uniform grid
     &y(nmaxa)        ! input function tab at non-uniform grid
      integer n_out   ! the number of output uniform grid points
      real*8 x_out(nmaxa),! output uniform grid
     &y_out(nmaxa)        ! output function tab at uniform grid
      real*8 ct(3,nmaxa-1) !work array

      
      call intrp_adp (nmax,x,y,ct,n_out,x_out,y_out)
     
      return
      end
        
      subroutine pack12S_nonregular(a,nmax,k,p,b,kmax,pmax)
c     put 1D array a(1...nmax)to 2d array b(1..k,1...p),
c       where b is dimensioned b(kmax,pmax)     
      implicit none
c-----input
      integer nmax  !maximal value of the dimension a
      real*8  a(nmax)
      
      integer kmax,pmax !dimensions of 2d array b, equal
                        !maximal values of the dimensions k,p  
      integer k,p(kmax)       !dimensions of 2d array of b
c-----output
      real*8 b(kmax,pmax)
c-----locals
      integer ik,ip,j,pmax_loc

      pmax_loc=0
      do ik=1,k
        if (p(ik).gt.pmax_loc) pmax_loc=p(ik)
      enddo
      write(*,*)'pmax_loc',pmax_loc
      write(*,*)'k',k
      j=0
      do ip=1,pmax_loc
         write(*,*)'ip',ip
         do ik=1,k
            write(*,*)'ip,ik,p(ik)',ip,ik,p(ik)
            if(ip.le.p(ik)) then
              j=j+1
              write(*,*)'ip,k,j',ip,k,j
              b(ik,ip)=a(j)
             endif
         enddo
      enddo


      return
      end 

      subroutine transform_genray_in_to_dat
c-------------------------------------------------------------------------
c     conversion input data from genray.in (MKSA) to genray.dat 
c     (mixed units) format.
c-------------------------------------------------------------------------
      include 'param.i'
c----------------------------------------
c     declaration and storage of all namelists variables
c-----------------------------------
      include 'cone_nml.i'
c      include 'emissa_nml.i'
      include 'grill_nml.i'
c      include 'ions_nml.i'
      include 'one_nml.i'
c      include 'onetwo_nml.i'
c      include 'output_nml.i'
c      include 'scatnper_nml.i'    
      include 'six_nml.i'  
c      include 'adj_nml.i'
cSAP090209
       include 'edge_prof_nml.i'
c--------------------------------------------------------------
c     declaration of namelist variables which
c     are used in subroutine dinit_mr and
c     are not used in any common block
c--------------------------------------------------------------  
c      include 'dinit_nml.i' 

c-----local
      integer kode,i,k
  
c-----transformation data /wave/ from genray.in to genray.dat form
      frqncy=frqncy*1.d-9    !from HZ to GHZ

      if (istart.eq.1) then
c---------------------------------------------------------
c        start point is outside the plasma, ECR case
c        the reading of the data for EC cone
c----------------------------------------------------------  
c--------transformation data /eccone/ from genray.in to genray.dat form 
         do i=1,ncone 
            powtot(i)=powtot(i)*1.d-6     !from Watt to MWatt
         enddo 
      else
        !istart= 2 or 3     
c-------------------------------------------------
c        start point is inside the plasma, LH or FW case
c        the reading of the data for LH grill
c--------------------------------------------------
         do i=1,ngrill
            powers(i)=powers(i)*1.d-6     !from Watt to MWatt
         enddo 
      endif !istart


c-----transformation data /output/ from genray.in to genray.dat form 
      max_plot_freq=max_plot_freq*1.d-9        !from HZ to GHZ

      if(idens.eq.0) then
         write(*,*)'0 Analytical radial profiles'
c--------transformation data /denprof/ from genray.in to genray.dat form
c        from 1/m**3 to 10**19/m**3
         do i=1,nbulk
            dense0(i)= dense0(i)*1.d-19
            denseb(i)= denseb(i)*1.d-19
         enddo     
      endif ! idens analytical


      if (idens.eq.1) then
c-----------------------------------------------------------------------
c        spline approximation of the density, temperature, zeff,
c        tpop=T_perp/T_parallel, vflow
c        radial profiles
c=====================================================================
c        read density profiles

cSAP080731
c         write(*,*)'nonuniform_profile_mesh=',nonuniform_profile_mesh

c------------------------------------------------------------------ 
         if(nonuniform_profile_mesh.eq.'disabled') then 
c           need not transormation           
         else 
c--------------------------------------------------------------------
c          read density 'profiles dentab_nonuniform_line' at
c          nonuniform mesh in row form
c-----------------------------------------------------------------
c------------transformation data /dentab_nonuniform_line/
c            from genray.in to genray.dat form 
c-----------------------------------------------------------------
c            density in prof2_uniform from 1/m**3 to 10**13/cm**3
c            density in dens1_nonuniform from 1/m**3 to 10**13/cm**3
c-----------------------------------------------------------------
             do i=1,nbulk
               do k=1,ndens                 
                  dens1_nonuniform(k,i)=dens1_nonuniform(k,i)*1.d-19
               enddo
             enddo

         endif !nonuniform_profile_mesh.eq.'disabled'


cSAP080731
c         write(*,*)'before dens1 trasnformation in_rto_dat'
c         write(*,*)'nbulk,ndens',nbulk,ndens

         do i=1,nbulk
	    do k=1,ndens
               dens1(k,i)=dens1(k,i)*1.d-19
cSAP080731
c           write(*,*)'after trans to dat i,k,dens1(k,i)',
c     &                                   i,k,dens1(k,i)

            enddo
         enddo
        
      endif !idens.eq.1

cSAP090209
c-------------------------------------------------------------
c     minimal density outside LCFS from /edge_prof_nml/
c------------------------------------------------------------
      dens_min_edge=dens_min_edge*1.d-19

      return
      end

      subroutine transform_input_data_to_MKSA
c-------------------------------------------------------------------------
c     Converts input data created by subroutine default_in in 
c     genray mixed units to MKSA, in the case that there is a
c     genray.in namelist file in the PWD.  Otherwise, does nothing.
c-------------------------------------------------------------------------
      include 'param.i'
c----------------------------------------
c     declaration and storage of all namelists variables
c-----------------------------------
      include 'cone_nml.i'
c      include 'emissa_nml.i'
      include 'grill_nml.i'
c      include 'ions_nml.i'
      include 'one_nml.i'
c      include 'onetwo_nml.i'
cSAP100921
      include 'output_nml.i'
c      include 'scatnper_nml.i'    
      include 'six_nml.i'  
c      include 'adj_nml.i'
cSAP090209
      include 'edge_prof_nml.i'
c-----local
      integer i_unit,i_genray_in_transformation,kode

c-----Check which input file is given
c     If genray.in  is given then it sets i_genray_in_transformation=1
c     If genray.dat is given then it sets i_genray_in_transformation=0
      i_unit=1
      i_genray_in_transformation=0

c      write(*,*)'subroutine transform_input_data_to_MKSA'

cSAP100914
c      open(i_unit,file = 'genray.in',status='old',iostat=kode)
      open(i_unit,file = 'genray.dat',status='old',iostat=kode)
cSAP100914
c      if (kode.eq.0) then 
      if (kode.ne.0) then
        open(i_unit,file='genray.in',status='old',iostat=kode)  
        if (kode.ne.0) then
            write(*,*)' subroutine transform_default_in_data_to_MKSA'
            write(*,*)' Neither genray.in or genray.dat are presented'
            stop         
         endif
        i_genray_in_transformation=1 !input data file genray.in
                                     !in MKSA system
      endif

cSAP100914
c      if (kode.ne.0) then
c         open(i_unit,file='genray.dat',status='old',iostat=kode)         
c         if (kode.ne.0) then
c            write(*,*)' subroutine transform_default_in_data_to_MKSA'
c            write(*,*)' Neither genray.in or genray.dat are presented'
c            stop         
c         endif
c      endif
      close(i_unit)

      write(*,*)'i_genray_in_transformation',
     &           i_genray_in_transformation

      if(i_genray_in_transformation.eq.1) then 
c--------for namelist /wave/ 
c        transformation data /wave/ from genray.dat to genray.in MKSA form 
         frqncy=frqncy*1.d9 !from GHZ to HZ

c--------for namelist /output/ 
c        transformation data /output/ from genray.dat to genray.in MKSA form 
       
         max_plot_freq=max_plot_freq*1.d9 !from GHZ to HZ

c--------for namelist /denprof/ 
c        transformation data /denprof/ from genray.dat to genray.in MKSA form 
         dense0(1)=dense0(1)*1.d19  !from 10**19 m**(-3)
         dense0(2)=dense0(2)*1.d19  !to          m**(-3)
         denseb(1)=denseb(1)*1.d19
         denseb(2)=denseb(2)*1.d19

c--------for namelist /grill/ 
c        transformation data /grill/ from genray.dat to genray.in MKSA form 
         powers(1)=powers(1)*1.d6 !from MWATT to WATT

c--------for namelist /eccone/ 
c        transformation data /eccone/ from genray.dat to genray.in MKSA form 
         powtot(1)=powtot(1)*1.d6 !from Mwatt to watt

c--------for namelist /dentab/ part
c        transformation density dens1 from genray.dat to genray.in MKSA form
         do i=1,nbulk
           do k=1,ndens 
             dens1(k,i)= dens1(k,i)*1.d19
           enddo
         enddo

c------------------------------------------------------
c        transformation of the density at nonuniform_profile_mesh 
c        from genray.dat to genray.in MKSA form
c-------------------------------------------------------------------
         do k=1,ndensa
            do j=1,nbulka  
               dens1_nonuniform(k,j)=dens1_nonuniform(k,j)*1.d19  !from 10**19/m**3 to 1/m**3
            enddo
         enddo

cSAP090209
c-------------------------------------------------------------------
c       transformation of thec minimal edge density (outside LCFS) 
c       from genray.dat to genray.in MKSA form 
c-------------------------------------------------------------------
        dens_min_edge=dens_min_edge*1.d19 !1/m**3

      endif

      return

      end
c     
c     
      integer function length_char(string)
c     Returns length of string, ignoring trailing blanks.
c     Uses the fortran intrinsic len().
c     The search for a non-blank character is from the end of the
c     string.  Thus, the last non-black character is found, and
c     embedded blanks are ignored.


      character*(*) string
      do i=len(string),1,-1
         if(string(i:i) .ne. ' ') goto 20
      enddo
 20   length_char=i
      return
      end
c     
c     
      integer function length_char1(string)
c     Returns length of string, ignoring characters from
c     the first blank character.
c     Uses the fortran intrinsic len().
      character*(*) string
      do i=1,len(string)
         if(string(i:i) .eq. ' ') goto 20
      enddo
 20   length_char1=i-1
      return
      end

      subroutine read_rz_launch_grill_nml(i_unit,kode)

c-----Reads grill launch points from namelist /rz_launch_grill_tab/
c     given in  genray.in or genray.dat input file
c
c     rlaunch(j,i),zlaunch(j,i) [meters]
c     and weight_power_launch(j,i),[non dimentional]
c     Here
c     i=1,ngrill
c       j=1,nthin(i)
c---------------------------------------------------------
c     Bellow are described input data for rz launch
c     which were set in namelist /grill/.
c
c     ngrill is a number of grills was set in namelist /grill/  
c
c     phigrill(i=1:ngrill) is a toroidal angle of grill, measured counter
c                          clockwise from horizontal through the
c                          magnetic axis (degrees).
c                          They were set in namelist /grill/.
c  
c     nthin(i=1:ngrill)    numbers of poloidal points at each grill
c                          They were set in namelist /grill/.
c  
c     total powers at each grill (i) is  powers(i=1,ngrill) .
c                          It were set in namelist /grill/.
c
c     with weight_power_launch(j,i) non dimensional      
c
c     weight_power_launch(j,i) should be normalized using
c     following formula:  
c     sum{j=1,ntin(j)}(weight_power_launch(j,i))=1.
c
c     The code will check the normalization condition.
c     If the normalization sum.ne.1 then
c     the code will renormalize weight_power_launch(j,i))
c
c     power at each ray witn index (j,i) is powers(i)*weight_power_launch(j,i)
c-----------------------------------------------------------
c     output: kode  < 0   the end of input file was detected
c                   = 0     reading has complited succefully
c                   > 0     an error has occurred
c-------------------------------------------------------------------
      implicit none

      include 'param.i'
      include 'grill_nml.i'

c-----input
      integer i_unit        ! is the number of opened input file 

c-----output
      integer kode          !  < 0   the end of input file was detected
                            !  =0     reading has complited succefully
                            !  >0     an error has occurred

c-----locals
      integer i,j,k
      real*8 sum

      real*8 prof_2d(nthinmax*3,ngrilla)
      namelist /rz_launch_grill_tab/prof_2d

      write(*,*)'in subroutine read_rz_launch_grill_nml'
  
      rewind(unit=i_unit) 

c--------------------------------------------------------------------
c     read namelist /rz_launch_grill_tab/
c     It is used at ilaunch=2
c
c     read r,z coodinates and power weights
c
c     The example for ngrill=2, ntin(1)=3, ntin(2)=4
C--------------------------------------------------------------
!&rz_launch_grill_tab!
!for grill with i=1
!
!      prof_2d(1,1)=
!      rlaunch(j=1,i=1), zlaunch(j=1,i=1), weight_power_launch(j=1,i=1),
!      rlaunch(j=2,i=1), zlaunch(j=2,i=1), weight_power_launch(j=2,i=1),
!      rlaunch(j=3,i=1), zlaunch(j=3,i=1), weight_power_launch(j=3,i=1)
!
! for grill with i=2
!
!      prof_2d(1,2)=
!      rlaunch(j=1,i=2), zlaunch(j=1,i=2), weight_power_launch(j=1,i=2),
!      rlaunch(j=2,i=2), zlaunch(j=2,i=2), weight_power_launch(j=2,i=2),
!      rlaunch(j=3,i=2), zlaunch(j=3,i=2), weight_power_launch(j=3,i=2),
!      rlaunch(j=4,i=2), zlaunch(j=4,i=2), weight_power_launch(j=4,i=2),
!&end
c-------------------------------------------------------------------
      call bcast(prof_2d(1,1),0.d0,ngrilla*nthinmax*3)
      read(i_unit,rz_launch_grill_tab,iostat=kode)
      call check_read(kode,'rz_launch_grill_tab')
      write(*,*)'after rz_launch_grill_tab kode',kode

      write(*,*)'ngrilla,nthinmax',ngrilla,nthinmax

c      do i=1,ngrill
c         write(*,*)'i',i,'nthin(i)',nthin(i)
c         do j=1,nthin(i)*3
c            write(*,*)'j,i,prof_2d(j,i)',j,i,prof_2d(j,i)           
c         enddo
c      enddo

      do i=1,ngrill
         write(*,*)
         write(*,*)'i',i,'nthin(i)',nthin(i)
         write(*,*)
         do j=1,nthin(i)
            rlaunch(j,i)=prof_2d((j-1)*3+1,i)
            zlaunch(j,i)=prof_2d((j-1)*3+2,i)
            weight_power_launch(j,i)=prof_2d((j-1)*3+3,i)

            write(*,*)'j,i,rlaunch(j,i)',j,i,rlaunch(j,i)
            write(*,*)'j,i,zlaunch(j,i)',j,i,zlaunch(j,i)
            write(*,*)'j,i,weight_power_launch(j,i)',
     &                 j,i,weight_power_launch(j,i)
            write(*,*)
         enddo
      enddo

c-----check normalization of Sum{j=1,ntin(i)}weight_power_launch(j,i)=1
c     If weight_power_launch normalization condition was non correct,
c     the code will renormalize weight_power_launch
      do i=1,ngrill

        sum=0.d0
        do j=1,nthin(i)
           sum=sum+ weight_power_launch(j,i)
        enddo

        if (dabs(sum-1.d0).gt.1.d-11) then
           write(*,*)
           write(*,*)'*******************************************'
           write(*,*)'in subroutine read_rz_launch_grill_nml'
           write(*,*)'for grill with number i=',i
          write(*,*)'weight_power_launch normalization was non correct'
           do j=1,nthin(i)
             write(*,*)'i,j,weight_power_launch(j,i)',
     &                  i,j,weight_power_launch(j,i)
           enddo

           write(*,*)'Sum{j=1,ntin(i)}weight_power_launch(j,i)=',sum
           write(*,*)'the code will renormalize weight_power_launch'

           do j=1,nthin(i)
             weight_power_launch(j,i)=weight_power_launch(j,i)/sum
           enddo

           write(*,*)'new renormalized weight_power_launch(j,i)'
           do j=1,nthin(i)
             write(*,*)'i,j,weight_power_launch(j,i)',
     &                  i,j,weight_power_launch(j,i)
           enddo
        endif

      enddo !i
c      stop 'read_rz_launch_grill_nml(i_unit,kode)'

      return
      end

      subroutine write_rz_launch_grill_nml(i_unit,kode)
c-----Writes grill launch points from namelist /rz_launch_grill_tab/
c     given in  genray.in or genray.dat input file
c
c     rlaunch(j,i),zlaunch(j,i) [meters]
c     and weight_power_launch(j,i),[non dimentional]
c     Here
c     i=1,ngrill
c       j=1,nthin(i)
c---------------------------------------------------------
c     Bellow are described input data for rz launch
c     which were set in namelist /grill/.
c
c     ngrill is a number of grills was set in namelist /grill/  
c
c     phigrill(i=1:ngrill) is a toroidal angle of grill, measured counter
c                          clockwise from horizontal through the
c                          magnetic axis (degrees).
c                          They were set in namelist /grill/.
c  
c     nthin(i=1:ngrill)    numbers of poloidal points at each grill
c                          They were set in namelist /grill/.
c  
c     total powers at each grill (i) is  powers(i=1,ngrill) .
c                          It were set in namelist /grill/.   
c     with weight_power_launch(j,i) non dimensional      
c
c     weight_power_launch(j,i) should be normalized using
c     following formula:  
c     sum{j=1,ntin(j)}(weight_power_launch(j,i))=1.
c
c     The code will check the normalization condition.
c     If the normalization sum.ne.1 then
c     the code will renormalize weight_power_launch(j,i))
c
c     power at each ray witn index (j,i) is powers(i)*weight_power_launch(j,i)
c-------------------------------------------------------------------------
      implicit none

      include 'param.i'
      include 'grill_nml.i'

c-----input
      integer i_unit        ! is the number of opened input file 
c-----output
      integer kode          !  < 0   the end of input file was detected
                            !  =0     reading has complited succefully
                            !  >0     an error has occurred


c-----locals
      integer i,j
      real*8 sum

      real*8 prof_2d(nthinmax*3,ngrilla)
      namelist /rz_launch_grill_tab/prof_2d

      write(*,*)'in subroutine write_rz_launch_grill_nml' 
      if(ilaunch.ne.2) return

      do i=1,ngrill
         do j=1,nthin(i)
            prof_2d((j-1)*3+1,i)=rlaunch(j,i)
            prof_2d((j-1)*3+2,i)=zlaunch(j,i)
            prof_2d((j-1)*3+3,i)=weight_power_launch(j,i)
         enddo
      enddo

      write(i_unit,rz_launch_grill_tab,iostat=kode)

      write(*,*)'in write_rz_launch_grill_nml'
      write(*,*)'after write rz_launch_grill_tab'

      return
      end

