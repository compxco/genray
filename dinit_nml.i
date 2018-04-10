c----------------------------------------
c     declaration of namelist variables which
c     are used in subroutine dinit_mr and
c     are not used in any common block
c-----------------------------------
      integer 
c     ndim1,                        ! /numercl/ it is set in one_nml.i
     &nj_tab(nbulka)                ! number of radial points of nonuniform
                                    ! profiles 

      real*8
     &prof(nbulka*ndensa),          !Working array for uniform and non-uniform
                                    !profiles
     &prof2(nbulka,ndensa),         !Working array for uniform and
                                    !and non-uniform profile          
     &prof_2d(ndensa,nbulka),       !Nonuniform table profile
     &radii_2d(ndensa,nbulka),      !given by lines at radii_2d mesh
     &prof1_uniform(ndensa)         !For uniform grid 1D TSC radial profiles 
c     &zeff1(ndensa)                 !or zeff unform radial profile
c-----------------------------------------------------------------------------
c     ndim1=6 (number of the ray tracing equations)
c     prmt1=tau initial= prmt(1)
c     prmt2=tau final=prmt(2)
c     prmt3=initial tau step=prmt(3)
c     prmt4=required accuracy=prmt(4)
c     prmt6=hprint=prmt(6) time period or poloidal distance [m] 
c                         for results output
c-----------------------------------------------------------------
c     namelists for all table plasma profiles at uniform
c     radial mesh written by columns:
c-----------------------------------------------------------------
c     namelist /dentab/ prof 
c     namelist /temtab/ prof
c     namelist /tpoptab/ prof
c     namelist /vflowtab/ prof
c     namelist /zeftab/ zeff1
c     namelist /EdcTSCtab/ prof
c     namelist /EdcTSCtab/ prof
c-----------------------------------------------------------------
c     namelists for all table plasma profiles at non uniform
c     radial mesh written by columns:
c-----------------------------------------------------------------
c     namelist /dentab_nonuniform/   nj_tab,prof,prof_radii
c     namelist /temtab_nonuniform/   nj_tab,prof,prof_radii
c     namelist /tpoptab_nonuniform/  nj_tab,prof,prof_radii
c     namelist /vflowtab_nonuniform / nj_tab,prof,prof_radii
c     namelist /zeftab_nonuniform /   nj_tab,prof,prof_radii
c-----------------------------------------------------------------
c     namelists for all table plasma profiles at non uniform
c     radial mesh written by lines:
c-----------------------------------------------------------------
c     namelist /dentab_nonuniform_line/   nj_tab,prof_2d,radii_2d
c     namelist /temtab_nonuniform_line/   nj_tab,prof_2d,radii_2d
c     namelist /tpoptab_nonuniform_line/  nj_tab,prof_2d,radii_2d
c     namelist /vflowtab_nonuniform_line / nj_tab,prof_2d,radii_2d
c     namelist /zeftab_nonuniform_line /   nj_tab,prof_2d,radii_2d
c------------------------------------------------------------------
