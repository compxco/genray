
c     The declarion of namelits data used in common /six_nml/
c
c     zeff is in the namelist:      namelist /zeftab/zeff1
c
c     Arrays dens1,temp1,tpop1,vflow1 are storages for plasma profiles
c
c     parameter (ndensa) It is given in param.i
c     parameter (nbulka) It is given in param.i
       
      real*8 
c-----/zeftab/
     &zeff1,dens1,temp1,tpop1,vflow1, !profiles at uniform radial mesh 
                                      !rhom(i=1,...,ndens)=(i-1)/(ndens-1)
c----calulated using dinit_nml.i data-------------------------------------
     &zeff1_nonuniform,dens1_nonuniform, ! profiles at nonuniform radial mesh 
     &temp1_nonuniform,tpop1_nonuniform, ! radii_...(nbulka,ndensa)
     &vflow1_nonuniform,                 !  
     &radii_nonuniform_zeff1,radii_nonuniform_dens1, !radii_...
     &radii_nonuniform_temp1,radii_nonuniform_tpop1, !
     &radii_nonuniform_vflow1
                        !
      integer nj_tab_zeff1,              !number of radial point
     &nj_tab_dens1,nj_tab_temp1,         !used for the given profile
     &nj_tab_tpop1,nj_tab_vflow1         ! 

      common/six_nml/
     &zeff1(ndensa),
     &dens1(ndensa,nbulka),
     &temp1(ndensa,nbulka),
     &tpop1(ndensa,nbulka),    
     &vflow1(ndensa,nbulka),
     &zeff1_nonuniform(ndensa,nbulka),
     &dens1_nonuniform(ndensa,nbulka),
     &temp1_nonuniform(ndensa,nbulka),
     &tpop1_nonuniform(ndensa,nbulka),
     &vflow1_nonuniform(ndensa,nbulka),               
     &radii_nonuniform_zeff1(ndensa,nbulka),
     &radii_nonuniform_dens1(ndensa,nbulka),
     &radii_nonuniform_temp1(ndensa,nbulka),
     &radii_nonuniform_tpop1(ndensa,nbulka),
     &radii_nonuniform_vflow1(ndensa,nbulka),
     &nj_tab_zeff1(nbulka),
     &nj_tab_dens1(nbulka),
     &nj_tab_temp1(nbulka),
     &nj_tab_tpop1(nbulka),
     &nj_tab_vflow1(nbulka)
c-------------------------------------------------------------------------
c  ndensa is the max number of points in arrays with the plasma density,
c           temperature and zeff.
c  ndens4=ndens+4 is the number of points in the arrays for the spline
c           approximation of the plasma density, temperature and zeff
c  nbulka  is the max number of plasma component =1 for the electron plasma 
c           must be .ge. nbulk
c-----------------------------------------------------------------------------

