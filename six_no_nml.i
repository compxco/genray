
c     Spline coefficients for the density ,temperature,tpop,vflow,zeff
c     radial profiles
c     parameter (ndensa) It is given in param.i
c     parameter (ndens4=ndens+4)
c
c     It has not namelist variables 
c
c     Following arrays have storage in six_nml.i:        
c     dens1(ndensa,nbulka),temp1(ndensa,nbulka),
c     tpop1(ndensa,nbulka),,vflow1(ndensa,nbulka),

c      real*8 
c     & trdens,
c     & rhom,
c     & densm,
c     & cxdens,
c    & trtempe,
c    & trzeff,
c    & cxtempe,
c     & cxzeff,
c     & cxdens1,
c     & cxtemp1,
c     & trdens1,
c     & trtemp1,
c     & trtpop1,
c     & cxtpop1,
c     & trtpop,
c     & cxtpop,
c     & trvflow1,
c     & cxvflow1,
c     & trvflow,
c     & cxvflow,
cSm061205 to use spline functions: coeff1 and terp1
c     & d2_dens_drho1,
c     & d2_temp_drho1, 
c     & d2_tpop_drho1,
c     & d2_vflow_drho1,
c     & d2_zeff_drho1

      real*8, pointer ::
     & cxdens1(:,:),         !(ndens4,nbulk)
     & cxtemp1(:,:),         !(ndens4,nbulk)
     & trdens1(:,:),         !(ndens4,nbulk)
     & trtemp1(:,:),         !(ndens4,nbulk)
     & trvflow1(:,:),        !(ndens4,nbulk)
     & cxvflow1(:,:),        !(ndens4,nbulk)
     & trvflow(:),           !(ndens4)
     & cxvflow(:),           !(ndens4)
     & trtpop1(:,:),         !(ndens4,nbulk)
     & cxtpop1(:,:),         !(ndens4,nbulk)
     & trtpop(:),            !(ndens4)
     & cxtpop(:),            !(ndens4)
     & trdens(:),            !(ndens4)
     & rhom(:),              !(ndens4)
     & densm(:),             !(ndens4)
     & cxdens(:),            !(ndens4)
     & trtempe(:),           !(ndens4)
     & trzeff(:),            !(ndens4)
     & cxtempe(:),           !(ndens4)
     & cxzeff(:),            !(ndens4)
     & d2_dens_drho1(:,:),   !(ndens4,nbulk)
     & d2_temp_drho1(:,:),   !(ndens4,nbulk)
     & d2_tpop_drho1(:,:),   !(ndens4,nbulk)
     & d2_vflow_drho1(:,:),  !(ndens4,nbulk)
     & d2_zeff_drho1(:)      !(ndens4)

      common/six_no_nml/
     & cxdens1,        
     & cxtemp1,        
     & trdens1,
     & trtemp1,
     & trvflow1,
     & cxvflow1,
     & trvflow,
     & cxvflow,
     & trtpop1,
     & cxtpop1,
     & trtpop,
     & cxtpop,
     & trdens,
     & rhom,
     & densm,
     & cxdens,
     & trtempe,
     & trzeff,
     & cxtempe,
     & cxzeff,
cSm061205 to use spline functions: coeff1 and terp1
     & d2_dens_drho1,
     & d2_temp_drho1, 
     & d2_tpop_drho1,
     & d2_vflow_drho1,
     & d2_zeff_drho1

c  ndensa is the max number of points in arrays with the plasma density,
c           temperature and zeff.
c  ndens4=ndens+4 is the number of points in the arrays for the spline
c           approximation of the plasma density, temperature and zeff
c  nbulka  is the max number of plasma component =1 for the electron plasma 
c           must be .ge. nbulk


