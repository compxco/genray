
c     The spline coefficients for the density ,temperature,tpop,vflow,zeff
c     radial profiles
c     parameter (ndensa) It is given in param.i
c     (ndens4=ndens+4)
       
      include 'six_nml.i'
      include 'six_no_nml.i'

c---------------------------------------------------------------------------
c  ndensa is the max number of points in arrays with the plasma density,
c           temperature and zeff.
c  ndens4 is the number of points in the arrays for the spline
c           approximation of the plasma density, temperature and zeff
c  nbulka  is the max number of plasma component =1 for the electron plasma 
c           must be .ge. nbulk
c--------------------------------------------------------------------------

