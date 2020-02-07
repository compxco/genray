c This include file originally for passing data from subroutines
c read_transport_prof to write_transport_prof (in partner.f).
c Now include in genray.f in order to set nj and pass it to 
c write_transport_prof for cases in which read_transport_prof
c is not called..

      !YuP[2020-01] Collected same-type variables into separate blocks,
      !for better alignment

      integer nspeca,nja
      parameter (nspeca=nbulka,nja=nra)     
c               nspeca is a maximum number of transport code species(e+ions)
c               nja is maximum length of transport code small radial mesh 

      integer nj                     !the number of radial points in 
                                     !transport profiles
      common/transport_prof_int/nj
      
      real*8 r_transport(nja) !transport code small radial mesh [cm]
      real*8 dr_transp
      common/transport_prof/ dr_transp,r_transport