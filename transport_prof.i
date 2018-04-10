c This include file originally for passing data from subroutines
c read_transport_prof to write_transport_prof (in partner.f).
c Now include in genray.f in order to set nj and pass it to 
c write_transport_prof for cases in which read_transport_prof
c is not called..
      integer nspeca,nja
      parameter (nspeca=nbulka,nja=nra)     
c               nspeca is a maximum number of transport code species(e+ions)
c               nja is maximum length of transport code small radial mesh 
      real*8 r_transport(nja) !transport code small radial mesh [cm]
      real*8 dr_transp
      integer nj                     !the number of radial points in 
                                     !transport profiles
      common/transport_prof/ nj,r_transport,dr_transp
