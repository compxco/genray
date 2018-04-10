
c The arrays for the calculation of the ripple magnetic field using 
c the direct calculations of the ripple field from the coils.
c These data are not used in the present genray version
c      parameter (ntor,nloop,npol)
      common/loopb/
     1          arbphi(ntor,nxeqda,nyeqda),
     1          arbz(ntor,nxeqda,nyeqda),arbr(ntor,nxeqda,nyeqda),
     1          arcosphl(nloop),arsinphl(nloop),arphloop(nloop),
     1          arcosph(ntor),arsinph(ntor),arphi(ntor),
     1          arcospl(npol),arsinpl(npol),artetpol(npol),
     1          arrloop(npol),arzloop(npol),
     1          arcospl_(npol),arsinpl_(npol),artetpl_(npol),
     1          arrloop_(npol),arzloop_(npol),
     1          ardelpol(npol),ardel_z(npol),ardel_r(npol)
c     ntor  is the number of points in the 'toroidal angle' mesh
c     nxeqda is the max number of poins in the 'r-radial coordinate' mesh
c     nyeqda is the max number of poins in the 'z -verticle coordinate' mesh
c     nloop is the number of the loops 
c     npol  is the number of points along the loop 
c     arbhi- toroidal magnetic field
c     arbz - z component of the magnetic field
c     arbr - major radius component of the magnetic field


