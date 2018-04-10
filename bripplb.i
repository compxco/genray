
c These variables are for calculation with the ripple magnetic field.
      common/bripplb/
     1  bripl_r0(nxmaxa,nymaxa),bripl_r1(nxmaxa,nymaxa),
     1  bripl_z0(nxmaxa,nymaxa),bripl_z1(nxmaxa,nymaxa),
     1  bripl_p0(nxmaxa,nymaxa),bripl_p1(nxmaxa,nymaxa),
     1  rippl(nxmaxa,nymaxa),
     1  rippl1(nxmaxa,nymaxa),ripplm(nxmaxa,nymaxa),
     1  ripplt(nxmaxa,nymaxa)
c the arrays for the ripple magnetic field
c b_phi(z,r,phi)=bripl_p0(z,r)+cos(phi)*bripl_p1(z,r)	-toroidal
c b_r(z,r,phi)  =bripl_r0(z,r)+sin(phi)*bripl_r1(z,r)	-r component
c b_z(z,r,phi)  =bripl_z0(z,r)+sin(phi)*bripl_z1(z,r)	-z component
c rippl(z,r)=bripl_p1(z,r)/bripl_p0(z,r)



