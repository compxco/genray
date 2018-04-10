
c the data for the spline coefficients for psi(r,z),feqd(psi),pres(psi),
c limiters :z_plus(r) and z_minus(r)
c They are created in equilib.f

c      parameter (nx4a=nxeqda+4,ny4a=nyeqda+4,nrya=ny4a)
c      parameter (nlim4=nlimit+4)
c      parameter nrya=max(nxeqda,nyeqda)+4

      real*8      rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin,
     1            tx,ty,txf,cx,cy,
     4            trlimp,trlimm,
     5            cxlimp,cxlimm,cxy,
     +		  tpres,cpres
      integer     nx,ny,ncx,ncy,ip,im
      common/five/rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin,
     1            tx(nx4a),ty(ny4a),txf(nx4a),cx(nx4a),cy(nrya),
     4            trlimp(nlim4),trlimm(nlim4),
     5            cxlimp(nlim4),cxlimm(nlim4),cxy(nx4a,ny4a),
     +		  tpres(nx4a),cpres(nx4a),
     6            nx,ny,ncx,ncy,ip,im
c
c  nxeqda must be >= nxeqd	number of the mech points along r in eqdsk.file
c  nyeqda must be >= nyeqd	number of the mech points along z in eqdsk.file
c

c   coordinates of Lackner rectangle:
c   rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin

c   coefficients for the spline approximation of
c   psi(r,z): tx(nx4a),ty(ny4a),cxy (cy is a working array)
c   feqd(psi):txf,cx
c   pres(psi):tpres(nx4a),cpres(nx4a)
c   z_limmiter plus(r) or (z_above): trlimp,cxlimp
c   z_limmiter minus(r) or (z_under): trlimm,cxlimm
