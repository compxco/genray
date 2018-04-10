c
c
c
      subroutine intinit(omega,n,npar,cdvt,vmax,v0,vpar0,vpar1,
     & vpar2,thet1,thet2,ires)
      implicit double precision (a-h,o-z)

c
c  For integration we require that at least one of the vpar1,vpar2
c  (vperp=0. ) points of the resonance curve lie
c  within +/- vmax.
c  Integration is carried out in the clockwise direction along only
c  that portion of the resonance curve which
c  lies inside vmax and includes that
c  vpar1 or vpar2 which has the smallest absolute value laying
c  inside vmax.
c     Normalized velocities are used.
c-----input omega=emega/omegac=1/y
c           n    is the number of EC harmonic
c           npar is the parallel refructive index
c           cdvt -is the velocity of the normalization (cm/sec)
c           vmax  is max

      double precision npar

c-----output 
c     resonace ellipse N_par**2<1
c     (vpar-vpar0)**2/v0**2+vperp**2/vmax*2=1
c     (v0/c)**2=((n*Y)**2-(1-npar**2))/(1-npar**2)**2
c     vpar0/c=Y*npar/(1=npar**2)    
c     vmax=v0*dsqrt(1-npar**2)
c     vmax is the max V_per at the ellipse
c     v0   is the 
c     vper1=vpar0-dsign(vpar0)*v0
c     vpar2=vpar0+dsign(vpar0)*v0
c
c     ires indicates the subroutine result:
c     ires= 1,  OK.  v0,vpar0,vpar1,vpar2,thet1,thet2 are found.
c           2,  npar**2.ge.1,  no results
c           3,  imaginary resonance curve.
c           4,  integration path not on grid.
c
c     Determine intersections of resonance curve with v=vmax and
c     work out the limits on the angles in the integration
c     theta1,theta 
c  
c---- local
      double precision npar2   !N_par**2
      double precision npar2m1 !N_par**1-1

c 
      pi=4.d0*dtan(1.d0)

      ires=1
      vpar1=0.d0
      vpar2=0.d0
      cdvt2=cdvt*cdvt
      npar2=npar*npar
      npar2m1=npar2-1.d0
c
      if(npar2m1.ge.0.d0)  ires=2 ! hyperbole case
      if(ires.ne.1)  return
c
      v02=(((n/omega)**2+npar2m1)/npar2m1**2)*cdvt*cdvt
      if(v02.le.0d0)  ires=3
c  If no real resonance
      if(ires.ne.1)  return
c
      v0=dsqrt(v02)
      vpar0=-n/omega*npar/npar2m1*cdvt
      s=dsign(1.d0,vpar0)
      vpar1=vpar0-s*v0
      vpar2=vpar0+s*v0
c
c  Discover relationship of resonance ellipse to v=vmax.
      if(dabs(vpar1).ge.vmax)  ires=4
c  If resonance curve at vperp=0. doesn't lie on grid:
      if(ires.ne.1)  return
c
c  Determine intersections of resonance curve with v=vmax and
c  work out the limits on the angles in the integration
c
      thet1=0.d0
      thet2=pi
      if(dabs(vpar2).lt.vmax) return
      gammax=dsqrt(1.d0+vmax*vmax/(cdvt*cdvt))
      vparr=cdvt*(gammax-n/omega)/npar
      vperr=dsqrt(vmax**2-vparr**2)
c
      if(s.gt.0.d0) then
        thet1=0.d0
        thet2=datan2(vperr/dsqrt(-npar2m1),vpar0-vparr)
      else
        thet2=pi
        thet1=datan2(vperr/dsqrt(-npar2m1),vpar0-vparr)
      endif
c
      if(thet1.lt.0.d0) stop 111
      if(thet1.gt.pi) stop 112
      if(thet2.lt.0.d0) stop 113
      if(thet2.gt.pi) stop 114
      if(thet2.le.thet1) stop 114
c
      return
      end
