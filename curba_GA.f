! YuP[2020-08-20] Renamed all b into b_GA, to avoid conflict with function b()

      subroutine TorGA_curba(rjpd, rjpd0, ratjpd, denom, aspct, enpar, 
     &tc, thtc, theta, elomom, lh, zeff, model, tol, n0, ig)
cProlog

      implicit none
c ----------------------------------------------------------------------
c     CURBA version 7.2  12/11/97, edited by Y.R. Lin-Liu
c     CURBA version 7.1, 11/21/90, R. H. Cohen, LLNL.
c     calculates  J/Pd for ECRH and lower hybrid
c     from numerical solution for Spitzer-Harm
c     problem, including trapped particles, electron-ion collisions,
c     and poloidal variation of collision operator.  J is flux-
c     surface-averaged current density.  Reference:
c     R.H. Cohen, Effect of Trapped Electrons on Current Drive,
c     Phys Fl. 30, 2442 (1987) and erratum in press.
c     Usage notes:
c     J/Pd normalized to e/(m nu vt), where nu =
c     4 pi n e**4 log lambda / m**2 vt**3, and vt = sqrt (2T/m).  Note
c     this normalization differs from that in Reference; J/Pd in
c     reference is factor of 2 larger.
c     When thtc .ne. 1, normalization uses  thermal speed and collision
c      frequency for bulk electrons.
c     When thtc < 0, calculates maximum efficiency for rising buckets.
c     Also provides values of Green's function g (Spitzer-Harm
c       distribution function including trapped-particle effects)
c      --see notes for versions 2.1 and 4.0.  Define angular and radial
c      parts H and F by:
c        g=<rt/R> EXP (-energy/T) c**4 F H/(4 nu vt**3 zhat)
c      where rt is mirror ratio from field minimum,
c      R is major radius, and zhat=(1+charge)/2.
c      For numerical reasons H is normalized differently than in reference;
c      it is larger by zhat.
c     User may wish to alter function FPP which calculates energy-dependent
c      part of Green's function and its derivitive.  See defining comments
c      in routine.
c    OUTPUTS:
c     rjpd is current drive efficiency calculated with trapped particles
c     rjpd0 is non-relativistic current-drive efficiency in limit of no
c     trapped particles, at vparallel=min on relativistic resonance curve.
c     ratjpd is the ratio of the above.
c     denom is normalized absorbed power (denominator in rjpd).
c   INPUTS:
c     enpar is k_parallel*c/wave freq.  Note enpar**2 < 1-elomom**2
c      implies no resonant particles; rjpd set to zero and rjpd0
c      determined for vparallel given by nonrelativistic
c      resonance condition.
c     tc is bulk electron temperature;
c     thtc is ratio of temperature along characteristic to that of bulk,
c      except for thtc < 0, -thtc is energy (in keV) of bucket rise;
c     aspct is inverse aspct ratio;
c     theta is poloidal angle at which power is absorbed (measured
c     from outside); for thtc<0, theta is poloidal angle
c      at which electrons become trapped in bucket; used to
c      calculate resonant energy given elomom and enpar.  Note
c      if calling program fixes minimum resonant energy and
c      calculates enpar, then theta has no significance to the
c      physics; rjpd depends only on Bmin*Y/B(theta), not on
c      theta or Y alone.  But it is useful to be able to specify
c      theta and Y separately in order to keep track of what higher
c      harmonics are doing.
c     ABS (lh) is power of e_perp in diffusion coeff; sign governs
c      power of p_parallel in diffusion coef: + gives p_par^0;
c      - gives p_par^2.  For ECRH, |lh| is
c      harmonic number (lh=1 for fundamental); + for E_-
c      contribution, - for E_parallel; + with lh-->lh+2 for E_+
c      component of electric field E.  For now, there is no provision
c      for the p_par^2 option with lh=0, as the compiler doesn't know
c      the difference betweeen +0 and -0.  if there is any interest,
c      send a message to user 313 and a revision including this option
c      will be created.
c     elomom is ABS (lh)*cyclotron freq over wave freq for Y<0,
c     interpret as evaluated at poloidal angle where electrons are first
c     trapped in bucket.
c     tol is absolute tolerence for DO2GBF integrator; starting in
c      version 1.1 the variables integrated by D02GBE are normalized
c       to be O(1), so tol can be viewed as relative error tolerence.
c     zeff is ion effective charge.
c     model: Absolute value of model selects
c      collisional model: 1 for full bounce av, 2 for square well
c      (numerical solution); 3 for analytic solution to square well.
c      negative model does parallel heating (lower hybrid, fast wave,..).
c     n0 is number of mesh points for gaussian integration of
c      Green's function; n0=2,4,6,8,10,12,16,20,24,32,48 or 64.  n0=128
c      gives 64 points over resonance if resonance has only single
c      passing-particle segment; 64 points on each if two segments.
c     20,24,32,48, or 64.
c     ig selects form for relativistic correction to green's function.
c      ig=1 uses model which agrees with Fisch thru first relativistic
c      correction to nonrelativistic limit and in relativistic limit;
c     ig=2 and 3 are older models; see comments in subroutine fupr.
c     version 1.1, 7/17/86, normalizes integrator variables better.
c     version 2.0, 7/23/86, adds lower hybrid
c     version 2.1 9/15/86, attempted to allow access to the h array on a
c      specified mesh of values of eta = mu*B0/energy.  Here h(1,) is
c      angular part H of generalized Spitzer-Harm function;
c      h(2,)=b*dh(1,)/d eta, where b = <sqrt (1-mirror
c      ratio * eta)>, and h(3,) was running integral of normalized current.
c      This never worked properly; there seems to be no way to force the
c      NAG routine d02gbE to maintain a fixed mesh while maintaining
c      accuracy (except for the case where the number is the maximum number
c      set in the dimension statements).  However, as of version 4.0
c      the interpolated h is directly available to users.
c     version 3.0, 10/12/86, allows different temperature for bulk electrons
c      and along resonance curve.
c     version 4.0, 10/13/86, integrates using gaussian quadrature in
c      EXP (-eperp) to get current, even for model=1,2.  As a consequence,
c      h is automatically available from calling program using h=hfun(eta,
c      hpr), where hpr=d h/d eta (output variable).  Also, h matrix is not
c      recomputed unless needed.
c      To calculate h's only, first call curba from your source program
c      with your values for aspct, zeff, model, and tol (a good
c      value is usually .002); all other input variables are irrelevant,
c      except tc should be nonzero.  then h is available as h=hfun(eta,hpr).
c     version 5.0, 10/15/86, is relativistic; enpar replaces epart
c      as input variable, and tc added.
c     version 5.1, 11/11/86, fixes bug so rjpd=0 when no resonant
c      particles, and adds denom to calling sequence.
c     version 5.2, 11/19/86, provides an improved relativistic correction
c      to Green's function and adds switch ig to select Green's function
c      models.
c     version 5.3, 11/21/86, adds a further-improved relativistic correction
c      to Green's function. ig=1 selects new model; 2 and 3 for older models.
c     version 5.4, 12/2/86, improves integration along resonance,
c      eliminating trapped-particle region from domain.
c     version 5.5, 2/16/87, eliminates bug in upper limit for power
c      integration
c    Version 6.0, 4/2/87, adds bucket-lift current-drive calculation,
c     for particles trapped at poloidal angle theta with specified elomom,
c     at zero perpendicular energy, and ending on the same characteristic
c     with energy increased by -thtc (keV).
c    Version 6.1, 12/6/87, corrects comments on version 2.1 and adds
c     capability of using model=3 (analytic Legendre function form
c     of Green's function) in accessing h via hfun from calling program.
c     Also removes possible crash in legfn1 when z**2=1 by substituting
c     limiting form for derivative of Legendre function.
c     Version 7.0, 3/7/88, adds capability of higher-order flr to
c      parallel heating calculation.  lh must now be set; diffusion
c      coefficient varies as eperp**|lh|.  Sign of lh determines
c      whether to weight diffusion coefficient with p_parallel^2;
c      negative lh adds the weighting.  Note that physically we
c      expect weighting with v_parallel^2, which is constant on
c      a Landau type resonance and hence irrelevant.  Hence one
c      will almost certainly want to set lh .ge. 0.
c     Version 7.1, 11/21/90, corrects bug in ffp, for ig=1.  In
c      expression for fupr, earlier version had -.5/bfac instead
c      of -.5*agrel/afac.  Also puts in sign of bucket lift current,
c      not tracked previously.
c------------------------------------------------------------------------
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8 
     & aspect, sgnj, enpar, charge, zeff, theta, rabs, 
     &rmax, etam, etamh, et0, thtcu, thtc, theth, tc, aspct, costh, rmin
     &, elom, elomom, expeu, enpar2, epst1, epst2, etmax, rjpd, rjpd0, 
     &ratjps, denom, expel, expe1, expe2, rint, qgauleg_GA, expeld,zhat, 
     &erise, etal, sgnj0, gammin, etau, dfh, h3fac, sigma, eparnr,
     &fjpd0_GA
     &, fjpd0h_GA, tol, epar, rl, asp1, ap0, alfsq, elooro, r3, eparsq, 
     &emin, gammax, ratjpd, ra0, elom2, agrel, bgrel, qgrel, tor
      integer kg, ig, model, iheat, modela, l, modl2, ifail, n, n0, ixo,
     & lh, lm1, lp1
c
      common /TorGA_fcncom/ aspect,epar,rl,rmin,rmax,rabs,etam,charge, 
     &asp1,zhat,ap0,alfsq,elooro,r3, eparsq,l,iheat,lm1,lp1,modl2
      common /TorGA_rela/ enpar2,elom,theth,emin,ra0,elom2,tor, gammin,
     &agrel,bgrel,qgrel,ixo,kg
c
cSAP080617
      real*8 fcn_GA,fcn

      external fcn_GA
      external fcd
c
c       

cSAP0890617
c      write(*,*)'TorGA_curba aspct,enpar',
c     & aspct,enpar
c      write(*,*)'tc, thtc, theta, elomom',tc, thtc, theta, elomom
c      write(*,*)'lh,zeff,model,tol,n0,ig',lh,zeff,model,tol,n0,ig

      kg = ig-2
      aspect = aspct
      sgnj = 1.0d0
      if (enpar .lt. 0) sgnj = -1.0d0
      charge = zeff
      costh = COS (theta)
      rabs = 1.0d0+aspect*costh
      rmax = 1.0d0+aspect
      rmin = 1.0d0-aspect
      etam = rmin/rmax
      etamh = etam
      et0 = 0.0d0
      thtcu = thtc

cSAP080617
c      write(*,*)'kg,aspect,sgnj,charge',kg,aspect,sgnj,charge
c      write(*,*)'costh,rabs,rmax,rmin',costh,rabs,rmax,rmin
c      write(*,*)'etam,etamh,et0,thtcu',etam,etamh,et0,thtcu

c
c     use th = tc if doing buckets
c
      if (thtc .le. 0.0d0) thtcu = 1.0d0
      theth = tc*thtcu/511.0d0

cSAP080617
c      write(*,*)' theth', theth

      if (model .lt. 0) then
      iheat = -1
      elom = 0
      modela = -3
      else
      iheat = 1
      elom = elomom
      modela = 3
      endif
      ixo = 0
      l = IABS (lh)
      !write(*,*)' TorGA_curba: lh=',lh
c
c     set flag to get ppar**(2*ixo) factor in diffusion coef.;
c     for ecrh, ixo=0 for E-, E+ and =1 for E_par contributions
c
      if (lh .lt. 0) ixo = 1
      modl2 = modela-2
      call TorGA_setpar
c
c     do current-drive integral:
c
      ifail = 0
c
c --- set limits of integration for expe = EXP (-epst), where epst=e-emin,
c     normalized to th.
c
      expeu = 1.0d0
      enpar2 = enpar*enpar
      call TorGA_limits (epst1, epst2, etmax)

cSAP080617
!      write(*,*)'after  TorGA_limits epst1, epst2, etmax',
!     & epst1, epst2, etmax
       
      if (etmax .le. -1.0d0 .or. ABS(enpar) .lt. 1.d-6) then
      rjpd=0.d0
      rjpd0=0.d0
      ratjps=1.d0
      denom=1.d0
      return
      endif
c
c      write(*,*)'etmax,thtc', etmax,thtc

      if (etmax .gt. 0.d0 .and. thtc .gt. 0.d0) then
      if (epst1 .gt. 0.d0 .and. epst2 .lt. etmax) then
c
c     ---calculate integral for resonance divided by trapped region:
c
        expel = EXP (-etmax)
        expe1 = EXP (-epst1)
        expe2 = EXP (-epst2)
        n = n0 / 2
        rint = qgauleg_GA(fcn_GA,expel,expe2,n,ifail)

!        write(*,*)'1,expel,expe2,n,ifail,rint',expel,expe2,n,ifail,rint

        rint = rint+qgauleg_GA(fcn_GA,expe1,expeu,n,ifail)

!        write(*,*)'2 expe1,expeu,n,ifail,rint', expe1,expeu,n,ifail,rint

      else
c
c     ---calculate integral for resonance which doesn't enter trapped region
c
        n = MIN0 (64, n0)
        if (epst2 .le. 0.0d0) expel = EXP (-etmax)
        if (epst2 .gt. 0.0d0 .and. epst1 .le. 0.0d0) expel = EXP (-MIN (
     &  epst2, etmax))
        if (epst2 .gt. 0.0d0 .and. epst1 .gt. 0.0d0) expel = EXP (-MIN (
     &  epst1, etmax))
        rint = qgauleg_GA(fcn_GA,expel,expeu,n,ifail)

!        write(*,*)'3 expel,expeu,n,ifail,rint', expel,expeu,n,ifail,rint

      endif ! (epst1 .gt. 0.d0 .and. epst2 .lt. etmax)

cSAP080617
c      write(*,*)'rint',rint

c
c     ---calculate denominator
c
      expeld = EXP (-etmax)
c%LL         call denomf(denom,expeld,etmax)
      !write(*,*)'expeld,expeu=',expeld,expeu
      denom=theth*qgauleg_GA(fcd,expeld,expeu,n,ifail)
      !write(*,*)'denom,zhat,ifail', denom,zhat,ifail

c%LL         PRINT *, denom, zdenom
      rjpd = -0.125d0*sgnj*thtcu*rmax*rint/(denom*zhat)

c      write(*,*)'1 rjpd', rjpd 

c
c bucket lift calculation:
      else if (etmax .gt. 0.d0 .and. thtc .le. 0) then 
      erise = -thtc/511.0d0
      etal = 0.0d0
      sgnj0 = SIGN (1.0d0,(gammin-elom)*enpar)
      gammax = gammin+erise
      etau = MIN (1.0d0,2.0d0*(rabs/rmax)*elom*erise/ (gammax*gammax-1.0
     &d0))
      call TorGA_getdfh (dfh, gammin, gammax, etal, etau)
      rjpd = -(511.0d0/(8.d0*zhat*tc))*sgnj0*rmax*dfh/erise

c      write(*,*)'2 rjpd', rjpd 

c
c set j/pd to zero if no resonance
      else 
      rjpd = 0.0d0
      gammin = 1.0d0
      endif ! (etmax .gt. 0.d0 .and. thtc .gt. 0.d0)
      
      
      h3fac = 4.0d0 * zhat / (5.0d0 + charge)
      rjpd = rjpd*h3fac

c      write(*,*)'3 rjpd', rjpd 

c
c --- calculate non-relativistic, zero-trapped-particle limits
c
      sigma = 1.0d0 - elomom
      eparnr = 0.5d0*(1.0d0-elom/gammin)**2/(theth*enpar2)
      sgnj0 = SIGN (1.0d0,(gammin-elom)*enpar)
      if (iheat .eq. +1) rjpd0 = sgnj0*fjpd0_GA(eparnr,l,sigma,zeff)
     &                          *thtcu
      if (iheat .eq. -1) rjpd0 = sgnj0*fjpd0h_GA(eparnr,zeff)*thtcu
      ratjpd = rjpd / rjpd0

c      write(*,*)' ratjpd,rjpd,rjpd0', ratjpd,rjpd,rjpd0

      return
c
      end subroutine TorGA_curba





      subroutine TorGA_setpar
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8 
     & rl, ra0, rmax, rabs, r3, asp1, aspect, zhat, 
     &charge, agrel, bgrel, ap0, rmin, eparsq, epar, alfsq, elooro, elom
     &, elom2, rlams, etam, alffac, pas, tor, enpar2, theth, emin, 
     &gammin, qgrel
      integer l, lm1, lp1, modl2, nc, nf, iheat, ixo, kg
c
c     set eta-independent parameters
ccSAP080617
      complex*16 
     &alpha,chalf,cwun,ci,zarg,cpas,q
      common /TorGA_fcncom/ aspect,epar,rl,rmin,rmax,rabs,etam,charge, 
     &asp1,zhat,ap0,alfsq,elooro,r3, eparsq,l,iheat,lm1,lp1,modl2
      common /TorGA_rela/ enpar2,elom,theth,emin,ra0,elom2,tor, gammin,
     &agrel,bgrel,qgrel,ixo,kg
      common /TorGA_mdl3com/alpha,rlams,pas
      data chalf/(-.5d0,0.0d0)/,cwun/(1.0d0,0.0d0)/,ci/(0.0d0,1.0d0)/
c
      rl=l
      ra0=rmax/rabs
      r3=3.0d0*ra0
      asp1=1.0d0-aspect
      zhat=(1.0d0+charge)/2.d0
      agrel=(1.0d0/(5.d0+charge))**2
      bgrel=1.0d0-agrel
      qgrel=1.0d0/bgrel
      ap0 = -0.5d0*sqrt (rmax/rmin)*rmax
      eparsq=epar*epar
      alfsq=-2.0d0*aspect/asp1
      lm1=l-1
      lp1=l+1
      elooro=elom/ra0
      elom2=elom*elom
c     parameters for model=3
      rlams = sqrt (1.0d0-etam)
      alffac=(4.d0/zhat)-.25d0
      if (alffac) 10,10,20
   10 alpha=chalf+cwun*sqrt (-alffac)
      go to 30
   20 alpha=chalf+ci*sqrt (alffac)
   30 zarg = dcmplx(rlams, 0.0d0)

cSAP080617
c      write(*,*)'TorGA_setpar modl2,nf ',modl2,nf 

cSAP080618
      nf=0
      if (modl2 .eq. 1) call TorGA_legfn(alpha,zarg,cpas,q,nc,nf)

cSAP080617
c      write(*,*)'after TorGA_legfn alpha,zarg,cpas,q,nc,nf',
c     &alpha,zarg,cpas,q,nc,nf

      pas = dble(cpas)

cSAP080617
c      write(*,*)'pas',pas

      tor=2.0d0/ra0
      return
c
      end





      subroutine TorGA_limits (epst1,epst2,etmax)
cProlog

      implicit none
c
c     calculate range of relativistic range of integration
c     etmax is emax-emin, normalized to thot.
c
c explicit type declaration 7/12/01 (RAJ)

cSAP080617
      real*8 
     &etmax0, gam1, gam2, gammin, etmax, etcrit, theth, 
     &gam3, gam4, etam, epst1, epst2, emin, enpar2, aspect, epar, rl, 
     &rmin, rmax, rabs, charge, asp1, zhat, ap0, alfsq, qgrel, r3, 
     &eparsq, elom, ra0, elom2, tor, agrel, bgrel, elooro
      integer iheat, l, modl2, lm1, lp1, ixo, kg
c
      common /TorGA_fcncom/ aspect,epar,rl,rmin,rmax,rabs,etam,charge, 
     &asp1,zhat,ap0,alfsq,elooro,r3, eparsq,l,iheat,lm1,lp1,modl2
      common /TorGA_rela/ enpar2,elom,theth,emin,ra0,elom2,tor, gammin,
     &agrel,bgrel,qgrel,ixo,kg
      data etmax0/20.d0/
c
      if (iheat) 60,60,5
c     ECH limits:
c     find intersections of resonance with ppar axis:
    5 call TorGA_gammas(gam1,gam2,0.0d0)
      if (gam1-1.0d0) 10,10,20
   10 if (gam2 .le. 1.0d0) go to 40
      gammin=gam2
      etmax=etmax0
      go to 30
   20 gammin=gam1
      etcrit=(gam2-gam1)/theth
      etmax=MIN (etmax0,etcrit)
c     find intersections of resonance with passing-trapped separatrix:
   30 call TorGA_gammas(gam3,gam4,etam)
      epst1=(gam3-gammin)/theth
      epst2=(gam4-gammin)/theth
      emin=(gammin-1.0d0)/theth
      return
c
c --- no resonance if elom+enpar2-1 < 0; set etmax negative as a flag.
c
   40 etmax=-1.0d0
      return
c     Limits for Landau resonance:
   60 etmax=etmax0
      epst2=-1.0d0
      gammin=1.0d0/sqrt (1.0d0-1.0d0/enpar2)
      emin=(gammin-1.0d0)/theth
      return
c
      end
      subroutine TorGA_getdfh (dfh, gammin, gammax, etal, etau)
cProlog

      implicit none
c
c     calculates difference of final and initial greens functions for
c     bucket lift from etal,gammin to etau,gammax
c
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8 
     & gammin, fl, fupr, gammax, fu, hl, hpr, etal, etau,
     & etam, hu, dfh, gam, eta, psq, gamsq, aspect, epar, rl, rmin, rmax
     &, rabs, charge, asp1, zhat, ap0, alfsq, elooro,r3, eparsq
      integer l, modl2, lm1, lp1, iheat
c
      common /TorGA_energy/gam,eta,psq,gamsq
      common /TorGA_fcncom/ aspect,epar,rl,rmin,rmax,rabs,etam,charge, 
     &asp1,zhat,ap0,alfsq,elooro,r3, eparsq,l,iheat,lm1,lp1,modl2
c
      call TorGA_setrel(gammin)
      call TorGA_ffp(fl,fupr)
      call TorGA_setrel(gammax)
      call TorGA_ffp(fu,fupr)
c      if (modl2) 10,10,20
c   10 hl=hfun(etal,hpr)
c      if (etau .ge. etam)  go to 30
c      hu=hfun(etau,hpr)
c      dfh=fu*hu-fl*hl
c      return
   20 call TorGA_geth3(hl,hpr,etal)
      if (etau .ge. etam) go to 30
      call TorGA_geth3(hu,hpr,etau)
      dfh=fu*hu-fl*hl
      return
   30 dfh=-fl*hl
      return
c
      end




      subroutine TorGA_legfn (v, z, p, q, nc, nf)
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617

      real*8 
     &  ac, acc, pi, sv, cv, eim, shv, chv, vr, vi, r1, r2
     &, r3, r4, rr, vri, eip, pisr
cSAP080617
      real*8 
     & z_imag, zn24
      integer nvv, ninter, n23, nff, nf, n24, nfrig, ncvg, nc
c
cSAP080617
      real*8 
     & qinf
cSAP080617
      complex*16
     & vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv,zz1,zz2,srz,
     &v,z,p,q, zzs,cisp,cism,pt,vvp,zzz, csqrtk_GA
c....&|--1---------2---------3---------4---------5---------6---------7-|
      common /TorGA_legbl/ vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv, zz1,
     &zz2,srz,acc,pisr,pi,r1,r2,r3,r4,ncvg,nfrig
c
      data ac /0.0000001d0/
      data qinf /1.79769313486231d+308/
c     uses c311bd, legv, leg1, legz, leror
c

cSAP080629: nf  initialization
      nf=0

cSDAP080617
c      write(*,*)'TorGA_legfn v, z, p, q, nc, nf',v, z, p, q, nc, nf
c      write(*,*)'pp',pp

      call TorGA_c311bd

c      write(*,*)'v,ac',v,ac

      acc = ac**2
      if (dble(v)+0.5d0) 16,17,17
   16 vv = -v - 1.0d0
      go to 18
   17 vv = v
   18 nvv = ninter (dble(vv))

cSAP080617
c      write(*,*)'18 vv,nvv',vv,nvv

   30 vvp = vv * pi

cSAP080617
c      write(*,*)'30 vvp',vvp

      go to (21, 22, 23,2 4, 25), nvv
   21 sv = 1.0d0
      go to 26
   23 sv = -1.0d0
   26 cv = 0.0d0
      go to 28
   22 cv = -1.0d0
      go to 27
   24 cv = 1.0d0
   27 sv = 0.0d0
      go to 28
   25 sv = SIN (dble(vvp))
      cv = COS (dble(vvp))
   28 eip = EXP ( dimag(vvp))
      eim = EXP (-dimag(vvp))
      shv = 0.5d0 * (eip-eim)
      chv = 0.5d0 * (eip+eim)
      svv = sv*chv+u*cv*shv
      cvv = cv*chv-u*sv*shv
      cisp = cvv+u*svv
      cism = cvv-u*svv

cSAP080617
c      write(*,*)'28 sv,cv',sv,cv
c      write(*,*)'eip,eim,vvp',eip,eim,vvp
c      write(*,*)'shv,chv',shv,chv
c      write(*,*)'svv,cvv,u',svv,cvv,u
c      write(*,*)'cisp,cism,z',cisp,cism,z
c      write(*,*)'pp',pp

      if (dble(z)) 9, 10, 11
    9 zz = -z
      n23 = 3
      nff = -nf

c      write(*,*)'9 zz,nff',zz,nff

      go to 12
   10 n24 = 4
      zz = z
      if (dimag(zz)) 13, 77, 13
   11 zz = z
      n23 = 2
      nff = nf
   12 n24 = 2
   13 zzs = zz**2 
      z1 = (1.0d0-zz) / 2.0d0

c      write(*,*)'13 zz,zzs,z1',zz,zzs,z1
c      write(*,*)'n23,nf',n23,nf

      if (dimag(zz)) 7, 6, 7
    6 if ( dble(z1)) 4, 70, 5
    4 nfrig = ISIGN (n23, nff)

c      write(*,*)'4 n23,nff, nfrig', n23,nff, nfrig

      go to 8
    5 nfrig = nff
      go to 8
    7 continue
      z_imag = dimag(zz)

c      write(*,*)'7 zz,z_imag',zz,z_imag

      zn24 = FLOAT (n24)

c      write(*,*)'n24,zn24',n24,zn24

      nfrig = SIGN (zn24, z_imag)

c      write(*,*)'nfrig',nfrig

    8 z2 = 1.0d0/zzs
      srz = csqrtk_GA (zzs-1.0d0, nfrig, 1)

c      write(*,*)'8 zzs,nfrig,srz',zzs,nfrig,srz

      zz1 = ( zz+srz) / (2.0d0*srz)
      zz2 = (-zz+srz) / (2.0d0*srz)

c      write(*,*)'zz,srz,zz1,zz2',zz,srz,zz1,zz2
c      write(*,*)' dble(z)',dble(z)

      if (dble(z)) 1, 2, 2
    1 srz = -srz
    2 zzz = zz
      zz = z
      vr = dble(vv)**2
      vi = dimag(vv)**2
      r1 = ABS (z1)
      r2 = ABS (z2)
      r3 = ABS (zz1)
      r4 = ABS (zz2)
c
      if (r1-r2)61,61,62
   61 rr = MAX (r3,r4)/0.045d0-19.5d0
      if (rr)66,66,65
   65 vri = vr+vi
      if (vri)165,67,165

c      write(*,*)'before 165 pp',pp

  165 if (rr**2-vri+vi**2/(2.0d0*vri)) 66, 67, 67
   66 call TorGA_legv

c      write(*,*)'66 after TorGA_legv pp',pp

      if (ncvg)166,97,166
  166 zz = zzz
      call TorGA_leg1

c      write(*,*)'166 after TorGA_leg1 pp',pp

      go to 90
   67 zz = zzz
      call TorGA_leg1

c      write(*,*)'67 after TorGA_leg1 pp',pp

      if (ncvg)167,90,167
  167 zz = z
      call TorGA_legv

c      write(*,*)'167 after TorGA_legv pp',pp

      go to 97
   62 if (vr+vi-16.0d0*r2**2)63,64,64
   64 call TorGA_legv

c      write(*,*)'64 after TorGA_legv pp',pp

      if (ncvg)164,97,164
  164 call TorGA_legz

c      write(*,*)'164 after TorGA_legz pp',pp

      go to 97
   63 call TorGA_legz

c      write(*,*)'63 after TorGA_legz pp',pp

      if (ncvg)163,97,163
  163 call TorGA_legv

c      write(*,*)'163 after TorGA_legv pp',pp

      go to 97
   90 if (dble(z)) 93,97,97
   93 pt = -2.0d0/pi*svv*qq
      if (nfrig) 94,95,96
   94 pp = cisp*pp+pt
 
c     write(*,*)'94 pp,cisp,pt',pp,cisp,pt

      qq = -cism*qq
      go to 97
   95 qq = -cvv*qq-pi/2.0d0*svv*pp
      pp = cvv*pp+pt

c      write(*,*)'95 pp,cvv,pt', pp,cvv,pt

      go to 97
   96 pp = cism*pp+pt

c      write(*,*)'96 pp,cism,pt',pp,cism,pt

      qq = -cisp*qq
      go to 97
c
   97 if (dble(v)+0.5d0) 91,92,92
   91 if (ABS (svv) .ne. 0.0d0) go to 98
      qq = qinf
      go to 92
   98 qq = (qq*svv-pi*cvv*pp)/svv
   92 p = pp

c      write(*,*)'92 p,pp',p,pp

      q = qq
      nc = ncvg
      return
   70 qq = qinf
c
      if ( dble(z)) 173, 173, 74
  173 if (dimag(v)) 71, 73, 71
   73 go to (71,72,71,74,71),nvv
   71 pp = qq

c      write(*,*)'71 pp,qq',pp,qq

      go to 82
   72 pp = -1.0d0
      go to 82
   74 pp = 1.0d0
      go to 82
c
   77 nfrig = nf
      call TorGA_legor

c      write(*,*)'77 after TorGA_legor pp',pp

   82 ncvg = 0
      go to 92
c
      end



      subroutine TorGA_setrel (gamma)
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8 
     & gam, gamma, gamsq, psq, eta
c
      common /TorGA_energy/ gam, eta, psq, gamsq
c
      gam = gamma
      gamsq = gam * gam
      psq = gamsq - 1.0d0
      return
c
      end

      subroutine TorGA_ffp (fu, fupr)
cProlog

      implicit none
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8 
     & cgrel, bfac, bgrel, psq, bfacq, qgrel, afac, agrel
     &, rutaf, fu, fupr, gam, gamsq, denom, eta, enpar2, elom, theth, 
     &emin, ra0, elom2, tor, gammin
      integer kg, ixo
c
c     energy-dependent part of Green's function and derivative.
c     the Green's function is defined as:
c       g=<B_t/R>c^4/4/nuv_t^3B_0)Fu(v/c)h(eta)
c     for kg+2=ig=1, fu=p^2(1-1/(1+bgrel*p^2)**qgrel)/sqrt (1+agrel*p^2)
c      with agrel=1/(5+Z)^2, bgrel=1-agrel, qgrel=1/bgrel.  This form
c      gives Green's function which is correct in nonrelativistic
c      limit; in limit of no trapped particles reproduces leading-order
c      (v^6/c^6) relativistic correction as well as ultrarelativistic
c      limit as given by Fisch.
c     for kg+2=ig=2, fu=gamma*v^4/c^4, (correct nonrelativistically; correct
c      energy dependence at large energy but wrong coefficient)
c     for kg+2=ig=3, fu=p^4/(1+cgrel*p^2)^(2/3) with
c      p=momentum/mc and cgrel=2/3 (correct to order v^6/c^6 and correct
c      energy dependence at large energy but wrong coefficient)
c
      common /TorGA_energy/gam,eta,psq,gamsq
      common /TorGA_rela/ enpar2,elom,theth,emin,ra0,elom2,tor, gammin,
     &agrel,bgrel,qgrel,ixo,kg
      data cgrel/.66666666666667d0/
c
      if (kg) 10,20,30
c
c --- ig=1 ------------------
c
   10 bfac=1.0d0+bgrel*psq
      bfacq=bfac**qgrel
      afac=1.0d0+agrel*psq
      rutaf = sqrt (afac)
      fu=psq*(1.0d0-1.0d0/bfacq)/rutaf
      fupr=2.0d0*gam*(fu*(1.0d0/psq-.5d0*agrel/afac)+psq/(rutaf*bfac*
     &bfacq))
      return
c
c --- ig=2 ------------------
c
   20 fu=psq*psq/(gamsq*gam)
      fupr=psq*(1.0d0+3.0d0/gamsq)/gamsq
      return
c
c --- ig=3 ------------------
c
   30 denom=1.0d0+cgrel*psq
      fu=psq*psq/denom**1.5d0
      fupr=fu*gam*(4.d0+cgrel*psq)/(psq*denom)
      return
c
      end

      subroutine TorGA_geth3 (hh, hpr, eta)
cProlog

      implicit none
c
c     gets angle part of green's function in Legendre-function (square-well)
c     approx.  Note that h calculated here needs to be multiplied by
c     4/(Z+5) to get h of UCRL-95813
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8 
     & rlam, eta, rpa, rppr, hh, rlams, pas, hpr
      integer nc, nf
c
cSAP080617
      complex*16 
     &alpha, zarg, pa, ppr, q
      common /TorGA_mdl3com/ alpha, rlams, pas
c
c      write(*,*)'TorGA_geth3 eta',eta

      rlam = sqrt (1.0d0-eta)
      zarg = dcmplx(rlam,0.0d0)

cSAP080617
c      write(*,*)' rlam,zarg',rlam,zarg

      call TorGA_legfn1(alpha,zarg,pa,ppr,q,nc,nf)

c      write(*,*)'after TorGA_legfn1 alpha,zarg,pa,ppr,q,nc,nf',
c     & alpha,zarg,pa,ppr,q,nc,nf

      rpa = dble(pa)
      rppr = dble(ppr)

c      write(*,*)'pa,rpa',pa,rpa
c      write(*,*)'ppr,rppr',ppr,rppr

      hh=-rlam+rlams*rpa/pas

c      write(*,*)'hh,rlam,rlams,rpa,pas',hh,rlam,rlams,rpa,pas

      hpr=(1.0d0-rlams*rppr/pas)*.5d0/rlam

c      write(*,*)'hpr,rlams,rppr,pas,rlam',hpr,rlams,rppr,pas,rlam

      return
c
      end



      subroutine TorGA_c311bd
cProlog

      implicit none
c
cSAP080617
      complex*16 
     &vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv,zz1,zz2,srz
c....&|--1---------2---------3---------4---------5---------6---------7-|
      common /TorGA_legbl/ vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv, zz1,
     &zz2,srz,acc,pisr,pi,r1,r2,r3,r4,ncvg,nfrig
c
c explicit type declaration added 7/12/01 RAJ
cSAP080617
      real*8 
     & pisr, pi, acc,r1, r2, r3, r4
      integer ibd, ncvg, nfrig
c
      data ibd/0/
c
      if (ibd .ne. 0) return
      ibd = 1
      pisr = 1.77245385090552d0
      pi = 3.141592653589793d0
      u = (0.0d0,1.0d0)
      return
c
      end



      subroutine TorGA_legv
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8 
     &acc, sgn, pisr, pi, r1, r2, r3, r4
      integer ncvg, nfrig, ncv
cSAP080617
      complex*16 
     & vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv,zz1,zz2,srz
c
c....&|--1---------2---------3---------4---------5---------6---------7-|
      common /TorGA_legbl/ vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv, zz1,
     &zz2,srz,acc,pisr,pi,r1,r2,r3,r4,ncvg,nfrig
c
cSAP080617
      complex*16 
     & f1, f2, ssrz, clogok_GA, csqrtk_GA, rgam_GA
c     uses hypgm, clogok_GA, csqrtk_GA, rgam_GA
c
      a=0.5d0
      c=vv+1.5d0

c      write(*,*)'in TorGA_legv a,c,zz1,zz2',a,c,zz1,zz2

      call TorGA_hypgm(a,a,c,zz1,f1,acc,ncvg)
cyup      write(*,*)'f1,ncvg,acc',f1,ncvg,acc
      
      call TorGA_hypgm(a,a,c,zz2,f2,acc,ncv )

c      write(*,*)'f2,ncv,acc',f2,ncv,acc
      
      ncvg=ncvg+2*ncv

    6 f1=f1*EXP ((vv+0.5d0)*clogok_GA(zz +srz,nfrig,3))

c      write(*,*)'zz,srz,nfrig,clogok_GA(zz +srz,nfrig,3)',
c     &zz,srz,nfrig,clogok_GA(zz +srz,nfrig,3)
c      write(*,*)'f1,zz,srz,nfrig',f1,zz,srz,nfrig

      f2=f2*EXP ((vv+0.5d0)*clogok_GA(zz -srz,-nfrig,3))

c      write(*,*)'zz,srz,nfrig,clogok_GA(zz -srz,-nfrig,3)',
c     &zz,srz,nfrig,clogok_GA(zz -srz,-nfrig,3)

c      write(*,*)'f2,zz,srz,nfrig',f2,zz,srz,nfrig

      a=1.5d0
      b_GA=1.0d0
      ssrz = csqrtk_GA (2.0d0*srz, nfrig, 2) * rgam_GA (vv, a, b_GA)
 
c      write(*,*)'srz, nfrig,csqrtk_GA (2.0d0*srz, nfrig, 2)',
c     &srz, nfrig,csqrtk_GA (2.0d0*srz, nfrig, 2)
c       write(*,*)'vv,a,b_GA,rgam_GA (vv, a, b_GA)',
c     &vv,a,b_GA,rgam_GA (vv, a, b_GA)
         
c      write(*,*)'ssrz,vv',ssrz,vv

      sgn=1.0d0
      if (dimag(zz)) 8, 14, 12
    8 sgn=-1.0d0
      go to 12
   14 sgn=SIGN (1.0d0,(FLOAT (nfrig)+0.5d0) * dble(zz))
   12 pp=(f1+sgn*u*f2)/(pisr*ssrz)

c      write(*,*)'12 pp,f1,sgn,u,f2,pisr,ssrz',pp,f1,sgn,u,f2,pisr,ssrz

      if (nfrig)11,10,11
   10 qq=0.5d0*pisr*(f2+sgn*u*f1)/ssrz
      go to 80
   11 qq=pisr*f2/ssrz
   80 return
c
      end




      subroutine TorGA_leg1
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8 
     & gamma, acc, el, crit_GA, pisr, pi, r1, r2, r3, r4
      integer ncvg, nfrig, l
cSAP080617
      complex*16 
     & vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv,zz1,zz2,srz
c
c....&|--1---------2---------3---------4---------5---------6---------7-|
      common /TorGA_legbl/ vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv, zz1,
     &zz2,srz,acc,pisr,pi,r1,r2,r3,r4,ncvg,nfrig
c
cSAP080617
      complex*16 
     & sigma, fac
     &, clogok_GA, psifun_GA
      data gamma /0.5772156649d0/
c     uses hypgm, clogok_GA, psifun_GA, crit_GA
c
      a=-vv
      b_GA=vv+1.0d0
      c=(1.0d0,0.0d0)
      call TorGA_hypgm(a,b_GA,c,z1,pp,acc,ncvg)
      sigma=(0.0d0,0.0d0)
      fac=(1.0d0,0.0d0)
      a = 0.5d0 * clogok_GA((zz+1.0d0)/(zz-1.0d0),nfrig,1) - gamma -
     & psifun_GA(b_GA)
      qq=a
      do 17 l=1,50
      el=l
      sigma=sigma+1.0d0/el
      fac=-fac*(vv+el)*(vv-el+1.0d0)*z1/(el*el)
      b_GA=(a+sigma)*fac
      qq=qq+b_GA
      if (crit_GA(qq,b_GA,acc))80,17,17
   17 continue
      ncvg=ncvg+2
   80 return
c
      end


      subroutine TorGA_legz
cProlog

      implicit none
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8 
     & accc, acc, pisr, pi, r1, r2, r3, r4
      integer ncvg, ncv, nfrig
c
cSAP080617
      complex*16 
     & vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv,zz1,zz2,srz
c....&|--1---------2---------3---------4---------5---------6---------7-|
      common /TorGA_legbl/ vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv, zz1,
     &zz2,srz,acc,pisr,pi,r1,r2,r3,r4,ncvg,nfrig
cSAP080617
      complex*16 
     & f1, f2, zv, clogok_GA, rgam_GA
c     uses clogok_GA, rgam_GA, hypgm, trdz
c
      gr=clogok_GA(2.0d0*zz,nfrig,2)
      a=1.5d0
      b_GA=1.0d0
      zv = rgam_GA (vv, a, b_GA) * EXP (vv*gr)
      a=vv/2.0d0+1.0d0
      b_GA=vv/2.0d0+0.5d0
      c=vv+1.5d0
      accc=acc/100.0d0
      call TorGA_hypgm(a,b_GA,c,z2,f1,accc,ncvg)
      f1=f1/(2.0d0*zz*zv)
      qq=pisr*f1
      if (ABS (cvv) - 0.001d0) 10, 10, 9
c
c   trdz expects gr = clogok_GA(2.0*zz,nfrig,2) but destroys contents
c
   10 call TorGA_trdz
      go to 80
    9 a=-vv/2.0d0
      b_GA=(1.0d0-vv)/2.0d0
      c=0.5d0-vv
      call TorGA_hypgm(a,b_GA,c,z2,f2,accc,ncv )
      ncvg=ncvg+2*ncv
      f2=f2*zv/(vv+0.5d0)
      pp=(f1*svv/cvv+f2)/pisr
   80 return
c
      end


      subroutine TorGA_legor
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8 
     & pi, pisr, acc, r1, r2, r3, r4
      integer nfrig, ncvg
cSAP080617
      complex*16 
     & vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv,zz1,zz2,srz
c
c....&|--1---------2---------3---------4---------5---------6---------7-|
      common /TorGA_legbl/ vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv, zz1,
     &zz2,srz,acc,pisr,pi,r1,r2,r3,r4,ncvg,nfrig
c
cSAP080617
      complex*16 
     & rgam_GA
c     uses rgam_GA
c
      a=0.5d0*(1.0d0+vv)
      b_GA = 0.0d0
      c=0.5d0
      gr = rgam_GA (a, b_GA, c)
      c=u*pi*vv*0.5d0
      a = exp ( c)
      b_GA = exp (-c)
      if (nfrig) 9, 10, 11
    9 qq=pisr*0.5d0*gr*u*a
      go to 12
   10 qq=u*pisr/4.0d0*gr*(a-b_GA)
      go to 12
   11 qq=-pisr/2.0d0*gr*b_GA*u
   12 pp=gr*(a+b_GA)/(2.0d0*pisr)
      return
c
      end



      subroutine TorGA_legfn1 (v,z,p,ppr,q,nc,nf)
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
      integer nc, nf, nc1
c
c     calculates ppr, the derivative of p, as well as p and q.
cSAP080617
      complex*16 
     & v,z,p,ppr,q,vm1,pm,qm,c1,zsq
      data c1/(1.0d0,0.0d0)/
c     uses legfn
c
      call TorGA_legfn(v,z,p,q,nc,nf)
      vm1=v-c1
      call TorGA_legfn(vm1,z,pm,qm,nc1,nf)
      zsq=z*z
      if (zsq .eq. c1) go to 5
      ppr=(v/(zsq-c1))*(z*p-pm)
      return
    5 ppr=p*v*(v+c1)/2.d0
      return
c
      end




      subroutine TorGA_hypgm (a, b_GA, c, z, h, acc, ncvg)
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8 
     & crit_GA, acc
      integer ncvg, i
c
cSAP080617
      complex*16 
     & a,b_GA,c,z,h,aa,bb,cc,zz,add,dd,hh
c     uses crit_GA
c
      ncvg=0
      zz=z
      hh=(1.0d0,0.0d0)
      aa=a
      bb=b_GA
      cc=c
      add=(1.0d0,0.0d0)
      dd=(1.0d0,0.0d0)
      do 21 i=1,50
      add=add*zz*aa/dd*bb/cc
      hh=hh+add
      if (crit_GA(hh,add,acc))3,4,4
    4 aa=aa+1.0d0
      bb=bb+1.0d0
      cc=cc+1.0d0
      dd=dd+1.0d0
   21 continue
      ncvg=1
    3 h=hh
      return
c
      end




      subroutine TorGA_trdz
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)SAP080617
cSAP080617 
      real*8
     & ek, crit_GA, pi, pisr, acc, r1, r2, r3, r4
      integer k, kk, i, ncvg, nfrig
cSAP080617
      complex*16 
     & vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv,zz1,zz2,srz
c
c....&|--1---------2---------3---------4---------5---------6---------7-|
      common /TorGA_legbl/ vv,zz,pp,qq,z1,z2,u,a,b_GA,c,gr,cvv,svv, zz1,
     &zz2,srz,acc,pisr,pi,r1,r2,r3,r4,ncvg,nfrig
cSAP080617
      complex*16 
     & rgam_GA, psifun_GA
     & ,ep, ep2, pie, zzs2, aa, bb, cc, fac, sum, add
c     uses rgam_GA, crit_GA, psifun_GA
c
c   trdz  expects gr set = clogok_GA(2.0*zz,nfrig,2)
c
      k = dble(vv+1.0d0)
      ek = k
      ep = vv+0.5d0-ek
      zzs2 = z2/4.0d0
      kk = k-1
      sum = 0.0d0
      if (k) 9, 9, 10
c
   10 b_GA = 0.5d0
      c = 1.0d0
      fac = rgam_GA (vv, b_GA, c) * EXP (vv*gr)
      sum = fac
      if (kk) 9, 11, 12
   12 aa = -vv-2.0d0
      bb = -vv-0.5d0
      do 2 i=1,kk
      bb = bb+1.0d0
      aa = aa+2.0d0
      fac = fac*(1.0d0+aa)*aa*zzs2/(bb*FLOAT (i))
      sum = sum+fac
      if (crit_GA(sum,fac,acc)) 51, 2, 2
    2 continue
c
      go to 11
c
    9 sum = 0.0d0
   51 continue
   11 pie = pi*ep
      ep2 = ep/2.0d0
      a = ek
      b_GA = 0.5d0-ep
      c = 1.0d0
      bb = 0.0d0
      fac = -(1.0d0-pie**2/3.0d0)*rgam_GA(a, b_GA, c)*rgam_GA(c,bb,-ep)/
     & (pi*EXP ((2.0d0*ek-vv)*gr))
      a = ek+0.5d0
      b_GA = 1.0d0-ep2
      c = ek+1.0d0+ep2
      add = 2.0d0 * psifun _GA(a) - psifun_GA (b_GA) - 
     &      psifun_GA (c) - 2.0d0 * gr
      gr = fac*add*(1.0d0+ep2*add)
      sum = sum+gr
      aa = ek-1.5d0-ep
      bb = -ep
      cc = ek
c
      do 22 i=1,50
      aa = aa+2.0d0
      bb = bb+1.0d0
      cc = cc+1.0d0
      fac = fac*(1.0d0+aa)*aa*zzs2/(bb*cc)
      add = add+2.0d0/(a+1.0d0)+2.0d0/a-1.0d0/b_GA-1.0d0/c
      a = a+2.0d0
      b_GA = b_GA+1.0d0
      c = c+1.0d0
      gr = fac*add*(1.0d0+ep2*add)
      sum = sum+gr
      if (crit_GA(sum,gr,acc)) 52, 22, 22
   22 continue
c
      ncvg = ncvg+4
   52 pp = sum/pisr

c      write(*,*)'TorGA_trdz pp,sum,pisr', pp,sum,pisr

      return
c
      end

!=======================================================================
!=======================================================================
      real*8 function fcn_GA (expe)
      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8
     &alfh, beta, signv, expe, fu, fupr, etam, hh, hpr, 
     &epar, rl, rmin, rmax, rabs, charge, asp1, zhat, ap0, alfsq, aspect
     &, elooro, r3, eparsq, eta, rlams, pas, gam, psq, gamsq
      integer l, modl2,lm1, lp1, iheat
c
c added 20 Jan 92 by Joe Freeman
cSAP080617
      complex*16
     & alpha 
      common /TorGA_fcncom/ aspect,epar,rl,rmin,rmax,rabs,etam,charge, 
     &asp1,zhat,ap0,alfsq,elooro,r3, eparsq,l,iheat,lm1,lp1,modl2
      common /TorGA_mdl3com/alpha,rlams,pas
      common /TorGA_energy/gam,eta,psq,gamsq
c

      !write(*,*)'in  function fcn_GA (expe)'

      call TorGA_coefc(alfh,beta,signv,expe)

cSAP080617
!      write(*,*)'fcn_GA: after TorGA_coefc alfh,beta',
!     &alfh,beta

      call TorGA_ffp(fu,fupr)

cSAP080617
       !write(*,*)'fcn_GA: after TorGA_ffp fu,fupr',fu,fupr

c      if (modl2) 10,10,20
c   10 hh=hfun(eta,hpr)
c      fcn=(alfh*fupr*hh+beta*fu*hpr)*signv
c      return
c
cSAP080617
c      write(*,*)'eta,etam',eta,etam

   20 if (eta-etam) 22,30,30
   22 call TorGA_geth3(hh,hpr,eta)

      !write(*,*)'fcn_GA: after 22 TorGA_geth3 hh,hpr,eta',hh,hpr,eta

      fcn_GA=(alfh*fupr*hh+beta*fu*hpr)*signv

cSAP080617
!      write(*,*)'fcn_GA: fcn_GA,alfh,fupr,hh,beta,fu,hpr,signv',
!     &fcn_GA,alfh,fupr,hh,beta,fu,hpr,signv
  
      return
   30 fcn_GA = 0.0d0
      return
      end function fcn_GA
!=======================================================================
!=======================================================================

cSAP080617
      real*8 function fcd (x)
cProlog

c
      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8
     & gamma, theth, emin, usq, uzsq, elom, uxsq, aspect,
     & epar, rl, rmin, rmax, rabs, etam, charge, asp1, zhat, ap0, alfsq,
     & elooro, r3, eparsq, ra0, elom2, tor, gammin, agrel, bgrel, qgrel,
     & enpar2
      integer l, modl2,lm1, lp1, iheat, ixo, kg
cSAP080617
      real*8 x
c
      common /TorGA_fcncom/ aspect,epar,rl,rmin,rmax,rabs,etam,charge, 
     &asp1,zhat,ap0,alfsq,elooro,r3, eparsq,l,iheat,lm1,lp1,modl2
      common /TorGA_rela/ enpar2,elom,theth,emin,ra0,elom2,tor, gammin,
     &agrel,bgrel,qgrel,ixo,kg
c
      gamma=1.d0+theth*(-log(x)+emin)
      usq=gamma**2-1.d0
      uzsq=(gamma-elom)**2/enpar2
      uxsq=usq-uzsq
      fcd=uxsq**l
      return
      end

cSAP080617
      real*8 function fjpd0_GA (epar, lh, sigma, zeff)
cProlog

      implicit none
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8
     & rutpi,erfcb_GA,sum,erffac, epar, rlfac, rlk, sum1,
     & ckj_GA, sigma, zeff, dz, dzp

      integer ifac_GA, k1, lh, k, j1, j


c
c     calculates nonrelativistic ecrh current drive efficiency for no
c     trapped particles.
c     epar is resonant eparallel/th; lh is harmonic number;
c     sigma is (lh*cyclotron freq/wave freq) - 1.
c
      data rutpi/1.772453851d0/
c
      sum = 0.0d0
      erffac=rutpi*erfcb_GA(0,dz,dzp,epar)/sqrt (epar)
      rlfac = FLOAT (ifac_GA(lh))
      do 20 k1=1,lh+1
      k=k1-1
      rlk=rlfac/(ifac_GA(k)*ifac_GA(lh-k))
      sum1 = 0.0d0
      do 10 j1=1,k1
      j=j1-1
      sum1=sum1+(3.0d0*ckj_GA(k,j)+sigma*ckj_GA(k1,j))*epar**(-j)
   10 continue
      sum1=sum1+sigma*ckj_GA(k1,k1)*epar**(-k1)
      sum=-sum+rlk*(epar*sum1+epar**(-k)*erffac*(3.0d0*epar*ckj_GA(k,k1)
     &+sigma*ckj_GA(k1,k1+1)))
   20 continue
   35 fjpd0_GA=epar**lh*sum/(ifac_GA(lh)*(5.d0+zeff))
      return
c
      end
cSAP080617
      real*8 function fjpd0h_GA (epar, zeff)
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ))
cSAP080617
      real*8
     &rutpi,erfcb_GA,rtepar,erffac, dz, dzp, zeff, epar

c
c --- calculate nonrelativistic lower hybrid efficiency for no trapped particles
c
      data rutpi/1.772453851d0/
c
      rtepar = sqrt (epar)
      erffac=rutpi*erfcb_GA(0,dz,dzp,epar)
      fjpd0h_GA=(epar*4.d0+1.5d0+(1.5d0*rtepar+.75d0/rtepar)*erffac)/
     &(5.d0+zeff)
      return
c
      end
      subroutine TorGA_gammas (gam1, gam2, eta)
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8
     & omre, ra0, eta, rutfac, enpar2, elom2, rut, denomi
     &, elom, gam1, gam2, gamd, emin, tor, gammin, agrel, bgrel, qgrel, 
     &theth
      integer ixo, kg
c
c     calculates roots of gamma(eta)
c
      common /TorGA_rela/ enpar2,elom,theth,emin,ra0,elom2,tor, gammin,
     &agrel,bgrel,qgrel,ixo,kg
c
      omre=1.0d0-ra0*eta
      rutfac=omre*enpar2*(omre*enpar2-1.0d0+elom2)
      if (rutfac .lt. 0) go to 40
      rut = sqrt (rutfac)
      denomi=1.0d0/(1.0d0-omre*enpar2)
      gam1=(elom-rut)*denomi
      gam2=(elom+rut)*denomi
c     arrange roots in increasing order:
      if (gam2 .ge. gam1) return
      gamd=gam2
      gam2=gam1
      gam1=gamd
      return
c
c     if no physical roots, returns gam1=gam2=-1.23456789e10
c
   40 gam1=-1.23456789d10
      gam2=-1.23456789d10
      return
c
      end


cSAP080617
      integer function ninter (x)
c      function ninter (x)

cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8
     & y, x

cSAP080617
c      integer ninter, n
      integer n
c
c --- This function was formerly named NINT, but changed on 13 May 1999
c --- by Joe Freeman to avoid a name conflict with the Fortran intrinsic
c --- function of the same name. The name henceforth is NINTER.
c
c     x = 1/2 of 4n+1 set ninter = 1
c     x = 1/2 of 4n+2 set ninter = 2
c     x = 1/2 of 4n+3 set ninter = 3
c     x = 1/2 of 4n+4 set ninter = 4
c     otherwise ..... set ninter = 5
c
      y = x + x
      n = INT (y)
      if (y-FLOAT (n)) 9, 10, 9
    9 ninter = 5
      go to 80
   10 ninter = MOD (n, 4)
      if (ninter) 12, 11, 80
   12 ninter = 4 + ninter
      go to 80
   11 ninter = 4
   80 return
c
      end

cSAP080617
      complex*16 function clogok_GA (z, n, m)
cProlog

      implicit none
ccSAP080617
      complex*16
     & z,sz
c explicit type declaration 7/12/02 RAJ
cSAP080617
      real*8
     & pi, s

      integer nf, n, m

      data pi /3.1415926535898d0/
c    
      sz = LOG (z)
      s = dble(sz)
      nf = n + 5
      go to (21,22,23) ,m
   21 go to (2,2,2,11,10,9,2,2,2),nf
   22 go to (2,11,2,2,2,2,2,9,2),nf
   23 go to (2,11,2,2,2,2,2,9,2),nf
    9 sz = dcmplx(s,-pi)
      go to 2
   10 sz = dcmplx(s,0.0d0)
      go to 2
   11 sz = dcmplx(s,pi)
    2 clogok_GA=sz
   80 return
c
      end

cSAP080617
      complex*16 function csqrtk_GA (z, n, m)
cProlog

      implicit none
c
cSAP080617
      complex*16
     & z, sz
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8 s 

      integer nf, n, m

c
      nf = n + 5
      sz = sqrt (z)
      s = ABS (dimag(sz))
      go to (21,22),m
   21 go to (9,2,2,9,11,11,2,2,11),nf
   22 go to (2,11,2,2,2,2,2,9,2),nf
    9 sz = dcmplx(0.0d0,-s)
      go to 2
   11 sz = dcmplx(0.0d0,s)
    2 csqrtk_GA = sz
      return
c
      end


cSAP080617
      complex*16 function rgam_GA (z, a, b_GA)
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8
     & xa, xb, ya, yb

      integer la, lb

c
      integer iotty
      data iotty /6/
c
cSAP080617
      complex*16
     &z, a, b_GA, za, zb, clogam_GA
c     uses clogam_GA
c
      za=z+a
      zb=z+b_GA
      xa = -dble(za)
      xb = -dble(zb)
      ya = dimag(za)
      yb = dimag(zb)
      if (xa .eq. xb .and. ya .eq. yb) go to 2
      la=0
      lb=0
      if (ya .ne. 0.d0 .or. xa .ne. ABS (AINT (xa))) la=1
      if (yb .ne. 0.d0 .or. xb .ne. ABS (AINT (xb))) lb=1
      if (la .eq. 1 .and. lb .eq. 1) go to 1
      rgam_GA = 0.0d0
      if (la .eq. 1 .and. lb .eq. 0) return
      write (iotty, 10) z, a, b_GA
   10 format (' RGAM ... is not defined or infinite for z = ',2e12.4, 
     &5h a = , 2e12.4, 8h b_GA = , 2e12.4)
      return
c
    2 rgam_GA = 1.0d0
      return
c
    1 rgam_GA = EXP (clogam_GA(za)-clogam_GA(zb))
      return
c
      end


cSAP080617
      complex*16 function psifun_GA (z)
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8   
     & pi, b_GA, x, a

      integer n, i

ccSAP080617
      complex*16
     &z, u, v, h, r
      dimension b_GA(6)
c
      integer iotty
      data iotty /6/
c
      data pi / 3.141592653589793d0/
      data b_GA / +8.333333333333333d-2, -8.333333333333333d-3, +3.
     &968253968253968d-3, -4.166666666666667d-3, +7.575757575757576d-3, 
     &-2.109279609279609d-2 /
c
      u = z
      x = dble(u)
      a = ABS (x)
      if (dimag(u) .eq. 0.0d0 .and. -a .eq. AINT (x)) go to 4
      if (x .lt. 0.0d0) u = -u
      v = u
      h = 0.0d0
      if (a .ge. 15.0d0) go to 3
      n = 14 - INT (a)
      h = 1.0d0 / v
      if (n .eq. 0) go to 2
      do 1 i = 1,n
      v=v+1.0d0
    1 h=h+1.0d0/v
    2 v=v+1.0d0
    3 r = 1.0d0 / v**2
      psifun_GA=LOG(v)- 0.5d0/v - 
     & r * (b_GA(1)+r*(b_GA(2)+r*(b_GA(3)+r*(b_GA(4)+r*(
     &b_GA(5)+r*(b_GA(6)+r*b_GA(1))))))) - h
      if (x .ge. 0.0d0) return
      h = pi * u
      psifun_GA = psifun_GA + 1.0d0 / u + pi * COS (h) / SIN (h)
      return
c
    4 write (iotty, 100) x
  100 format (' PSIFUN ... argument is non-positive integer = ',f20.2)
      psifun_GA = 0.0d0
      return
c
      end

cSAP080617
      real*8 function crit_GA (sum, del, accs)
cProlog

      implicit none
c
cSAP080617
      complex*16
     & sum, del
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8 accs

      crit_GA=dble(del)**2 + dimag(del)**2 - accs * (dble(sum)**2 +dimag
     &(sum)**2)
      return
c
      end

      subroutine TorGA_coefc (alfh, beta, signv, expe)
cProlog

      implicit none
c
c     calculates coefficients for current-drive integral.
c
c explicit type declaration 7/12/02 RAJ
cSAP080617
      real*8
     & et, e, emin, signv, gam, theth, elom, gamsq, psq, 
     &gme, gmesq, pparsq, enpar2, eta, ra0, alfh, beta, tor, expe, 
     &aspect, epar, rl, rmin, rmax, rabs, etam, charge, asp1, zhat, ap0,
     & alfsq, elooro, r3, eparsq, elom2, gammin, agrel, bgrel, qgrel
      integer l, iheat,lm1, lp1, kg, modl2, ixo
c
      common /TorGA_fcncom/ aspect,epar,rl,rmin,rmax,rabs,etam,charge, 
     &asp1,zhat,ap0,alfsq,elooro,r3, eparsq,l,iheat,lm1,lp1,modl2
      common /TorGA_rela/ enpar2,elom,theth,emin,ra0,elom2,tor, gammin,
     &agrel,bgrel,qgrel,ixo,kg
      common /TorGA_energy/ gam,eta,psq,gamsq
c
      et = -LOG (expe)
      e=et+emin
      signv = 1.0d0
      gam=1+theth*e
      if (gam .lt. elom) signv=-1.0d0
      gamsq=gam*gam
      psq=gamsq-1.0d0
      gme=gam-elom
      gmesq=gme*gme
      pparsq=gmesq/enpar2
      eta=(1.0d0-pparsq/psq)/ra0
!      write(*,*)' TorGA_coefc: expe,et,psq-pparsq,l=',
!     &                         expe,et,psq-pparsq,l
   40 if (l) 44,44,42
   42 alfh=(psq-pparsq)**l
      go to 46
   44 alfh=1.0d0
   46 if (ixo .eq. 1) alfh = gmesq*alfh
      if (iheat) 60,55,55
   55 beta=tor*alfh*gme*(-1.0d0+gam*gme/(psq*enpar2))/psq
      return
   60 beta=-2.0d0*alfh*eta*gam/psq
      return
c
      end subroutine TorGA_coefc


cSAP080617
      real*8 function erfcb_GA (j, z, zp, zhi)
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
      real*8
     &zero, pii, pin, dd, a, b_GA, c, d, e, p, a2, b2, c2, 
     &d2, e2, zzz, t, z, zp, zar, zhi, zfb

      integer ivar, j

c
c     normalization*EXP (zhi)*Integral(t**m*EXP (-t^2)) from
c     sqrt (zhi) to infty;  normalization is such that erfcb(,,,0)=1.
c     if j>0, need to dimension z, zp in calling program; then
c     elements of z, zp are integrals for all m .le. 2*(j+1)/2.
c     erfcb is complementary error function; z is above function;
c     zp is derivitive.
c
cSAP080617
      complex*16
     & zarg, cim
      dimension z(*), zp(*)
      data ivar/0/, zero/0.0d0/
      data pii, pin, dd, a, b_GA, c, d, e, p, a2, b2, c2, d2, e2/1.
     &12837916709551d0 , 0.56418958354776d0 , 0.66666666666667d0 , 0.
     &24598329392237d0 , 0.22260371336736d0 , 0.12873395368514d0 , 0.
     &27281029618672d0 , -.11611455108396d0 , 0.43599403657113d0 , 0.
     &49196658784475d0 , 0.66781114010207d0 , 0.51493581474056d0 , 1.
     &3640514809336d0 , -.69668730650377d0 /
      data cim/(0.0d0,1.0d0)/
c
      zzz = zhi
      if (zzz .le. zero) go to 2
      if (ivar .gt. 0) go to 1
      zarg = sqrt (zzz)
      t = 1.0d0/(1.0d0+p*zarg)
      zfb = t*(a+t*(a+t*(b_GA+t*(c+t*(d+t*e)))))
      erfcb_GA = zfb
      if (j .eq. 0) return
      z(1) = zfb+pii*zarg
      if (j .eq. 2) return
      zp(1) = 0.5d0*(pii-p*t**2*(a+t*(a2+t*(b2+t*(c2+t*(d2+e2*t))))))/
     &zarg
      if (j .lt. 3) return
      z(2)=z(1)+pii*zarg*zzz*dd
      if (j .eq. 4) return
      zp(2)=zp(1)+pii*zarg
      return
c
    1 zar = sqrt (zzz)
      zarg=cim*zar
c***  zfb=zprime(zarg,zeta,zdprime)
c***  zfb=-zeta*cim*pin
      erfcb_GA=zfb
      if (j .eq. 0) return
      z(1) = zfb+pii*zar
      zp(1)=zfb
      if (j .le. 2) return
      z(2)=z(1)+dd*pii*zzz*zar
      zp(2)=z(1)
      return
c
    2 zfb = EXP (zzz)
      erfcb_GA = zfb
      z(2)=zfb
      z(1)= z(2)
      zp(2)=zfb
      zp(1)=zp(2)
      return
c
      end


cSAP080617
      integer function ifac_GA (l)
c      function ifac_GA (l)
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
cSAP080617
c      integer i, l, ifac_GA
      integer i, l

c
      ifac_GA=1
      if (l .lt. 2) return
      do 15 i=2,l
   15 ifac_Ga=ifac_GA*i
      return
c
      end

cSAP080617
      real*8   function ckj_GA (k, j)
c      function ckj_GA (k, j)
cProlog

      implicit none
c
c     calculates (.5+k)!/(.5+k-j)!
c
cSAP080617
c      real*8
c     & ckj_GA
      integer i,j,k
      ckj_GA = 1.0d0
      if (j) 5, 5, 10
    5 return
   10 do i=1,j
      ckj_GA=ckj_GA*(1.5d0+k-i)
      enddo
      return
c
      end

cSAP080617
      complex*16 function clogam_GA (z)
cProlog

      implicit none
ccSAP080617
      complex*16
     & z, v, h, r
      dimension b_GA(10)
c
      integer iotty
      data iotty /6/
c
c explicit type declaration 7/12/02 RAJ
cSAP080617
      real*8
     & pi, b_GA, x, t, f, c, d, a, e

      integer i, n

      data pi /3.141592653589793d0/
      data b_GA /+8.33333333333333d-2, -2.77777777777778d-3, +7.
     &93650793650794d-4, -5.95238095238095d-4, +8.41750841750842d-4, -1.
     &91752691752692d-3, +6.41025641025641d-3, -2.95506535947712d-2, +1.
     &79644372368831d-1, -1.39243221690590d+0/
c
      x = dble(z)
      t = dimag(z)
      if (-ABS (x) .eq. AINT (x) .and. t .eq. 0.0d0) go to 5
      f = ABS (t)
      v = dcmplx(x,f)
      if (x .lt. 0.0d0) v = 1.0d0 - v
      h = (0.0d0, 0.0d0)
      c = dble(v)
      if (c .ge. 7.0d0) go to 3
      n = 6 - INT (c)
      h=v
      d = dimag(v)
      a = ATAN2 (d,c)
      if (n .eq. 0) go to 2
      do 1 i = 1,n
      c=c+1.0d0
      v = dcmplx(c,d)
      h=h*v
    1 a = a + ATAN2 (d,c)
    2 h = dcmplx(0.5d0*LOG (dble(h)**2 + dimag(h)**2),a)
      v = v + 1.0d0
    3 r = 1.0d0 / v**2
      clogam_GA=0.918938533204673d0+(v-0.5d0)*LOG(v)-v+
     & (b_GA(1)+r*(b_GA(2)+r*
     &(b_GA(
     &3)+ r*(b_GA(4)+r*(b_GA(5)+r*(b_GA(6)+r*(b_GA(7)
     & +r*(b_GA(8)+r*(b_GA(9)+r*b_GA(10))))))))
     &)) /v-h
      if (x .ge. 0.0d0) go to 4
c
      a = AINT (x) - 1.0d0
      c=pi*(x-a)
      d=pi*f
      e = EXP (-2.0d0*d)
      f = SIN (c)
      e=d+0.5d0*LOG (e*f**2+0.25d0*(1.0d0-e)**2)
      f = ATAN2 (COS (c)*TANH (d),f)-a*pi
      clogam_GA = 1.144729885849400d0 - dcmplx(e,f)-clogam_GA
c
    4 if (t .lt. 0.0d0) clogam_GA = dconjg(clogam_GA)
      return
c
    5 write (iotty, 100) x
      clogam_GA = (0.0d0, 0.0d0)
      return
  100 format (' clogam ... argument is non-positive integer = ',f20.2)
c
      end

cSAP080617
      real*8 function qgauleg_GA (func, a, b_GA, n, ifail)
c      function qgauleg_GA (func, a, b_GA, n, ifail)
cProlog

c
      implicit none
c
cSAP080617
      real*8
c     & qgauleg_GA,func,a,b_GA
     & func,a,b_GA
      integer n,ifail
      external func
c
c     returns n-point gauss-legendre quadrature of 'func'
c     from a to b.  n should be less than nmax =64
c ----------------------------------------------------------------------
c
      integer nmax
      parameter (nmax=64)
      integer stderr
      parameter (stderr=6)
      integer nn,ns,j
cSAP080617
      real*8
     & x(nmax),w(nmax),xx(nmax),y(nmax),ss,xm,xr
      character errmsg*80
      save ns,x,w,errmsg
      data ns/-1/
      data errmsg/'WARNING: from qgauleg_GA: n is set to nmax =64'/
c

c      write(*,*)'in qgauleg_GA   a, b_GA, n', a, b_GA, n
      
      nn=n

c      write(*,*)'nmax,nn',nmax,nn

      if (nn .gt. nmax) then
      nn=nmax
      ifail=1
      write (stderr,'(a)')errmsg
      endif

c      write(*,*)'ns,nn',ns,nn

      if (ns .ne. nn) then
      ns=nn
      call gauleg_GA(x,w,nn)

c      write(*,*)'after gauleg_GA'
c      write(*,*)'x',x
c      write(*,*)'w',w

      endif
c
      xm=0.5d0*(b_GA+a)
      xr=0.5d0*(b_GA-a)

cyup      write(*,*)'xm,xr',xm,xr

      do j=1,nn
      xx(j)=xm+xr*x(j)
      y(j)=func(xx(j))

      !write(*,*)'qgauleg_GA: j,xx(j),y(j)',j,xx(j),y(j)

      enddo
      ss=0.d0
      do j=1,nn
      ss=ss+w(j)*y(j)

c      write(*,*)'j,w(j),y(j),ss',j,w(j),y(j),ss

      enddo
      qgauleg_GA=ss*xr

c      write(*,*)'qgauleg_GA,ss,xr',qgauleg_GA,ss,xr

      return
c
      end function qgauleg_GA



      subroutine gauleg_GA (x, w, n)
cProlog

c
      implicit none
c
      integer n
cSAP080617
      real*8
     & x(n),w(n)
c     returns the abscissas and weights for n-point
c     gauss-legendre integration
c ----------------------------------------------------------------------
cSAP080617
      real*8
     &pi,eps
cSAP
c      parameter (pi=3.141592654d0,eps=3.d-14)
      parameter (eps=3.d-14)
      integer i,j,m
cSAP080617
      real*8
     & z,p1,p2,p3,pp,z1
c

c      write(*,*)'gauleg_GA'

cSAP080917
      pi=4.d0*atan(1.d0)

      m=(n+1)/2
      do i=1,m
      z=cos(pi*(i-0.25d0)/(n+0.5d0))

cSAP080916
c      write(*,*)'i,z',i,m

    1 continue
      p1=1.d0
      p2=0.d0
      do j=1,n
      p3=p2
      p2=p1
      p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
      enddo

cSAP080616
c      write(*,*)'p1,p2,p3',p1,p2,p3

      pp=n*(z*p1-p2)/(z*z-1.d0)
      z1=z
      z=z1-p1/pp

cSAP080617
c      write(*,*)'pp,z,z1',pp,z,z1

      if (abs(z-z1) .gt. eps) go to 1
      x(i)=-z
      x(n+1-i)=z
      w(i)=2.d0/((1.d0-z*z)*pp*pp)
      w(n+1-i)=w(i)

cSAP080617
c      write(*,*)'x(i),x(n+1-i)',x(i),x(n+1-i)
c      write(*,*)'w(i),w(n+1-i)',w(i),w(n+1-i)

      enddo
c
      return
c
      end
