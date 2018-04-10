      subroutine curba (rjpd,rjpd0,ratjpd,denom,aspct,enpar,tc,thtc,
     1                   theta,elomom,lh,zeff,model,tol,n0,ig)
c
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
c     4 pi n e**4 log lambda / m**2 vt**3, and vt = SQRT (2T/m).  Note
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
c      h(2,)=b*dh(1,)/d eta, where b = <SQRT (1-mirror
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
c
* NAG
      integer nvar,ngmax
      parameter (nvar=2,ngmax=200)
      integer lw,liw
      parameter (lw=nvar*(3*nvar*ngmax+4*nvar+1)+ngmax*(7*nvar+2))
      parameter (liw=2*nvar*ngmax+nvar+ngmax)
      integer ip,ir,iy,iwork(liw),ier
      real par(5),work(lw),abt(nvar)
      dimension c(2,2),d(2,2),gam(2),xm(200),h(2,200)
      common /bdryc/c,d,gam
      data mnp/200/,etmax0/20./
* NAG dimension c(2,2),d(2,2),gam(2),xm(200),h(2,200),w(12000),iw(1600)
* NAG data mnp/200/,lw/12000/,liw/1600/,etmax0/20./
      common /fcncom/ aspect,epar,rl,l,rmin,rmax,rabs,etam,charge,
     .                modl2,asp1,zhat,ap0,alfsq,lm1,lp1,elooro,r3,
     .                eparsq,iheat
      common /hfcom / modela,h3fac,etamh
      common /gval  / h,xm,inp
      common /rela  / enpar2,elom,theth,emin,ra0,elom2,ixo,tor,
     .                gammin,agrel,bgrel,qgrel,kg

* NAG external fcnf,fcng,bazd01,fcn
      external fcnf,fcng,fcn,fcni,fcnj,fcnb
* NAG
c
      write(*,*)'curba denom,aspct,enpar,tc,thtc',
     &denom,aspct,enpar,tc,thtc
      write(*,*) 'theta,elomom,lh,zeff,model,tol,n0,ig',
     &            theta,elomom,lh,zeff,model,tol,n0,ig

      kg=ig-2
      aspect=aspct
      sgnj=1.0
      if (enpar .lt. 0) sgnj=-1.0
      charge=zeff
      costh=COS (theta)
      rabs=1.0+aspect*costh
      rmax=1.0+aspect
      rmin=1.0-aspect
      etam=rmin/rmax
      etamh=etam
      et0 = 0.0
      thtcu=thtc
c
c     use th=tc if doing buckets
c
      if (thtc .le. 0.0) thtcu=1.0
      theth=tc*thtcu/511.0
      if (model) 5,7,7
    5 iheat=-1
      elom=0
      modela=-model
      go to 8
    7 iheat=1
      elom=elomom
      modela=model
    8 ixo=0
      l = IABS (lh)
c
c     set flag to get ppar**(2*ixo) factor in diffusion coef.;
c     for ecrh, ixo=0 for E-, E+ and =1 for E_par contributions
c
      if (lh .lt. 0) ixo=1
      modl2=modela-2
      call setpar
c
c --- Skip h calculation if nothing is new:
c
c      write(*,*)'in curba before ijump aspct,zeff,modela,tol,inp',
c     1                                aspct,zeff,modela,tol,inp

      ijump=kjump(aspct,zeff,modela,tol,inp)
      if (ijump .eq. 1)  go to 20
      call setcdg(c,d,gam)
      tolt=tol
      if (modela-3) 10,20,20
   10 n=2
      np=inp
      if (np) 18,18,12
   12 do 14 i=1,np-1
        xm(i)=(i-1)*(etam-et0)/(np-1)
   14   continue
      xm(np)=etam-et0
   18 ifail = 000
* NAG
      if (c(2,1) .eq. 0) then
         ip=1
         ir=10
      else
         ip=0
         ir=1
      end if
c
      iy=2
      par(1)=1.0
      par(2) = 0.0
      par(3) = 0.0
      par(4)=1.0            ! input initial grid
      par(5)=1.0            ! linear problem
      call dvcpr (n,fcni,fcnj,fcnb,et0,etam,mnp,np,ip,ir,tolt,xm,h,iy,
     .            abt,par,work,iwork,ier)
c
c     determine array for green's function
c
* NAG   call d02gbE (et0,etam,n,tolt,fcnf,fcng,c,d,gam,mnp,xm,h,np,w,
* NAG.               lw,iw,liw,ifail)
c
c     set up spline for calculation of interpolated values of h:
c
      call setara(h,xm,np)
c
c     do current-drive integral:
c
   20 ifail=0
      h3fac=4.0*zhat/(5.+charge)
c
c     CAN COMMENT OUT STATEMENTS FROM HERE TO RETURN if YOU ONLY WANT
c     GREEN'S FN, NOT CURRENT-DRIVE EFFICIENCY.
c --- set limits of integration for expe = EXP (-epst), where epst=e-emin,
c     normalized to th.
c
      expeu=1.0
      enpar2=enpar*enpar
      call limits(epst1,epst2,etmax)
c
c     set j/pd to zero if no resonance:
c
      if (etmax) 30,30,25
   25 if (thtc .lt. 0)  go to 45
      if (epst1 .gt. 0. .and. epst2 .lt. etmax)  go to 27
      n = MIN0 (64, n0)
c
c --- calculate integral for resonance which doesn't enter trapped region
c
      if (epst2 .le. 0.0) expel = EXP (-etmax)
      if (epst2 .gt. 0.0 .and. epst1 .le. 0.0)
     .  expel = EXP (-MIN (epst2, etmax))
      if (epst2 .gt. 0.0 .and. epst1 .gt. 0.0)
     .  expel = EXP (-MIN (epst1, etmax))
* NAG rint = d01baE(bazd01,expel,expeu,n,fcn,ifail)
      rint = qgauleg(fcn,expel,expeu,n,ifail)
      go to 29
c
c --- calculate integral for resonance divided by trapped region:
c
   27 expel = EXP (-etmax)
      expe1 = EXP (-epst1)
      expe2 = EXP (-epst2)
      n = n0 / 2
* NAG rint=d01baE(bazd01,expel,expe2,n,fcn,ifail)
      rint=qgauleg(fcn,expel,expe2,n,ifail)
* NAG rint=rint+d01baE(bazd01,expe1,expeu,n,fcn,ifail)
      rint=rint+qgauleg(fcn,expe1,expeu,n,ifail)
c
c --- calculate denominator
c
   29 expeld = EXP (-etmax)
      call denomf(denom,expeld,etmax)
      rjpd = -0.125*sgnj*thtcu*rmax*rint/(denom*zhat)
   28 if (modela .eq. 3) rjpd=rjpd*h3fac
      go to 35
c     set J/Pd=0 if no resonance:
   30 rjpd = 0.0
      gammin=1.0
   35 sigma=1.0-elomom
c
c --- calculate non-relativistic, zero-trapped-particle limits
c
      eparnr = 0.5*(1.0-elom/gammin)**2/(theth*enpar2)
      sgnj0=SIGN (1.0,(gammin-elom)*enpar)
      if (iheat .eq. 1) rjpd0=sgnj0*fjpd0(eparnr,l,sigma,zeff)*thtcu
      if (iheat .eq. -1) rjpd0=sgnj0*fjpd0h(eparnr,zeff)*thtcu
      ratjpd=rjpd/rjpd0
      return
c
c     bucket lift calculation:
c
   45 erise=-thtc/511.0
      etal = 0.0
      sgnj0=SIGN (1.0,(gammin-elom)*enpar)
      gamma1x=gammin+erise
      etau=MIN (1.0,2.0*(rabs/rmax)*elom*erise/(gamma1x*gamma1x-1.0))
      call getdfh(dfh,gammin,gamma1x,etal,etau)
      rjpd=-(511.0/(8.*zhat*tc))*sgnj0*rmax*dfh/erise
      go to 28
c
      end


      function fjpd0 (epar,lh,sigma,zeff)
c
c     calculates nonrelativistic ecrh current drive efficiency for no
c     trapped particles.
c     epar is resonant eparallel/th; lh is harmonic number;
c     sigma is (lh*cyclotron freq/wave freq) - 1.
c
      data rutpi/1.772453851/
      sum = 0.0
      erffac=rutpi*erfcb(0,dz,dzp,epar)/SQRT (epar)
      rlfac = FLOAT (ifac(lh))
      do 20 k1=1,lh+1
        k=k1-1
        rlk=rlfac/(ifac(k)*ifac(lh-k))
        sum1 = 0.0
        do 10 j1=1,k1
           j=j1-1
           sum1=sum1+(3.0*ckj(k,j)+sigma*ckj(k1,j))*epar**(-j)
   10      continue
        sum1=sum1+sigma*ckj(k1,k1)*epar**(-k1)
        sum=-sum+rlk*(epar*sum1 + epar**(-k)*erffac*(3.0*epar*
     .      ckj(k,k1)+sigma*ckj(k1,k1+1)))
   20   continue
   35 fjpd0=epar**lh*sum/(ifac(lh)*(5.+zeff))
      return
c
      end


      subroutine denomf (denom,expel,etmax)
c
c --- calculates normalized absorbed power
c
      common /fcncom/ aspect,epar,rl,l,rmin,rmax,rabs,etam,charge,
     .                modl2,asp1,zhat,ap0,alfsq,lm1,lp1,elooro,r3,
     .                eparsq,iheat
      common /rela  / enpar2,elom,theth,emin,ra0,elom2,ixo,tor,
     .                gammin,agrel,bgrel,qgrel,kg
c
      epsmin=emin*theth
      epst=etmax*theth
      epsmax=epst+epsmin
      ololdi = 0.0
      oldi=theth*(1.0-expel)
      ri=oldi
      ail=1.0-1.0/enpar2
      elm1=elom-1.0
      bil=2.0*(1.0+elm1/enpar2)
      cil=-elm1*elm1/enpar2
      capbu=(2.0*ail*epsmax+bil)*expel
      capbl=2.0*ail*epsmin+bil
c
c     set initial values for iteration:
c
      rj=theth*(theth-(epst+theth)*expel)
      rk=theth*(-epst*epst*expel+2.0*rj)
c
c     begin iteration---
c
      if (l .eq. 0)  go to 40
      do 30 j=1,l
         if (j-1) 10,10,15
   10    bfac=capbl-capbu
         go to 20
   15    bfac = 0.0
   20    ri=j*theth*theth*(bfac+(j-1)*(bil*bil-4.0*
     .    ail*cil)*ololdi+2.0*(2*j-1)*ail*oldi)
         if (ixo) 28,28,22
   22    rj=theth*((2*j+1)*ri-j*capbl*rj)
         rk=theth*(2.0*(j+1)*rj-capbl*j*rk)
   28    ololdi=oldi
         oldi=ri
   30 continue
   40 fac=gammin-elom
      denom=(1-ixo)*ri+ixo*(fac*fac*ri+2.0*fac*rj+rk)
      return
c
      end

      function fcn (expe)
c
      complex alpha                     ! added 20 Jan 92 by Joe Freeman
      common /fcncom/ aspect,epar,rl,l,rmin,rmax,rabs,etam,charge,
     .                modl2,asp1,zhat,ap0,alfsq,lm1,lp1,elooro,r3,
     .                eparsq,iheat
      common /mdl3com/alpha,rlams,pas
      common /energy/gam,eta,psq,gamsq
      call coefc(alfh,beta,signv,expe)
      call ffp(fu,fupr)
      if (modl2) 10,10,20
   10 hh=hfun(eta,hpr)
      fcn=(alfh*fupr*hh+beta*fu*hpr)*signv
      return
   20 if (eta-etam) 22,30,30
   22 call geth3(hh,hpr,eta)
      fcn=(alfh*fupr*hh+beta*fu*hpr)*signv
      return
   30 fcn = 0.0
      return
c
      end

      subroutine fcni (n, x, y, yp)
c
c --- uses subroutine 'fcnf' and 'fcng'
c
      implicit none
c
      integer    n
      real       x, y(n), yp(n)
      integer    nvar
      parameter (nvar = 2)
      real       f(nvar,nvar), g(nvar)
c
      call fcnf (x, f)
      call fcng (x, g)
      yp(1) = f(1,1)*y(1)+f(1,2)*y(2)+g(1)
      yp(2) = f(2,1)*y(1)+f(2,2)*y(2)+g(2)
      return
c
      end

      subroutine fcnj (n, x, y, pd)
c
c --- uses subroutine 'fcnf'
c
      implicit none
c
      integer    n
      real       x, y(n), pd(n,n)
      integer    nvar
      parameter (nvar = 2)
      real       f(nvar,nvar)
      integer    i, j
c
      call fcnf (x, f)
      do j=1,nvar
        do i=1,nvar
          pd(i,j) = f(i,j)
        end do
      end do
      return
c
      end

      subroutine limits (epst1,epst2,etmax)
c
c     calculate range of relativistic range of integration
c     etmax is emax-emin, normalized to thot.
c
      common /fcncom/aspect,epar,rl,l,rmin,rmax,rabs,etam,charge,
     . modl2,asp1,zhat,ap0,alfsq,lm1,lp1,elooro,r3,eparsq,
     . iheat
      common /rela/enpar2,elom,theth,emin,ra0,elom2,ixo,tor,
     . gammin,agrel,bgrel,qgrel,kg
      data etmax0/20./
      if (iheat) 60,60,5
c     ECH limits:
c     find intersections of resonance with ppar axis:
    5 call gamma1s(gam1,gam2,0.0)
      if (gam1-1.0) 10,10,20
   10 if (gam2 .le. 1.0)  go to 40
      gammin=gam2
      etmax=etmax0
      go to 30
   20 gammin=gam1
      etcrit=(gam2-gam1)/theth
      etmax=MIN (etmax0,etcrit)
c     find intersections of resonance with passing-trapped separatrix:
   30 call gamma1s(gam3,gam4,etam)
      epst1=(gam3-gammin)/theth
      epst2=(gam4-gammin)/theth
      emin=(gammin-1.0)/theth
      return
c
c --- no resonance if elom+enpar2-1 < 0; set etmax negative as a flag.
c
   40 etmax=-1.0
      return
c     Limits for Landau resonance:
   60 etmax=etmax0
      epst2=-1.0
      gammin=1.0/SQRT (1.0-1.0/enpar2)
      emin=(gammin-1.0)/theth
      return
c
      end

      subroutine setpar
c
c     set eta-independent parameters
c
      complex alpha,chalf,cwun,ci,zarg,cpas,q
      common /fcncom/aspect,epar,rl,l,rmin,rmax,rabs,etam,charge,
     . modl2,asp1,zhat,ap0,alfsq,lm1,lp1,elooro,r3,eparsq,iheat
      common /rela/enpar2,elom,theth,emin,ra0,elom2,ixo,tor,
     . gammin,agrel,bgrel,qgrel,kg
      common /mdl3com/alpha,rlams,pas
      data chalf/(-.5,0.0)/,cwun/(1.0,0.0)/,ci/(0.0,1.0)/
c
      rl=l
      ra0=rmax/rabs
      r3=3.0*ra0
      asp1=1.0-aspect
      zhat=(1.0+charge)/2.
      agrel=(1.0/(5.+charge))**2
      bgrel=1.0-agrel
      qgrel=1.0/bgrel
      ap0 = -0.5*SQRT (rmax/rmin)*rmax
      eparsq=epar*epar
      alfsq=-2.0*aspect/asp1
      lm1=l-1
      lp1=l+1
      elooro=elom/ra0
      elom2=elom*elom
c     parameters for model=3
      rlams = SQRT (1.0-etam)
      alffac=(4./zhat)-.25
      if (alffac) 10,10,20
   10 alpha=chalf+cwun*SQRT (-alffac)
      go to 30
   20 alpha=chalf+ci*SQRT (alffac)
   30 zarg = CMPLX (rlams, 0.0)
      nf=0   !BH080604, added for initialization, following SAP
             !BH080604, Should nf possibly be something else?
             !BH080504, Pletzer chose to initialize nf = 0 in legfn.
      if (modl2 .eq. 1)  call legfn(alpha,zarg,cpas,q,nc,nf)
      pas = REAL (cpas)
      tor=2.0/ra0
      return
c
      end

      subroutine getdfh (dfh,gammin,gamma1x,etal,etau)
c
c     calculates difference of final and initial greens functions for
c     bucket lift from etal,gammin to etau,gamma1x
c
      common /energy/gam,eta,psq,gamsq
      common /fcncom/aspect,epar,rl,l,rmin,rmax,rabs,etam,charge,
     . modl2,asp1,zhat,ap0,alfsq,lm1,lp1,elooro,r3,eparsq,
     . iheat
      call setrel(gammin)
      call ffp(fl,fupr)
      call setrel(gamma1x)
      call ffp(fu,fupr)
      if (modl2) 10,10,20
   10 hl=hfun(etal,hpr)
      if (etau .ge. etam)  go to 30
      hu=hfun(etau,hpr)
      dfh=fu*hu-fl*hl
      return
   20 call geth3(hl,hpr,etal)
      if (etau .ge. etam)  go to 30
      call geth3(hu,hpr,etau)
      dfh=fu*hu-fl*hl
      return
   30 dfh=-fl*hl
      return
c
      end

      function fjpd0h (epar, zeff)
c
c --- calculate nonrelativistic lower hybrid efficiency for no trapped particles
c
      data rutpi/1.772453851/
      rtepar = SQRT (epar)
      erffac=rutpi*erfcb(0,dz,dzp,epar)
      fjpd0h=(epar*4.+1.5+(1.5*rtepar+.75/rtepar)*erffac)/(5.+zeff)
      return
c
      end


      subroutine fcnb (n, ya, yb, f)
c
      implicit none
c
      integer        n
      real           ya(n), yb(n), f(n)
      real           c(2,2), d(2,2), gam(2)
      common /bdryc/ c,      d,      gam
c
      f(1) = c(1,1)*ya(1)+c(1,2)*ya(2)+d(1,1)*yb(1)+d(1,2)*yb(2)-gam(1)
      f(2) = c(2,1)*ya(1)+c(2,2)*ya(2)+d(2,1)*yb(1)+d(2,2)*yb(2)-gam(2)
      return
c
      end


      function kjump (aspct,zeff,modela,tol,inp)
c
      kjump=0
      if ( aspct .eq. aspcto .and. zeff .eq. zeffo .and.
     .    modela .eq. modelo .and.  tol .eq. tolo  .and.
     .                              inp .eq. inpo)  kjump = 1
      aspcto=aspct
      zeffo=zeff
      modelo=modela
      tolo=tol
      inpo=inp
      return
c
      end


      subroutine setara (h,xm,np)
c
c cubspl     from testlib             version   3          05/19/77
c [1]  This version was obtained from Gary Smith, attached to a message
c [1]  dated 11/29/90, to Tim Luce.  There are corrections relative to
c [1]  previous versions.  Incorporated in TORAY (BobH:1-21-91).
c     This is the CURBA3 relativistic, finite aspect ratio
c      current-drive package.  See comments at beginning of
c      subroutine CURBA for use.  To compile this package
c      as a standalone code under Basis, run COSMOS CCLCUR.  Need to
c      read from CFS, user 313, directory bcur3, the files CURBA3,
c      VAR.CUR, CCLCUR and PACK.IN.  From directory legendb, file
c      LIBLEG. Also, need to run lib BASIS;
c      then type lines: x ctl; n ctl; x ctl pack.ctl ctl.o.
c      Note in current BASIS, need to type "package cur" at first
c       cur> prompt to get going; also the old "generate noplot"
c       is replaced by typing "ctlplot=no" once; then just
c       "generate".
c
c     sets up array of inp values of h(1,), the pitch-angle-dependent part
c     of the Spitzer-Harm function with trapped particles.  inp set to 50 by
c     data declaration.  Input: aspct is r/R, zeff is effective charge.
* NAG dimension h(2,200),xm(200),yh(200),wrk(1216),rk(204),c(204)
* NAG common /splicm/rk,c,etmax,npd
* NAG data lck/204/, lwrk/1216/
* NAG npd=np
* NAG do 10 i=1,np
* NAG   yh(i)=h(1,i)
c**10   continue
* NAG etmax=xm(np)
* NAG call e01baE(np,xm,yh,rk,c,lck,wrk,lwrk,ifail)
* NAG return
c
      dimension h(2,200),xm(200),yh(200),rk(204),c(204)
      dimension wk1(200,200),wk2(200,200),wkr(200)
      dimension ipvt(200)
      common /splicm/rk,c,etmax,npd
c
      npd=np
      do i=1,np
        yh(i) = h(1,i)
      end do
      etmax=xm(np)
      call bsplcof(np,xm,yh,rk,c,200,wk1,wk2,wkr,ipvt,ifail)
      return
c
      end


      subroutine setcdg (c,d,gam)
c
      dimension c(2,2),d(2,2),gam(2),cd(2,2),dd(2,2),gamd(2)
      data cd/3*0.0,1.0/, dd/1.0,3*0./, gamd/0.0,1.0/
      do 20 i=1,2
         gam(i)=gamd(i)
         do 20 j=1,2
            c(i,j)=cd(i,j)
            d(i,j)=dd(i,j)
   20 continue
      call coefg(ca,cb,0.0)
      c(2,1)=ca
      return
c
      end

      subroutine fcnf (eta,f)
c
      dimension f(2,2)
      common /fcncom/aspect,epar,rl,l,rmin,rmax,rabs,etam,charge,
     .               modl2,asp1,zhat,ap0,alfsq,lm1,lp1,elooro,r3,eparsq,
     .               iheat
c
      do 5 i=1,2
      do 5 j=1,2
        f(i,j) = 0.0
    5 continue
      call coefg(ca,cb,eta)
      f(1,2)=1.0/cb
      if (ABS (eta)-1.0e-8) 10,10,20
   10 f(2,1) = -0.5*ap0
      f(2,2) = -0.5*ca/cb
      return
   20 f(2,2)=-1.0/eta
      f(2,1)=ca*f(2,2)
      return
c
      end

      subroutine fcng (eta,g)
c
      dimension g(2)
c
      g(1) = 0.0
      if (ABS (eta)-1.0e-8) 10,10,20
   10 g(2) = 0.0
      return
   20 g(2)=1.0/eta
      return
c
      end


      subroutine setrel (gamma1)
c
      common /energy/ gam, eta, psq, gamsq
c
      gam   = gamma1
      gamsq = gam * gam
      psq   = gamsq - 1.0
      return
c
      end

      subroutine geth3 (hh,hpr,eta)
c
c     gets angle part of green's function in Legendre-function (square-well)
c     approx.  Note that h calculated here needs to be multiplied by
c     4/(Z+5) to get h of UCRL-95813
c
      complex alpha,zarg,pa,ppr,q
      common /mdl3com/alpha,rlams,pas
      rlam = SQRT (1.0-eta)
      zarg = CMPLX (rlam,0.0)
      call legfn1(alpha,zarg,pa,ppr,q,nc,nf)
      rpa  = REAL (pa)
      rppr = REAL (ppr)
      hh=-rlam+rlams*rpa/pas
      hpr=(1.0-rlams*rppr/pas)*.5/rlam
      return
c
      end

      subroutine coefc (alfh,beta,signv,expe)
c
c     calculates coefficients for current-drive integral.
c
      common /fcncom/ aspect,epar,rl,l,rmin,rmax,rabs,etam,charge,
     .                modl2,asp1,zhat,ap0,alfsq,lm1,lp1,elooro,r3,
     .                eparsq,iheat
      common /rela  / enpar2,elom,theth,emin,ra0,elom2,ixo,tor,
     .                gammin,agrel,bgrel,qgrel,kg
      common /energy/ gam,eta,psq,gamsq
c
      et = -LOG (expe)
      e=et+emin
      signv = 1.0
      gam=1+theth*e
      if (gam .lt. elom) signv=-1.0
      gamsq=gam*gam
      psq=gamsq-1.0
      gme=gam-elom
      gmesq=gme*gme
      pparsq=gmesq/enpar2
      eta=(1.0-pparsq/psq)/ra0
   40 if (l) 44,44,42
   42 alfh=(psq-pparsq)**l
      go to 46
   44 alfh=1.0
   46 if (ixo .eq. 1)  alfh = gmesq*alfh
      if (iheat) 60,55,55
   55 beta=tor*alfh*gme*(-1.0+gam*gme/(psq*enpar2))/psq
      return
   60 beta=-2.0*alfh*eta*gam/psq
      return
c
      end

      function ckj (k,j)
c
c     calculates (.5+k)!/(.5+k-j)!
c
      ckj = 1.0
      if (j)  5, 5, 10
    5 return
   10 do i=1,j
        ckj=ckj*(1.5+k-i)
      end do
      return
c
      end


      subroutine gamma1s (gam1,gam2,eta)
c
c     calculates roots of gamma1(eta)
c
      common /rela/ enpar2,elom,theth,emin,ra0,elom2,ixo,tor,
     .              gammin,agrel,bgrel,qgrel,kg
c
      omre=1.0-ra0*eta
      rutfac=omre*enpar2*(omre*enpar2-1.0+elom2)
      if (rutfac .lt. 0)  go to 40
      rut = SQRT (rutfac)
      denomi=1.0/(1.0-omre*enpar2)
      gam1=(elom-rut)*denomi
      gam2=(elom+rut)*denomi
c     arrange roots in increasing order:
      if (gam2 .ge. gam1)  return
      gamd=gam2
      gam2=gam1
      gam1=gamd
      return
c
c     if no physical roots, returns gam1=gam2=-1.23456789e10
c
   40 gam1=-1.23456789e10
      gam2=-1.23456789e10
      return
c
      end


      subroutine coefg (ca,cb,eta)
c
c     calculates coefficients for current drive Green's function:
c      ca=-av((B/Bo)/(zhat*vpar/v));  cb=av(vpar/v), where av is ds/B
c      average.
c
      data pi/3.14159265/
      common /fcncom/ aspect,epar,rl,l,rmin,rmax,rabs,etam,charge,
     .                modl2,asp1,zhat,ap0,alfsq,lm1,lp1,elooro,r3,
     .                eparsq,iheat
c
      if (modl2) 5,30,30
    5 b=1.0-eta-eta*aspect
      if (ABS ((b-aspect)/b) - 1.0e-8) 10,10,20
   10 cb = (rmax * ASIN (SQRT (2.0*aspect/rmax))
     .             + SQRT (asp1*2.0*aspect))/pi
      ca = 0.0
      return
   20 bpasp=b+aspect
      bmasp=b-aspect
      rksq=-alfsq*(1.0-b)/bpasp
      denom=asp1*bpasp
      rutdnm = SQRT (denom)
      rfcoef=2.0*b/bpasp
  130 rdcoef=-rksq/3.
      en=2.0*aspect/bpasp
      en3=en/3.
      rjcoef=en3*bmasp*(1.0+b)/denom
      em=1.0-rksq
      en1=1.0-en
      ifl1=0
      ifl2=0
      ifl3=0
* NAG rf=s21bbE(0.0,em,1.0,ifl1)
* NAG rd=s21bcE(0.0,em,1.0,ifl2)
* NAG rj=s21bdE(0.0,em,1.0,en1,ifl2)
      rf=carlsnf(0.0,em,1.0,ifl1)
      rd=carlsnd(0.0,em,1.0,ifl2)
      rj=carlsnj(0.0,em,1.0,en1,ifl2)
      avfac=rfcoef*rf+rdcoef*rd+rjcoef*rj
      cb=rutdnm*avfac/pi
      ca=-(2.0*rmax/(zhat*pi*rutdnm))*(rmin*rf+en3*bmasp*rj)
      return
   30 ruteta = SQRT (1.0-eta)
      ca=-1.0/(zhat*ruteta)
      cb=ruteta
      return
c
      end
      function ifac (l)
c
      ifac=1
      if (l .lt. 2)  return
      do 15 i=2,l
   15 ifac=ifac*i
      return
c
      end

      subroutine legfn (v, z, p, q, nc, nf)
c
      common /legbl/ vv,zz,pp,qq,acc,ncvg,z1,z2,cvv,svv,nfrig,u,
     .               pisr,pi,a,b,c,gr,r1,r2,r3,r4,zz1,zz2,srz
      complex        vv,zz,pp,qq,z1,z2,u,a,b,c,gr,cvv,svv,zz1,zz2,srz
      complex        v,z,p,q, zzs,cisp,cism,pt,vvp,zzz, csqrtk
c
      data ac   /0.0000001/
c      data qinf /1.79769313486231e+308/
      data qinf /1.79769313486231e+38/
c
      call c311bd
      acc = ac**2
      if (REAL (v)+0.5) 16,17,17
   16 vv  = -v - 1.0
      go to 18
   17 vv  = v
   18 nvv = NINT (REAL (vv))
   30 vvp = vv * pi
      go to (21, 22, 23,2 4, 25), nvv
   21 sv  =  1.0
      go to 26
   23 sv  = -1.0
   26 cv  =  0.0
      go to 28
   22 cv  = -1.0
      go to 27
   24 cv  =  1.0
   27 sv  =  0.0
      go to 28
   25 sv   = SIN (REAL (vvp))
      cv   = COS (REAL (vvp))
   28 eip  = EXP ( AIMAG (vvp))
      eim  = EXP (-AIMAG (vvp))
      shv  = 0.5 * (eip-eim)
      chv  = 0.5 * (eip+eim)
      svv  = sv*chv+u*cv*shv
      cvv  = cv*chv-u*sv*shv
      cisp = cvv+u*svv
      cism = cvv-u*svv
      if (REAL (z))  9, 10, 11
    9 zz   = -z
      n23  = 3
      nff  = -nf
      go to 12
   10 n24  = 4
      zz   = z
      if (AIMAG (zz))  13, 77, 13
   11 zz   = z
      n23  = 2
      nff  = nf
   12 n24  = 2
   13 zzs  = zz**2
      z1   = (1.0-zz) / 2.0
      if (AIMAG (zz))  7, 6, 7
    6 if ( REAL (z1))  4, 70, 5
    4 nfrig = ISIGN (n23, nff)
      go to 8
    5 nfrig = nff
      go to 8
    7 nfrig = SIGN (FLOAT (n24), AIMAG (zz))
    8 z2    = 1.0/zzs
      srz   = csqrtk (zzs-1.0, nfrig, 1)
      zz1   = ( zz+srz) / (2.0*srz)
      zz2   = (-zz+srz) / (2.0*srz)
      if (REAL (z))  1, 2, 2
    1 srz   = -srz
    2 zzz   =  zz
      zz    =  z
      vr =  REAL (vv)**2
      vi = AIMAG (vv)**2
      r1 =  CABS (z1)
      r2 =  CABS (z2)
      r3 =  CABS (zz1)
      r4 =  CABS (zz2)
c
      if (r1-r2)61,61,62
   61 rr = MAX (r3,r4)/0.045-19.5
      if (rr)66,66,65
   65 vri = vr+vi
      if (vri)165,67,165
  165 if (rr**2-vri+vi**2/(2.0*vri))  66, 67, 67
   66 call legv
      if (ncvg)166,97,166
  166 zz = zzz
      call leg1
      go to 90
   67 zz = zzz
      call leg1
      if (ncvg)167,90,167
  167 zz = z
      call legv
      go to 97
   62 if (vr+vi-16.0*r2**2)63,64,64
   64 call legv
      if (ncvg)164,97,164
  164 call legz
      go to 97
   63 call legz
      if (ncvg)163,97,163
  163 call legv
      go to 97
   90 if (REAL (z)) 93,97,97
   93 pt = -2.0/pi*svv*qq
      if (nfrig) 94,95,96
   94 pp =  cisp*pp+pt
      qq = -cism*qq
      go to 97
   95 qq = -cvv*qq-pi/2.0*svv*pp
      pp =  cvv*pp+pt
      go to 97
   96 pp =  cism*pp+pt
      qq = -cisp*qq
      go to 97
c
   97 if (REAL (v)+0.5) 91,92,92
   91 if (CABS (svv) .ne. 0.0)  go to 98
      qq = qinf
      go to 92
   98 qq = (qq*svv-pi*cvv*pp)/svv
   92 p  = pp
      q  = qq
      nc = ncvg
      return
   70 qq = qinf
c
      if ( REAL (z))  173, 173, 74
  173 if (AIMAG (v))   71,  73, 71
   73 go to (71,72,71,74,71),nvv
   71 pp = qq
      go to 82
   72 pp = -1.0
      go to 82
   74 pp = 1.0
      go to 82
c
   77 nfrig = nf
      call legor
   82 ncvg = 0
      go to 92
c
      end

      function erfcb (j,z,zp,zhi)
c
c     normalization*EXP (zhi)*Integral(t**m*EXP (-t^2)) from
c     SQRT (zhi) to infty;  normalization is such that erfcb(,,,0)=1.
c     if j>0, need to dimension z, zp in calling program; then
c     elements of z, zp are integrals for all m .le. 2*(j+1)/2.
c     erfcb is complementary error function; z is above function;
c     zp is derivitive.
c
      complex     zarg, cim
      dimension   z(*), zp(*)
      data        ivar/0/, zero/0.0/
      data
     .   pii,      pin,      dd,       a,        b,        c,        d,
     .   e,        p,        a2,       b2,       c2,       d2,
     .   e2/1.12837916709551  ,        0.56418958354776  ,
     .   0.66666666666667  , 0.24598329392237  , 0.22260371336736  ,
     .   0.12873395368514  , 0.27281029618672  , -.11611455108396  ,
     .   0.43599403657113  , 0.49196658784475  , 0.66781114010207  ,
     .   0.51493581474056  , 1.3640514809336   , -.69668730650377  /
      data        cim/(0.0,1.0)/
c
      zzz = zhi
      if (zzz .le. zero)  go to 2
      if (ivar .gt. 0)  go to 1
      zarg = SQRT (zzz)
      t = 1.0/(1.0+p*zarg)
      zfb = t*(a+t*(a+t*(b+t*(c+t*(d+t*e)))))
      erfcb = zfb
      if (j .eq. 0)  return
      z(1) = zfb+pii*zarg
      if (j .eq. 2)  return
      zp(1) = 0.5*(pii-p*t**2*(a+t*(a2+t*(b2+t*(c2+t*(d2+e2*t))))))/zarg
      if (j .lt. 3)  return
      z(2)=z(1)+pii*zarg*zzz*dd
      if (j .eq. 4)  return
      zp(2)=zp(1)+pii*zarg
      return
c
    1 zar = SQRT (zzz)
      zarg=cim*zar
****  zfb=zprime(zarg,zeta,zdprime)
****  zfb=-zeta*cim*pin
      erfcb=zfb
      if (j .eq. 0)  return
      z(1) = zfb+pii*zar
      zp(1)=zfb
      if (j .le. 2)  return
      z(2)=z(1)+dd*pii*zzz*zar
      zp(2)=z(1)
      return
c
    2 zfb = EXP (zzz)
      erfcb = zfb
      z(2)=zfb
      z(1)= z(2)
      zp(2)=zfb
      zp(1)=zp(2)
      return
c
      end

      subroutine ffp (fu,fupr)
c
c     energy-dependent part of Green's function and derivative.
c     the Green's function is defined as:
c       g=<B_t/R>c^4/4/nuv_t^3B_0)Fu(v/c)h(eta)
c     for kg+2=ig=1, fu=p^2(1-1/(1+bgrel*p^2)**qgrel)/SQRT (1+agrel*p^2)
c      with agrel=1/(5+Z)^2, bgrel=1-agrel, qgrel=1/bgrel.  This form
c      gives Green's function which is correct in nonrelativistic
c      limit; in limit of no trapped particles reproduces leading-order
c      (v^6/c^6) relativistic correction as well as ultrarelativistic
c      limit as given by Fisch.
c     for kg+2=ig=2, fu=gamma1*v^4/c^4, (correct nonrelativistically; correct
c      energy dependence at large energy but wrong coefficient)
c     for kg+2=ig=3, fu=p^4/(1+cgrel*p^2)^(2/3) with
c      p=momentum/mc and cgrel=2/3 (correct to order v^6/c^6 and correct
c      energy dependence at large energy but wrong coefficient)
c
      common /energy/gam,eta,psq,gamsq
      common /rela/enpar2,elom,theth,emin,ra0,elom2,ixo,tor,
     . gammin,agrel,bgrel,qgrel,kg
      data cgrel/.66666666666667/
      if (kg) 10,20,30
c
c --- ig=1 ------------------
c
   10 bfac=1.0+bgrel*psq
      bfacq=bfac**qgrel
      afac=1.0+agrel*psq
      rutaf = SQRT (afac)
      fu=psq*(1.0-1.0/bfacq)/rutaf
      fupr=2.0*gam*(fu*(1.0/psq-.5*agrel/afac)+psq/(rutaf*bfac*bfacq))
      return
c
c --- ig=2 ------------------
c
   20 fu=psq*psq/(gamsq*gam)
      fupr=psq*(1.0+3.0/gamsq)/gamsq
      return
c
c --- ig=3 ------------------
c
   30 denom=1.0+cgrel*psq
      fu=psq*psq/denom**1.5
      fupr=fu*gam*(4.+cgrel*psq)/(psq*denom)
      return
c
      end




      complex function csqrtk (z,n,m)
c
      complex  z, sz
c
      nf=n+5
      sz = SQRT (z)
      s  = ABS (AIMAG (sz))
      go to (21,22),m
   21 go to (9,2,2,9,11,11,2,2,11),nf
   22 go to (2,11,2,2,2,2,2,9,2),nf
    9 sz = CMPLX (0.0,-s)
      go to 2
   11 sz = CMPLX (0.0,s)
    2 csqrtk = sz
      return
c
      end

      subroutine legor
c
      common /legbl/ vv,zz,pp,qq,acc,ncvg,z1,z2,cvv,svv,nfrig,u,
     .               pisr,pi,a,b,c,gr,r1,r2,r3,r4,zz1,zz2,srz
      complex        vv,zz,pp,qq,z1,z2,u,a,b,c,gr,cvv,svv,zz1,zz2,srz
c
      complex        rgam
c
      a=0.5*(1.0+vv)
      b = 0.0
      c=0.5
      gr = rgam (a, b, c)
      c=u*pi*vv*0.5
      a = exp ( c)
      b = exp (-c)
      if (nfrig)  9, 10, 11
    9 qq=pisr*0.5*gr*u*a
      go to 12
   10 qq=u*pisr/4.0*gr*(a-b)
      go to 12
   11 qq=-pisr/2.0*gr*b*u
   12 pp=gr*(a+b)/(2.0*pisr)
      return
c
      end

      subroutine legz
c
      complex        vv,zz,pp,qq,z1,z2,u,a,b,c,gr,cvv,svv,zz1,zz2,srz
      common /legbl/ vv,zz,pp,qq,acc,ncvg,z1,z2,cvv,svv,nfrig,u,
     .               pisr,pi,a,b,c,gr,r1,r2,r3,r4,zz1,zz2,srz
c
      complex        f1, f2, zv, clogok, rgam
c
      gr=clogok(2.0*zz,nfrig,2)
      a=1.5
      b=1.0
      zv = rgam (vv, a, b) * EXP (vv*gr)
      a=vv/2.0+1.0
      b=vv/2.0+0.5
      c=vv+1.5
      accc=acc/100.0
      call hypgm(a,b,c,z2,f1,accc,ncvg)
      f1=f1/(2.0*zz*zv)
      qq=pisr*f1
      if (CABS (cvv) - 0.001)  10, 10, 9
c
c   trdz expects gr = clogok(2.0*zz,nfrig,2) but destroys contents
c
   10 call trdz
      go to 80
    9 a=-vv/2.0
      b=(1.0-vv)/2.0
      c=0.5-vv
      call hypgm(a,b,c,z2,f2,accc,ncv )
      ncvg=ncvg+2*ncv
      f2=f2*zv/(vv+0.5)
      pp=(f1*svv/cvv+f2)/pisr
   80 return
c
      end

      function hfun (eta,hpr)
c
c     evaluates spline set up by setara to find angle-dependent part of
c     Spitzer-Harm function and its derivative.
c     NOTE: h is odd in vparallel.
c           For vparallel < 0 reverse the signv of indicated output.
c
      common /splicm/rk,c,etmax,npd
      common /hfcom/modela,h3fac,etamh
      dimension rk(204),c(204),s(4)
      data left/0/
c
      if (modela-3)  5, 15, 5
    5 if (eta .gt. etmax)  go to 10
* NAG call e02bcE(npd+4,rk,c,eta,left,s,ifail)
      call bsplint (rk,c,npd,4,eta,s,2)
      hpr=s(2)
      hfun=s(1)
      return
   10 hfun = 0.0
      hpr = 0.0
      return
   15 if (eta .gt. etamh)  go to 10
      call geth3(hh,hhpr,eta)
      hfun=h3fac*hh
      hpr=hhpr*h3fac
      return
c
      end

      subroutine leg1
c
      common /legbl/vv,zz,pp,qq,acc,ncvg,z1,z2,cvv,svv,nfrig,u
     .               ,pisr,pi,a,b,c,gr,r1,r2,r3,r4,zz1,zz2,srz
      complex vv,zz,pp,qq,z1,z2,u,a,b,c,gr,cvv,svv,zz1,zz2,srz
c
      complex    sigma, fac
      complex    clogok, psifun
      data gamma1 /0.5772156649/
c
      a=-vv
      b=vv+1.0
      c=(1.0,0.0)
      call hypgm(a,b,c,z1,pp,acc,ncvg)
      sigma=(0.0,0.0)
      fac=(1.0,0.0)
      a = 0.5 * clogok((zz+1.0)/(zz-1.0),nfrig,1) - gamma1 - psifun(b)
      qq=a
      do 17 l=1,50
      el=l
      sigma=sigma+1.0/el
      fac=-fac*(vv+el)*(vv-el+1.0)*z1/(el*el)
      b=(a+sigma)*fac
      qq=qq+b
      if (crit(qq,b,acc))80,17,17
   17 continue
      ncvg=ncvg+2
   80 return
c
      end

      subroutine c311bd
c
      common /legbl/ vv,zz,pp,qq,acc,ncvg,z1,z2,cvv,svv,nfrig,u,
     .               pisr,pi,a,b,c,gr,r1,r2,r3,r4,zz1,zz2,srz
      complex vv,zz,pp,qq,z1,z2,u,a,b,c,gr,cvv,svv,zz1,zz2,srz
      data ibd/0/
c
      if (ibd .ne. 0)  return
      ibd  = 1
      pisr = 1.7724 53850 90552
      pi   = 3.141592653589793
      u    = (0.0,1.0)
      return
c
      end

      subroutine legfn1 (v,z,p,ppr,q,nc,nf)
c
c     calculates ppr, the derivative of p, as well as p and q.
      complex v,z,p,ppr,q,vm1,pm,qm,c1,zsq
      data c1/(1.0,0.0)/
c
      call legfn(v,z,p,q,nc,nf)
      vm1=v-c1
      call legfn(vm1,z,pm,qm,nc1,nf)
      zsq=z*z
      if (zsq .eq. c1)  go to 5
      ppr=(v/(zsq-c1))*(z*p-pm)
      return
    5 ppr=p*v*(v+c1)/2.
      return
c
      end

      subroutine legv
c
      common /legbl/ vv,zz,pp,qq,acc,ncvg,z1,z2,cvv,svv,nfrig,u,
     .               pisr,pi,a,b,c,gr,r1,r2,r3,r4,zz1,zz2,srz
      complex        vv,zz,pp,qq,z1,z2,u,a,b,c,gr,cvv,svv,zz1,zz2,srz
c
      complex        f1, f2, ssrz, clogok, csqrtk, rgam
c
      a=0.5
      c=vv+1.5
      call hypgm(a,a,c,zz1,f1,acc,ncvg)
      call hypgm(a,a,c,zz2,f2,acc,ncv )
      ncvg=ncvg+2*ncv
    6 f1=f1*EXP ((vv+0.5)*clogok(zz +srz,nfrig,3))
      f2=f2*EXP ((vv+0.5)*clogok(zz -srz,-nfrig,3))
      a=1.5
      b=1.0
      ssrz = csqrtk (2.0*srz, nfrig, 2) * rgam (vv, a, b)
      sgn=1.0
      if (AIMAG (zz))  8, 14, 12
    8 sgn=-1.0
      go to 12
   14 sgn=SIGN (1.0,(FLOAT (nfrig)+0.5) * REAL (zz))
   12 pp=(f1+sgn*u*f2)/(pisr*ssrz)
      if (nfrig)11,10,11
   10 qq=0.5*pisr*(f2+sgn*u*f1)/ssrz
      go to 80
   11 qq=pisr*f2/ssrz
   80 return
c
      end

      complex function rgam (z, a, b)
c
      integer iotty
      data    iotty /6/
c
      complex  z, a, b, za, zb, clogam
c
      za=z+a
      zb=z+b
      xa = -REAL (za)
      xb = -REAL (zb)
      ya = AIMAG (za)
      yb = AIMAG (zb)
      if (xa .eq. xb .and. ya .eq. yb)  go to 2
      la=0
      lb=0
      if (ya .ne. 0. .or. xa .ne. ABS (AINT (xa))) la=1
      if (yb .ne. 0. .or. xb .ne. ABS (AINT (xb))) lb=1
      if (la .eq. 1 .and. lb .eq. 1)  go to 1
      rgam = 0.0
      if (la .eq. 1 .and. lb .eq. 0)  return
      write  (iotty, 10)  z, a, b
   10 format (' RGAM ... is not defined or infinite for z = ',2e12.4,
     .          5h a = , 2e12.4, 5h b = , 2e12.4)
      return
c
    2 rgam = 1.0
      return
c
    1 rgam = EXP (clogam(za)-clogam(zb))
      return
c
      end

      subroutine trdz
c
      common /legbl/ vv,zz,pp,qq,acc,ncvg,z1,z2,cvv,svv,nfrig,u,
     .               pisr,pi,a,b,c,gr,r1,r2,r3,r4,zz1,zz2,srz
      complex        vv,zz,pp,qq,z1,z2,u,a,b,c,gr,cvv,svv,zz1,zz2,srz
      complex        rgam, psifun
      complex        ep, ep2, pie, zzs2, aa, bb, cc, fac, sum, add
c
c   trdz  expects gr set = clogok(2.0*zz,nfrig,2)
c
      k  = REAL (vv+1.0)
      ek = k
      ep = vv+0.5-ek
      zzs2 = z2/4.0
      kk = k-1
      sum = 0.0
      if (k)  9, 9, 10
c
   10 b = 0.5
      c = 1.0
      fac = rgam (vv, b, c) * EXP (vv*gr)
      sum = fac
      if (kk)  9, 11, 12
   12 aa = -vv-2.0
      bb = -vv-0.5
      do 2 i=1,kk
        bb = bb+1.0
        aa = aa+2.0
        fac = fac*(1.0+aa)*aa*zzs2/(bb*FLOAT (i))
        sum = sum+fac
        if (crit(sum,fac,acc)) 51, 2, 2
    2 continue
c
      go to 11
c
    9 sum = 0.0
   51 continue
   11 pie = pi*ep
      ep2 = ep/2.0
      a = ek
      b = 0.5-ep
      c = 1.0
      bb = 0.0
      fac = -(1.0-pie**2/3.0) * rgam (a, b, c) * rgam (c, bb, -ep)
     .                        / (pi*EXP ((2.0*ek-vv)*gr))
      a = ek+0.5
      b = 1.0-ep2
      c = ek+1.0+ep2
      add = 2.0 * psifun (a) - psifun (b) - psifun (c) - 2.0 * gr
      gr = fac*add*(1.0+ep2*add)
      sum = sum+gr
      aa = ek-1.5-ep
      bb = -ep
      cc = ek
c
      do 22 i=1,50
        aa = aa+2.0
        bb = bb+1.0
        cc = cc+1.0
        fac = fac*(1.0+aa)*aa*zzs2/(bb*cc)
        add = add+2.0/(a+1.0)+2.0/a-1.0/b-1.0/c
        a = a+2.0
        b = b+1.0
        c = c+1.0
        gr = fac*add*(1.0+ep2*add)
        sum = sum+gr
        if (crit(sum,gr,acc)) 52, 22, 22
   22 continue
c
      ncvg = ncvg+4
   52 pp = sum/pisr
      return
c
      end

      subroutine hypgm (a,b,c,z,h,acc,ncvg)
c
      complex a,b,c,z,h,aa,bb,cc,zz,add,dd,hh
c
      ncvg=0
      zz=z
      hh=(1.0,0.0)
      aa=a
      bb=b
      cc=c
      add=(1.0,0.0)
      dd=(1.0,0.0)
      do 21 i=1,50
      add=add*zz*aa/dd*bb/cc
      hh=hh+add
      if (crit(hh,add,acc))3,4,4
    4 aa=aa+1.0
      bb=bb+1.0
      cc=cc+1.0
      dd=dd+1.0
   21 continue
      ncvg=1
    3 h=hh
      return
c
      end

      complex function clogok (z,n,m)
c
      complex z,sz
      data pi/3.1415926535898/
      sz =  log (z)
      s  = REAL (sz)
      nf=n+5
      go to (21,22,23) ,m
   21 go to (2,2,2,11,10,9,2,2,2),nf
   22 go to (2,11,2,2,2,2,2,9,2),nf
   23 go to (2,11,2,2,2,2,2,9,2),nf
    9 sz = CMPLX (s,-pi)
      go to 2
   10 sz = CMPLX (s,0.0)
      go to 2
   11 sz = CMPLX (s,pi)
    2 clogok=sz
   80 return
c
      end


      function crit (sum,del,accs)
c
      complex sum,del
c
      crit =           real (del)**2 + AIMAG (del)**2
     .       - accs * (REAL (sum)**2 + AIMAG (sum)**2)
      return
c
      end




      complex function psifun (z)
c
      complex    z, u, v, h, r
      dimension  b(6)
c
      integer iotty
      data    iotty /6/
c
      data pi /  3.14159 26535 89793/
      data b  / +8.33333 33333 33333e-2, -8.33333 33333 33333e-3,
     .          +3.96825 39682 53968e-3, -4.16666 66666 66667e-3,
     .          +7.57575 75757 57576e-3, -2.10927 96092 79609e-2 /
c
      u = z
      x = REAL (u)
      a =  ABS (x)
      if (AIMAG (u) .eq. 0.0 .and. -a .eq. AINT (x))  go to 4
      if (x .lt. 0.0)  u = -u
      v = u
      h = 0.0
      if (a .ge. 15.0)  go to 3
      n = 14 - INT (a)
      h = 1.0 / v
      if (n .eq. 0)  go to 2
      do 1 i = 1,n
      v=v+1.0
    1 h=h+1.0/v
    2 v=v+1.0
    3 r = 1.0 / v**2
      psifun = LOG (v) - 0.5/v
     .  - r * (b(1)+r*(b(2)+r*(b(3)+r*(b(4)+r*(b(5)+r*(b(6)+r*b(1)))))))
     .  - h
      if (x .ge. 0.0)  return
      h      = pi * u
      psifun = psifun  +  1.0 / u  +  pi * COS (h) / SIN (h)
      return
c
    4 write  (iotty, 100) x
  100 format (' PSIFUN ... argument is non-positive integer = ',f20.2)
      psifun = 0.0
      return
c
      end

      complex function clogam (z)
c
      complex   z, v, h, r
      dimension b(10)
c
      integer iotty
      data    iotty /6/
c
      data pi /3.14159 26535 89793/
      data b  /+8.33333 33333 3333e-2, -2.77777 77777 7778e-3,
     .         +7.93650 79365 0794e-4, -5.95238 09523 8095e-4,
     .         +8.41750 84175 0842e-4, -1.91752 69175 2692e-3,
     .         +6.41025 64102 5641e-3, -2.95506 53594 7712e-2,
     .         +1.79644 37236 8831e-1, -1.39243 22169 0590e+0/
c
      x =  REAL (z)
      t = AIMAG (z)
      if (-ABS (x) .eq. AINT (x) .and. t .eq. 0.0)  go to 5
      f = ABS (t)
      v = CMPLX (x,f)
      if (x .lt. 0.0)  v = 1.0 - v
      h = (0.0, 0.0)
      c = REAL (v)
      if (c .ge. 7.0)  go to 3
      n = 6 - INT (c)
      h=v
      d = AIMAG (v)
      a = ATAN2 (d,c)
      if (n .eq. 0)  go to 2
      do 1 i = 1,n
      c=c+1.0
      v = CMPLX (c,d)
      h=h*v
    1 a = a + ATAN2 (d,c)
    2 h = CMPLX (0.5*LOG (REAL (h)**2 + AIMAG (h)**2),a)
      v = v + 1.0
    3 r = 1.0 / v**2
      clogam=0.918938533204673+(v-0.5)*LOG (v)-v+(b(1)+r*(b(2)+r*(b(3)+
     . r*(b(4)+r*(b(5)+r*(b(6)+r*(b(7)+r*(b(8)+r*(b(9)+r*b(10))))))))))
     . /v-h
      if (x .ge. 0.0)  go to 4
c
      a = AINT (x) - 1.0
      c=pi*(x-a)
      d=pi*f
      e = EXP (-2.0*d)
      f = SIN (c)
      e=d+0.5*LOG (e*f**2+0.25*(1.0-e)**2)
      f = ATAN2 (COS (c)*TANH (d),f)-a*pi
      clogam = 1.144729885849400 - CMPLX (e,f)-clogam
c
    4 if (t .lt. 0.0)  clogam = conjg(clogam)
      return
c
    5 write (iotty, 100) x
      clogam = (0.0, 0.0)
      return
  100 format (' clogam ... argument is non-positive integer = ',f20.2)
c
      end

