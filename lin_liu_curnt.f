      subroutine TorGA_curgap(rjpd,rjpd0,ratjpd,denom,eps,npara,nperp  
     &,omode,cefldx,cefldy,cefldz,tebulk,thtc,thetap,yy,lh,zeffin,model
     &,tol,ngauss,ig)
c-----This subroutine was transformed from Fortran 90 to Fortran 77
c      use precision_mod
c      use real_mod
      implicit none

      REAL*8 rjpd,rjpd0,ratjpd,denom
      REAL*8 eps,npara,nperp,omode,tebulk,thtc,thetap,yy,tol,zeffin
      COMPLEX*16 cefldx,cefldy,cefldz
      INTEGER lh,model,ngauss,ig
      INCLUDE 'globcd.h'
      INCLUDE 'globcd1.h'
      INCLUDE 'globcd2.h'
      REAL*8 TorGa_zgauleg
      EXTERNAL TorGa_zgauleg
      EXTERNAL TorGa_funxjs,TorGa_funxjz,TorGa_alpha
!------------------------------------------------------------------------
!     CURGAP(version 1.0) Y.R. Lin-Liu, 08/16/04
!     This routine is an improved version and the replacement of CURGAC.
!     The code calculates a normalized ECCD efficiency using the Green's
!     function formulation[1,2]. The improvement is the use of the exact
!     polarization-dependent rf quasi-linear diffusion operator[3]in 
!     evaluating the current drive efficiency.
!
!     References:
!     1. Y.R. Lin-Liu, V.S. Chan, and R. Prater,"Electron cyclotron
!        current drive efficiency in general tokamak gemetry," 
!        Phys. Plasmas 10, 4064 (2003).
!     2. R.H. Cohen, "Effect of trapped electrons on current drive,"
!        Phys. Fluids 30, 2442 (1987); 31, 421 (1988).
!     3. R.W. Harvey and M.G. McCoy, "The CQL3D Fokker-Plank code,"
!        GA report GA-A20978 (1992).
!
!     OUTPUTS:
!     rjpd = Cohen's normalized current drive efficiency (j/Pd),
!            where j is the flux-surface-average of the parallel driven 
!            current density (in units of e(*n_e)*Vt) and Pd is the 
!            absorbed power density ( in units of nu*n_e*m*Vt^2).
!            Here n_e = electron density
!                 Vt  = SQRT(2*T/m)
!                 nu  = e^4*n_e*Log(Lambda)/(4*pi*eps0^2*m^2*Vt^3)
!     rjpd0 = the equivalent rjpd in the absence of trapped electrons
!     ratjpd = rjpd/rjpd0
!     denom  = the denomintor in evaluation of rjpd
!
!     INPUTS:
!     eps = inverse aspect ratio of the flux surface of interest
!     npara = parallel index of refraction
!     nperp = Real part of perpendicular index of refraction
!     omode = +1.0 for O-mode and -1.0 for x-mode  [YuP: Not used?]
!     cefldx = the x-component of the wave electric field (COMPLEX)
!     cefldy = the y-component                            (COMPLEX)
!     cefldy = the z-component                            (COMPLEX)
!     tebulk = Electron temperature in keV's
!     thtc   = an obsolete variable
!     thetap = poloidal angle at which power is absorbed in radians
!              0: outborad; pi: inboard
!     lh = the cyclotron harmonic number
!     yy = lh*omega_c/omega (y in Refs [1] and [2])
!     zeffin = effective ion charge
!     model < 5  gives rjpd in CURGAC
!           = 5  gives rjpd using the exact polarization-dependent
!                rf diffusion operator
!           > 5  gives rjpd using the polarization-dependent rf diffusion
!                operator with but small gyro-radius expansion 
!     tol = The relative error tolerence in numerical integration is set
!           to be MAX (tol, 1.0E-6).
!     ngauss = number of points in the Gaussian quadrature (ngauss = 64
!              is recommended.)
!     ig     = 1 : the relativistic Green's function
!            = 0 : the non-relativistic approximation
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!     Comments from older versions of CURGAC
!     revisions (02/13/99):
!     improve numerical integeration of ft and hcap
!     introduce ismall for flaging small eps
!     introduce new routine "getftrap"
!     combine interation and integrant routines
!     revisions (02/10/99):
!     introduce nghalf and use TorGa_zgauleg
!     remove the bug in cxi2
!
!     revisions (02/01/99):
!     set ng=min(ngauss,64)
!     remove ptfctr option
!     if model .eq. 4 will calculate the collsionality correction of CD
!     add new variables: rjpdc, nustar, ccf
!     add a new global variable: eta
!
!------------------------------------------------------------------------
      REAL*8 rteps,ss
      REAL*8 etmax,epst1,epst2
      REAL*8 expt,expb,exp1,exp2,rint,rint0
      INTEGER ng,nghalf,ifail
      INTEGER jromb
      INTEGER ismall
!
      pi=dacos(-1.d0)
!
      if (dabs(npara) .lt. small)then   ! don't bother 
         rjpd=0.d0
         rjpd0=0.d0
         ratjpd=1.d0
         denom=1.d0
         return
      endif
!
!     set up the global variables for the current-drive package
!     numerical control parameters:
      tolval=dmax1(tol,1.d-6)
      nghalf=MIN(ngauss,64)/2
      ng=nghalf*2
      modelv=model
      igv=ig
!
!     variables of wave characteristics:
      yval=yy
      nharm=lh
      omodev=omode  ! [YuP: Not used?]
      enz=npara
      enzsq=enz*enz
      sgnnz=dsign(1.0d0,npara)
      enperp=nperp
      cez=cefldz
      ceplus=cefldx+(0.d0, 1.d0)*cefldy
      ceminus=cefldx-(0.d0, 1.d0)*cefldy
!
!     geometric variables:
      call TorGa_ceqmdl(eps,thetap)
      rteps=sqrt(abs(1.d0-hav))
      if (rteps .lt. small)then
         ismall=1
         cxi2=-0.125d0
         cxi4=pi**2/8.d0-1.d0
         ft=1.46d0*rteps
         fc=1.d0-ft
      else
         ismall=0
         cxi2=-0.25d0*(hsqav-hav**2)/(1.d0-hav)**2
         cxi4=chrtav**2/(1.d0-hav)-1.d0
         call TorGa_getftrap
      endif
!
!     kinetic variables:
      tau=tebulk/tmass     ! tmass given in globcd.h
      zeff=zeffin
      zrat=(zeff+1.d0)/fc
      etcutoff=20.d0    ! upper bound for (gamma-1)/tau integration
!
!     setup integration limits
      call TorGa_getlims(etmax,epst1,epst2)
      if (ismall .eq. 1)then          ! del with samll eps case
         if (etmax .gt. 0.d0)then
            ifail=0
            expb=dEXP(-etmax)
            expb=1.d0
            denom=tau*TorGa_zgauleg(TorGa_alpha,expb,expt,nghalf,2,
     &                              ifail)
            rint0=TorGa_zgauleg(TorGa_funxjz,expb,expt,nghalf,2,ifail)
            rjpd0=-0.5d0*hav*rint0/denom
            rjpd=rjpd0
            ratjpd=1.d0
         else
            denom=0.d0
            rjpd0=0.d0
            rjpd=0.d0
            ratjpd=1.d0
         endif
         return
      endif
!
      if (etmax .gt. 0.d0)then
!        evaluate the normalized absorbed power
!        call pdnorm(denom,tau,yval,enzsq,nharm,gammin,gammax)
         expt=1.d0
         expb=dexp(-etmax)
         denom=tau*TorGa_zgauleg(TorGa_alpha,expb,expt,nghalf,2,ifail)
         rint0=TorGa_zgauleg(TorGa_funxjz,expb,expt,nghalf,2,ifail)
         if (epst1 .gt. 0.d0 .and. epst2 .lt. etmax)then
            expb=dexp(-etmax)
            exp1=dexp(-epst1)
            exp2=dexp(-epst2)
            ifail=0
            rint=TorGa_zgauleg(TorGa_funxjs,expb,exp2,nghalf,1,ifail)
            rint=rint+TorGa_zgauleg(TorGa_funxjs,exp1,expt,nghalf,1,
     &                              ifail)
         else
            ifail=0
            if(epst2.le.0.d0)expb=dexp(-etmax)
            if(epst2.gt.0.d0.and.epst1.le.0.d0)expb=
     &                                         dexp(-MIN(epst2,etmax))
            if(epst2.gt.0.d0.and.epst1.gt.0.d0)expb=
     &                                         dexp(-MIN(epst1,etmax))
            rint=TorGa_zgauleg(TorGa_funxjs,expb,expt,nghalf,2,ifail)
!%PR         write (6,*)rint0,rint
         endif
         rjpd=-0.5d0*hav*rint/denom
         rjpd0=-0.5d0*hav*rint0/denom
         ratjpd=rint/rint0
      else
         denom=0.d0
         rjpd=0.d0
         rjpd0=0.d0
         ratjpd=1.d0
      endif
!
      return
!
      end
      real*8 FUNCTION TorGa_alpha(x)
!
c      USE precision_mod
c      USE real_mod
      IMPLICIT NONE
      REAL*8 x
      REAL*8 TorGa_bessj
      INCLUDE 'globcd.h'
      INCLUDE 'globcd1.h'
!-------------------------------------------------------------------------
      REAL*8 gamma,usq,u,uz,uxsq,uperp
      REAL*8 aperp,xjz,xjp,xjm
      COMPLEX*8 cuel
      INTEGER np,nm
!
      gamma=-log(x)*tau+gammin 
      usq=dABS(gamma**2-1.d0)
      u=dSQRT(usq)
      uz=(gamma-yval)/enz      
      uxsq=usq-uz**2
      uperp=dSQRT(dABS(uxsq))
!
      aperp=nharm*enperp*uperp/yval
!
      IF (modelv == 5) THEN
         xjz=TorGa_bessj(nharm,aperp)
         np=nharm+1
         xjp=TorGa_bessj(np,aperp)
         nm=nharm-1
         xjm=TorGa_bessj(nm,aperp)
         cuel=cez*uz*xjz+0.5d0*uperp*(ceminus*xjm+ceplus*xjp)
         TorGa_alpha=ABS(cuel)**2
      ELSE IF (modelv < 5)THEN
         TorGa_alpha=uxsq**nharm
      ELSE
         xjm=1.d0-0.25d0*aperp**2/nharm
         xjz=0.5d0*aperp/nharm
         xjp=0.25d0*aperp**2/(nharm*(nharm+1))
         cuel=cez*uz*xjz+0.5d0*uperp*(ceminus*xjm+ceplus*xjp)
         TorGa_alpha=uxsq**(nharm-1)*ABS(cuel)**2
      ENDIF
!
      RETURN
      END FUNCTION TorGa_alpha

      real*8 function TorGa_funxjz(x)

c      use precision_mod
c      use real_mod
      implicit none
      REAL*8 x
      REAL*8 TorGa_alpha
      include 'globcd.h'
      include 'globcd2.h'
! N Bertelli 15 July 2014
      include 'param.i'
      include 'one_nml.i'
!-------------------------------------------------------------------------
      REAL*8 alpha
      REAL*8 gamma,usq,u,uz,uxsq
      REAL*8 fcap,fprime
      REAL*8 hksi,hcap
!N Bertelli 14 July 2014 from emp
      REAL*8 tebulk
!
      alpha=TorGa_alpha(x)
!
      gamma=-log(x)*tau+gammin 
      usq=abs(gamma**2-1.d0)
      u=dsqrt(usq)
      uz=(gamma-yval)/enz      
      uxsq=usq-uz**2
!
      hksi=hloc/hsqav
      hcap=uz/u*hksi
!
      zrtmp=zeff+1.d0
      fctmp=1.d0
! N Bertelli 14 July 2014 from emp
      if (ieffic_mom_cons == 0) then
      call TorGa_getfcap(u,fcap,fprime)
      elseif (ieffic_mom_cons == 1) then
      tebulk = tau*tmass
      call green_func_emp(igv, tebulk, zeff, fc, u, gamma, fcap, fprime)
      endif
!
      TorGa_funxjz=alpha*(fprime*hcap+(enz-gamma*uz/usq)*fcap*hksi/u)
!
      return
      end

      real*8 function TorGa_funxjs(x)

c      use precision_mod
c      use real_mod
      implicit none

      REAL*8 x
      REAL*8 TorGa_alpha
      include 'globcd.h'

! N Bertelli 15 July 2014
      include 'param.i'
      include 'one_nml.i'
!-------------------------------------------------------------------------
      REAL*8 alpha
      REAL*8 gamma,usq,u,uz,uxsq,lambda
      REAL*8 sgnup,betayr,hcap,hprime,hksi,fcap,fprime
!N Bertelli 14 July 2014 from emp
      REAL*8 tebulk
!
      alpha=TorGa_alpha(x)
!
      gamma=-dlog(x)*tau+gammin 
      usq=dabs(gamma**2-1.d0)
      u=dsqrt(usq)
      uz=(gamma-yval)/enz      
      uxsq=usq-uz**2
      lambda=uxsq/usq/hloc
!
      sgnup=dsign(1.0d0,uz)
!%LL      beta=2./hloc*(gamma-yval)/usq*(gamma*(gamma-yval)/(usq*enzsq)-1.)
      betayr=(enz-gamma*uz/usq)
      call TorGa_gethcap(lambda,hcap,hprime)
      hksi=-2.d0/hloc*(uz/u)*hprime
! N Bertelli 14 July 2014 from emp
      if (ieffic_mom_cons == 0) then
      call TorGa_getfcap(u,fcap,fprime)
      elseif (ieffic_mom_cons == 1) then
      tebulk = tau*tmass
      call green_func_emp(igv, tebulk, zeff, fc, u, gamma, fcap, fprime)
      endif
!
!%LL      TorGa_funcds=sgnup*alpha*(fprime*hcap+beta*fcap*hprime)
      TorGa_funxjs=sgnup*alpha*(fprime*hcap+betayr*fcap*hksi/u)
!
      return
      end

      real*8 FUNCTION TorGA_bessj0(x)
c      USE precision_mod
      IMPLICIT NONE
!
      REAL*8 x
      REAL*8 ax,xx,z
      REAL*8 p1,p2,p3,p4,p5
      REAL*8 q1,q2,q3,q4,q5
      REAL*8 r1,r2,r3,r4,r5,r6
      REAL*8 s1,s2,s3,s4,s5,s6
      REAL*8 y
      SAVE p1,p2,p3,p4,p5
      SAVE q1,q2,q3,q4,q5
      SAVE r1,r2,r3,r4,r5,r6
      SAVE s1,s2,s3,s4,s5,s6
!
!234567890123456789012345678901234567890123456789012345678901234567890123
!-----------------------------------------------------------------------&
      DATA p1,p2,p3,p4,p5/1.E0,-.1098628627E-2                    
     &,.2734510407E-4,-.2073370639E-5,.2093887211E-6/
      DATA q1,q2,q3,q4,q5/-.1562499995E-1,.1430488765E-3          
     &,-.6911147651E-5,.7621095161E-6,-.934945152E-7/         
      DATA r1,r2,r3,r4,r5,r6/57568490574.E0,-13362590354.E0       
     &,651619640.7E0,-11214424.18E0,77392.33017E0             
     &,-184.9052456E0/
      DATA s1,s2,s3,s4,s5,s6/57568490411.E0,1029532985.E0         
     &,9494680.718E0,59272.64853E0,267.8532712E0,1.E0/
!-----------------------------------------------------------------------&
!
      IF (dABS(x).lt. 8.d0)THEN
        y=x**2
        TorGA_bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))              
     &               /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      ELSE
        ax=dABS(x)
        z=8.d0/ax
        y=z**2
        xx=ax-.785398164d0
        TorGA_bessj0=dsqrt(.636619772d0/ax)*(dcos(xx)*                   
     &              (p1+y*(p2+y*(p3+y*(p4+y*p5))))                      
     &               -z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      ENDIF
!
      RETURN
      END FUNCTION TorGA_bessj0

      real*8 FUNCTION TorGA_bessj1(x)
c      USE precision_mod
c      USE real_mod
      IMPLICIT NONE
!
      REAL*8 x
      REAL*8 ax,xx,z
      REAL*8 p1,p2,p3,p4,p5
      REAL*8 q1,q2,q3,q4,q5
      REAL*8 r1,r2,r3,r4,r5,r6
      REAL*8 s1,s2,s3,s4,s5,s6
      REAL*8 y
      SAVE p1,p2,p3,p4,p5
      SAVE q1,q2,q3,q4,q5
      SAVE r1,r2,r3,r4,r5,r6
      SAVE s1,s2,s3,s4,s5,s6
!
!234567890123456789012345678901234567890123456789012345678901234567890123
!-----------------------------------------------------------------------&
      DATA r1,r2,r3,r4,r5,r6/72362614232.E0,-7895059235.E0        
     &,242396853.1E0,-2972611.439E0,15704.48260E0              
     &,-30.16036606E0/
      DATA s1,s2,s3,s4,s5,s6/144725228442.E0,2300535178.E0        
     &,18583304.74E0,99447.43394E0,376.9991397E0,1.E0/
      DATA p1,p2,p3,p4,p5/1.E0,.183105E-2,-.3516396496E-4      
     &,.2457520174E-5,-.240337019E-6/
      DATA q1,q2,q3,q4,q5/.04687499995E0,-.2002690873E-3          
     &,.8449199096E-5,-.88228987E-6,.105787412E-6/
!-----------------------------------------------------------------------&
!
      IF (dABS(x).lt.8.d0) THEN
        y=x**2
        TorGA_bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/           
     &           (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      ELSE
        ax=dABS(x)
        z=8.d0/ax
        y=z**2
        xx=ax-2.356194491d0
        TorGA_bessj1=dsqrt(.636619772d0/ax)*                            
     &        (dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))                   
     &        -z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*dsign(1.d0,x)
      ENDIF
!
      RETURN
      END FUNCTION TorGA_bessj1

      real*8 FUNCTION TorGA_bessj(n,x)
c      USE precision_mod
c      USE real_mod
      IMPLICIT NONE
!
      INTEGER n,IACC
      REAL*8 x,BIGNO,BIGNI
      PARAMETER (IACC=40,BIGNO=1.d10,BIGNI=1.d-10)
!     USES bessj0,bessj1
      REAL*8 TorGA_bessj0,TorGA_bessj1
      EXTERNAL TorGA_bessj0,TorGA_bessj1
!-----------------------------------------------------------------------&
!234567890123456789012345678901234567890123456789012345678901234567890123
!-----------------------------------------------------------------------&
      INTEGER j,jsum,m
      REAL*8 ax,bj,bjm,bjp,sum,tox
!
      IF (n.lt.0) THEN
         WRITE(6,"(A)")'bad argument n in bessj'
         STOP
      ENDIF
!
      IF (n.eq.0) THEN
         TorGA_bessj=TorGA_bessj0(x)
         RETURN
      ENDIF
!
      IF (n.eq.1) THEN
         TorGA_bessj=TorGA_bessj1(x)
         RETURN
      ENDIF
!
      ax=ABS(x)
      IF (ax.eq.0.d0) THEN
        TorGA_bessj=0.d0
      ELSE IF(ax.gt.FLOAT(n)) THEN
        tox=2.d0/ax
        bjm=TorGA_bessj0(ax)
        bj=TorGA_bessj1(ax)
        DO j=1,n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
        ENDDO
        TorGA_bessj=bj
      ELSE
        tox=2.d0/ax
        m=2*((n+INT(dSQRT(dFLOAT(IACC*n))))/2)
        TorGA_bessj=0.d0
        jsum=0
        sum=0.d0
        bjp=0.d0
        bj=1.d0
        DO j=m,1,-1
          bjm=j*tox*bj-bjp
          bjp=bj
          bj=bjm
          IF (dABS(bj).gt.BIGNO) THEN
            bj=bj*BIGNI
            bjp=bjp*BIGNI
            torGA_bessj=TorGA_bessj*BIGNI
            sum=sum*BIGNI
          ENDIF
          IF (jsum.ne.0)sum=sum+bj
          jsum=1-jsum
          IF (j.eq.n)TorGA_bessj=bjp
        ENDDO
!
        sum=2.d0*sum-bj
        TorGA_bessj=TorGA_bessj/sum
      ENDIF

      IF (x.lt.0.d0.and.mod(n,2).eq.1)TorGA_bessj=-TorGA_bessj
!
      RETURN
      END FUNCTION TorGA_bessj

c------------------------------------------------------------------
!     function TorGa_zgauleg
!     uses:
!     use precision_mod
!     use real_mod
!     TorGa_mgauleg +
!     TGaLib_stop   +
!
!     subroutine TorGa_mgauleg
!     uses:
!     use precision_mod
!     use real_mod
!
!     subroutine TGaLib_stop
!     uses:
!     use precision_mod
!     use real_mod
!     GETARG
!     INDEX
!     TGaLib_ishell +
!     EXIT
!
!     integer function TGaLib_ishell 
!     uses:
!     use precision_mod
!     use real_mod
!     SYSTEM
!
!     subroutine TorGa_ceqmdl(eps,thetap)
!     uses:
!     use precision_mod
!     use real_mod
!
!     subroutine TorGA_curga
!     uses:
!     use precision_mod
!     use real_mod
!     INCLUDE 'globcd.h'
!     INCLUDE 'globcd1.h'
!     INCLUDE 'globcd2.h'
!     TorGa_zgauleg +
!     TorGa_funxjs  +
!     TorGa_funxjz  +
!     TorGa_alpha  +
!     TorGa_ceqmdl +
!     TorGa_getftrap +
!     TorGa_getlims  +
!
!     FUNCTION TorGa_alpha
!     uses:
!     USE precision_mod
!     USE real_mod
!     INCLUDE 'globcd.h'
!     INCLUDE 'globcd1.h'
!     TorGa_bessj     +
!
!     function TorGa_funxjz
!     uses:
!     use precision_mod
!     use real_mod
!     TorGa_getfcap +
!     include 'globcd.h'
!     include 'globcd2.h'
!     TorGa_alpha +
!
!     function TorGa_funxjs
!     uses:
!     use precision_mod
!     use real_mod
!     include 'globcd.h'
!     TorGa_gethcap +
!     TorGa_alpha +
!
!
!     FUNCTION TorGA_bessj(n,x)
!     uses:
!     USE precision_mod
!     USE real_mod
!     TorGA_bessj0 +
!     TorGA_bessj1 +
!
!
!     FUNCTION TorGA_bessj1(x)
!     uses:
!     USE precision_mod
!     USE real_mod
!
!     FUNCTION TorGA_bessj0(x)
!     uses:
!     USE precision_mod
!     USE real_mod
!
!     subroutine TorGa_getftrap
!     uses:
!     use precision_mod
!     use real_mod
!     include 'globcd.h'
!     include 'globcd2.h'
!     TorGa_mqromb1 +
!     TorGa_qftint  +
!
!     subroutine TorGa_getlims
!     uses:
!     use precision_mod
!     use real_mod
!     TorGa_gamsrc +
!
!     subroutine TorGa_getfcap
!     uses:
!     use precision_mod
!     use real_mod
!     include 'globcd.h'
!     include 'globcd2.h'
!     TorGa_mqromb1 +
!     TorGa_fcapint +
!
!     subroutine TorGa_gethcap(z,hcap,hprime)
!     uses:
!     use precision_mod
!     use real_mod
!     TorGa_mqromb1 +
!     TorGa_hcapint +
!
!     SUBROUTINE TorGa_mqromb1
!     uses: 
!     use precision_mod
!     use real_mod
!     TorGa_trapzd +
!     TorGa_polint +
!
!     SUBROUTINE TorGa_trapzd
!     uses:
!     use precision_mod
!     use real_mod
!
!     SUBROUTINE TorGa_polint
!     uses:
!     use precision_mod
!     use real_mod
!
!     function TorGa_qftint
!     uses:
!     use real_mod
!     implicit none
!
!     subroutine TorGa_gamsrc
!     uses:
!     use precision_mod
!     use real_mod
!     include 'globcd.h'
!
!     function TorGa_fcapint
!     uses:
!     use precision_mod
!     use real_mod
!
!      function TorGa_hcapint(s)
!     uses:
!     use precision_mod
!     use real_mod
!     include 'globcd.h'
!------------------------------------------
!     external  
!     GETARG,
!     INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
!     EXIT
!     TGaLib_ishell is in portlib.f90
!--------------------------------------------------------------------------------  

      real*8 function TorGa_zgauleg (func, a, b, n, multi, ifail)

c      use precision_mod
c      use real_mod
      implicit none

      REAL*8 func,a,b
      integer n,multi,ifail
      external func
!
!     returns (n*multi)-point gauss-legendre quadrature of 'func'
!     from a to b.  (n*multi) should be less than nmax =64
!     multi = 1 or 2
! ----------------------------------------------------------------------
!
      integer nmax,nmaxh
      parameter (nmax=64,nmaxh=nmax/2)
      integer stderr
      parameter (stderr=6)
      integer nh,nn,ns,j
      REAL*8 x(nmax),w(nmax),xx(nmax),y(nmax),ss,xm,xr
      REAL*8 xh(nmaxh),wh(nmaxh)
      save ns,x,w,xh,wh
      data ns/-1/
!
      nh=n
      nn=2*nh
      if (nh .gt. nmaxh) then
         nh=nmaxh
         nn=nmax
         ifail=1
         write (stderr,'(a)')'WARNING: TorGa_zgauleg: (n) set to 32 '
      endif
      if (ns .ne. nh) then
         ns=nh
         call TorGa_mgauleg(xh,wh,nh)
         call TorGa_mgauleg(x,w,nn)
      endif
!
      xm=0.5d0*(b+a)
      xr=0.5d0*(b-a)
      if (multi .eq. 1)then
         do j=1,nh
            xx(j)=xm+xr*xh(j)
            y(j)=func(xx(j))
         enddo
         ss=0.d0
         do j=1,nh
            ss=ss+wh(j)*y(j)
         enddo
         TorGa_zgauleg=ss*xr
      else if (multi .eq. 2)then
         do j=1,nn
            xx(j)=xm+xr*x(j)
            y(j)=func(xx(j))
         enddo
         ss=0.d0
         do j=1,nn
            ss=ss+w(j)*y(j)
         enddo
         TorGa_zgauleg=ss*xr  
      else
         write (stderr,'(a)')'TorGa_zgauleg:  dimensional error'
         call TGaLib_stop('TorGa_zgauleg', 0)
      endif
      return
!
      end


!--------------------------------------------------------------------------
!     Gauss Legendre package
!--------------------------------------------------------------------------
      subroutine TorGa_mgauleg (x, w, n)
c      use precision_mod
c      use real_mod
      implicit none
      integer n
      REAL*8 x(n),w(n)
!     returns the abscissas and weights for n-point
!     gauss-legendre integration
! ----------------------------------------------------------------------
      REAL*8 eps
      parameter (eps=3.d-14)
      integer i,j,m
      REAL*8 pi,z,p1,p2,p3,pp,z1
!
      pi=dacos(-1.d0)
!     pi = 3.1415926536E0d0

      m=(n+1)/2
      do i=1,m
         z=dcos(pi*(i-0.25d0)/(n+0.5d0))
 1       continue
         p1=1.d0
         p2=0.d0
         do j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
         enddo
         pp=n*(z*p1-p2)/(z*z-1.d0)
         z1=z
         z=z1-p1/pp
         if (dabs(z-z1) .gt. eps)  go to 1
         x(i)=-z
         x(n+1-i)=z
         w(i)=2.d0/((1.d0-z*z)*pp*pp)
         w(n+1-i)=w(i)
      enddo
!
      return
      end
!


      subroutine TGaLib_stop (message, status)
!
! --- print specified message and then exit with the specified status --
!

c      use precision_mod
c      use real_mod
      implicit none
!
      integer   TGaLib_ishell, status, position
      character message*(*), program_name*128, command*267
      external  GETARG
!
      if (status .eq. 0) then
        write (*, '(/ '' TGaLib_stop in '', a /)')
     &  message
      else
        PAUSE 'Fix getarg call'
cBH080103, Fortran 2003 replacement        call GETARG (0, program_name)
c        call GET_COMMAND_ARGUMENT (0, program_name)
        position = INDEX (program_name, 'onetwo')
        if (position .ne. 0) then
          position = INDEX (program_name(position:), '/')
          if (position .eq. 0) then
            command = 'rm -f bpltfil eqpltfil ; '              //       
     &                'touch outone qikone runlog '            //       
     &                      'trpltfil trpltout trpltout.nc ; ' //       
     &                'mv -f outone _outone ; '                   //    
     &                'mv -f qikone _qikone ; '                   //    
     &                'mv -f runlog _runlog ; '                   //    
     &                'mv -f trpltfil _trpltfil ; '               //    
     &                'mv -f trpltout _trpltout ; '               //    
     &                'mv -f trpltout.nc _trpltout.nc ; '         //    
     &                'echo "\n OUTPUT IS INVALID"'
            if (TGaLib_ishell (command) .ne. 0) then
              stop 'Last-chance abort in subroutine TGaLib_stop'
            end if
          end if
        end if
        write (*, '(/ '' TGaLib_stop in '',
     & a / '' EXIT status is'', i4 /)')
     &  message, status
      end if
      call EXIT (status)
!
      end

      integer function TGaLib_ishell (command)
!
! --- replacement for this function in UNICOS --------------------------
!

c      use precision_mod
c      use real_mod
      implicit none
!
      integer       SYSTEM
      character*(*) command
!
      if (SYSTEM (command) .eq. 0) then
        TGaLib_ishell =  0
      else
        TGaLib_ishell = -1
      end if
      return
!
      end
      subroutine TorGa_ceqmdl(eps,thetap)

c      use precision_mod
c      use real_mod
      implicit none

      REAL*8 eps,thetap
      include 'globcd.h'
!     return geometrical quantities of circular model equilibrium with
!     given eps and thetap
!--------------------------------------------------------------------------
!     real pi
!
      pi=dacos(-1.d0)
!     pi = 3.1415926536d0

      hloc=(1.d0-eps)/(1d0+eps*dcos(thetap))
      href=(1d0-eps)/(1d0+eps)
      hav=1d0-eps
      hsqav=(1d0-eps)**2/dsqrt(1d0-eps**2)
      chrtav=(dsqrt(2d0*eps*(1d0-eps))+(1d0+eps)*   
     &        dasin(dsqrt(2d0*eps/(1d0+eps))))/pi
!
!     cxi2=-0.25*(hsqav-hav**2)/(1.d0-hav)**2
!      cxi2=-0.25*(1.d0-eps)**2/(dsqrt(1.d0-eps**2)+1.d0-eps**2)
!      cxi4=chrtav**2/(1.d0-hav)-1.d0
!
      return
      end



      subroutine TorGa_getftrap

c      use precision_mod
c      use real_mod
      implicit none

      include 'globcd.h'
      include 'globcd2.h'
      REAL*8 hbar,c2,c4
      common /cmbftrap/hbar,c2,c4
!---------------------------------------------------------------------------
      REAL*8 tolft
      parameter (tolft=1.d-6)
      REAL*8 ap,ss
      integer jromb
      external TorGa_qftint
!
      hbar=hav
      c2=cxi2
      c4=cxi4
!
      ap=0.5d0*dsqrt(1.d0-hav)*hsqav/hav**2*     
     &    (2.d0+hav)-hsqav/hav**2+1.d0
!
      call TorGa_mqromb1(TorGa_qftint,0.d0,1.d0,ss,tolft,jromb)
      ft=ap-0.75d0*dsqrt(1.d0-hav)*hsqav*ss
      fc=1.d0-ft
!
      return
      end


      subroutine TorGa_getlims(etmax,epst1,epst2)

c      use precision_mod
c      use real_mod
      implicit none

      REAL*8 etmax,epst1,epst2
!  real gammin,gammax,etmax,epst1,epst2
!  real etcutoff
!  real tau,yval,enzsq,hloc
      include 'globcd.h'
!------------------------------------------------------------------------
      REAL*8 gam1,gam2,gam3,gam4,ettmp
      REAL*8 xisqc
!
!     find the intersections of resonance curve with u_perp=0
      call TorGa_gamsrc(gam1,gam2,1.d0)
!
      if (gam2 .le. 1.d0)then     ! no resonance set etmax negative
         gammin=1.d0
         gammax=1.d0
         etmax=-1.d0
         epst1=-1.d0
         epst2=-1.d0
         return
      endif
!
!     in the case of gam2 > 1.0
      if (gam1 .le. 1.d0)then
         gammin=gam2
         etmax=etcutoff
         gammax=etmax*tau+gammin
      else
         gammin=gam1
         ettmp=(gam2-gam1)/tau
         etmax=DMIN1(etcutoff,ettmp)
         gammax=etmax*tau+gammin
      endif
!
!     find the intersections of resonance curve with passing- 
!     trapped separatrix
      xisqc=1.d0-hloc
      call TorGa_gamsrc(gam3,gam4,xisqc)
      epst1=(gam3-gammin)/tau
      epst2=(gam4-gammin)/tau
      return
!
      end


      subroutine TorGa_getfcap(u,fcap,fprime)

c      use precision_mod
c      use real_mod
      implicit none

      REAL*8 u,tol
      REAL*8 fcap,fprime
       include 'globcd.h'
       include 'globcd2.h'
!----------------------------------------------------------------------------
      REAL*8 uv,rhocap,gamma,ss
      integer jromb
      common /cmbqfc/uv,rhocap,gamma
      external TorGa_fcapint
!
      if (u .le. 0.d0)then
         fcap=0.d0
         fprime=0.d0
         return
      endif
!
      if (igv .eq. 0)then
         fcap=1.d0/(fc*(zrat+4.d0))*u**4
         fprime=1.d0/(fc*(zrat+4.d0))*(4.d0*dsqrt(1.d0+u**2)*u**2)
      else
         uv=u
         rhocap=zrat
         gamma=dsqrt(1.d0+u*u)
         tol=tolval
         call TorGa_mqromb1(TorGa_fcapint,0.d0,1.d0,ss,tol,jromb)
!%PR         write (6,'(a,e12.4,i5)')'getfcap:  /u/jromb = ',u,jromb
         fcap=u**4/fc*ss
         fprime=(u/gamma)**2/fc-rhocap*fcap/u**2
      endif
!
      return
      end


      subroutine TorGa_gethcap(z,hcap,hprime)

c      use precision_mod
c      use real_mod
      implicit none

      REAL*8 z,hcap,hprime
      include 'globcd.h'
!----------------------------------------------------------------------------
      REAL*8 s,ss,ap
      integer jromb
      external TorGa_hcapint
!
      s=z*(1.d0-hav)/(1.d0-z*hav)
      hprime=-0.5d0/dsqrt(1.d0-z*hav)/dsqrt(1.d0+s**2* 
     &     (cxi2*(1.d0-s**2)+cxi4*s**2) )
      if (z .ge. 1.0d0)then
         hcap=0.d0
         return
      endif
!
      ap=dsqrt(1.d0-hav)/hav*(1.d0/dsqrt(1.d0-(1.d0-s)*hav)-1.d0)
      call TorGa_mqromb1(TorGa_hcapint,s,1.d0,ss,tolval,jromb)
!%PR      write (6,'(a,e12.4,i5)')'gethcap:   /z/jromb = ',z,jromb
      hcap=0.5d0*dsqrt(1.d0-hav)*ss+ap
!
      return
      end

      SUBROUTINE TorGa_mqromb1(func,a,b,ss,eps,jt)
c      use precision_mod
c      use real_mod
      implicit none

      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,ss, eps
      integer jt
      EXTERNAL func
      PARAMETER (JMAX=200, JMAXP=JMAX+1, K=5, KM=K-1)
!CU    USES TorGa_polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      if (a .eq. b)then
         ss=0.d0
         return
      endif
      h(1)=1.d0
      do 11 j=1,JMAX
         jt=j
        call TorGa_trapzd(func,a,b,s(j),j)
!        write (6,*)j
        if (j.ge.K) then
          call TorGa_polint(h(j-KM),s(j-KM),K,0.0d0,ss,dss)
          if (dabs(dss).le.eps*dabs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d0*h(j)
11    continue
      write (6,'(a)') 'WARNING:   too many steps in mqromb'
      END

      SUBROUTINE TorGa_trapzd(func,a,b,s,n)
c      use precision_mod
c      use real_mod
      implicit none

      INTEGER n
      REAL*8 a,b,s
      REAL*8 func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      END
!------------------------------------------------------------------------
      SUBROUTINE TorGa_polint(xa,ya,n,x,y,dy)
c      use precision_mod
c      use real_mod
      implicit none

      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0)then 
             write (6,'(a)') 'WARNING:  failure in TorGa_polint'
             return
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

      real*8 function TorGa_qftint(s)

c      use precision_mod
c      use real_mod
      implicit none

      REAL*8 s
      REAL*8 hbar,c2,c4
      common /cmbftrap/hbar,c2,c4     
!----------------------------------------------------------------------
      REAL*8 dd
      dd=c2*(1.d0-s**2)+c4*s**2
      TorGa_qftint=-dd*s**3/dsqrt(1.d0-hbar+s*hbar)**5/  
     &     (dsqrt(1.d0+dd*s**2)*(1.d0+dsqrt(1.d0+dd*s**2)))
!
      return
      end

      subroutine TorGa_gamsrc(gam1,gam2,xisq)

c      use precision_mod
c      use real_mod
      implicit none

      REAL*8 gam1,gam2
      REAL*8 xisq
      include 'globcd.h'
!     return the gamma values on the resonace curve for given yval,
!     enzsq, and xisq
!     yval=nharm*(omega_c)/omega
!     enzsq=(n_parallel)^2
!     xisq=(u||/u)^2
!-----------------------------------------------------------------------
      REAL*8 rutsq,dd,rut,gtmp
!
      rutsq=xisq*enzsq*(yval**2-1.d0+enzsq*xisq)
!
!     if no physical roots, returns gam1=gam2=-1.d10
      if (rutsq .lt. 0.d0)then  ! no physcal roots
         gam1=-1.d10
         gam2=-1.d10

         return
      endif
!
      dd=1.d0/(1.d0-xisq*enzsq)
      rut=dsqrt(rutsq)
      gam1=dd*(yval-rut)
      gam2=dd*(yval+rut)
!     arrange roots in increasing order
      if (gam2 .lt. gam1)then
         gtmp=gam2
         gam2=gam1
         gam1=gtmp
      endif
!
      return
      end




      real*8 function TorGa_fcapint(y)

c      use precision_mod
c      use real_mod
      implicit none 

      REAL*8 y
      REAL*8 uv,rhocap,gamma
      common /cmbqfc/uv,rhocap,gamma
!-------------------------------------------------------------------------
      REAL*8 gamy
!
!     gamma=dsqrt(1.d0+uv**2)
      gamy=dsqrt(1.d0+(uv*y)**2)
      TorGa_fcapint=y**(rhocap+3.d0)*((1.d0+gamma)/
     &              (1.d0+gamy))**rhocap/gamy**3
!
      return
      end

      real*8 function TorGa_hcapint(s)
 
c      use precision_mod
c      use real_mod
      implicit none

      REAL*8 s
      include 'globcd.h'
!     return the integrand for evaluating hcap
!-------------------------------------------------------------------------
      REAL*8 qq
      qq=s**2*(cxi2*(1.d0-s**2)+cxi4*s**2)
      TorGa_hcapint=-qq/dsqrt(1.d0-hav+s*hav)**3  
     &     /(dsqrt(1.d0+qq)*(1.d0+dsqrt(1.d0+qq)) )
!
      return
      end



! N Bertelli 14 July 2014:following 3 lines commented
!#######################################################################
!!!      INCLUDE 'const_and_precisions.f90'
!!!      INCLUDE 'config_ext.f90'
!!!      INCLUDE 'green_func_ext.f90'
!#######################################################################

      SUBROUTINE green_func_emp(igv,Te,Zeff,fc,u,gam, K,dKdu) 
!=======================================================================
! Author:  N.B.Marushchenko
! March 2010: prepared for TORBEAM
!
! Returns the 'velocity part' of the Green's function
! Actually, is only the shell for recalling the original subroutines.
!=======================================================================
! INPUTS:
!  igv      - switcher for the models (emp definitions)
!  adj_appr - switcher for the models (nm definitions)
!  Te       - temperature [keV]
!  Zeff     - effective charge
!  fc       - fraction of circulating particles
!  u        - p/sqrt(2mT)
!  gam      - relativistic factor
!
! OUTPUTS:
!   K   - Spitzer's function
!  dKdu = dK/du, i.e. its derivative over normalized momentum
!=======================================================================
! USE precision_mod         ! emp
! USE real_mod              ! emp
!---
      USE const_and_precisions  ! NM
      USE green_func_ext        ! NM
!---
      IMPLICIT NONE
!--- remove it later! ---
      INTEGER, PARAMETER :: p_ = 8
!---
      INTEGER,          INTENT(in)  :: igv
      REAL(p_),         INTENT(in)  :: Te,Zeff,fc,u,gam
      REAL(p_),         INTENT(out) :: K,dKdu
!--- internal variables
      REAL(wp_) :: SS1,ne1,Te1,Zeff1,fc1,u1,q1,gam1
      REAL(wp_) :: K1,dKdu1
      CHARACTER(Len=1) :: adj_appr(6)
      LOGICAL, SAVE :: first =.true.
!=======================================================================
!--- Spitzer function definitions ---
      adj_appr(1) = 'l'         ! collisionless limit
! adj_appr(1) = 'c'         ! collisional (classical) limit, w/o trap. part.
      adj_appr(2) = 'm'         ! momentum conservation
! adj_appr(2) = 'h'         ! high-speed limit
!--- Green's function (pitch-angle part) ---
      adj_appr(3) = 'l'         ! DO NOT CHANGE!
!--- Spitzer function approach: relativistic or non-relat. ---
      adj_appr(4) = 'r'         ! relativistic formulation
! adj_appr(4) = 'n'         ! non-relativistic formulation
!--- Spitzer function (method) ---
      adj_appr(5) = 'v'         ! DO NOT CHANGE!
!---
      adj_appr(6) = 'i'         ! DO NOT CHANGE!
!=======================================================================
      Te1   = Te
      Zeff1 = Zeff
      fc1   = fc
      q1    = u
      u1    = u/sqrt(2*Te/mc2_)
      gam1  = gam
!---
      IF (first) THEN
       CALL Setup_SpitzFunc(adj_appr)
       first =.false.
      ENDIF
!---
      CALL SpitzFuncCoeff(Te1,Zeff1,fc1)
      CALL GenSpitzFunc  (Te1,Zeff1,fc1,u1,q1,gam1, K1,dKdu1)
! CALL GreenFunction (SS,Te1,Zeff1,b,bav,b2av,ft, &
!                            gam,qq,qpar,lamb, X,dXdw,dXdqpar)
!---
      K   =  K1
      dKdu = dKdu1
!emp: added for consistency with curgap by Lin-Liu (where u_LL=p/mc):
!     getfcap in curgap returns (vth/c)^4*K/fc and 
!     (vth/c)^4*gamma/u_LL*d(K/fc)/du_LL
      K = K/fc*(2*Te/mc2_)**2
      dKdu = dKdu/fc/u1*gam*(2*Te/mc2_)
!=======================================================================
      RETURN
      END SUBROUTINE green_func_emp

!#######################################################################
