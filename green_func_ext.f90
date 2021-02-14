!#######################################################################

 MODULE green_func_ext

!#######################################################################
!
! The module contains few subroutines which are requested to calculate 
! the current drive value by adjoint approach
!
!#######################################################################
 USE const_and_precisions
! USE rayinfo_data
! USE config
! USE kinetic
! USE mymathlib
!------- necessary for gfunc (solver SYNCH of S.Kasilov)
! USE store_gfunc_mod, only : prop
!-------
 USE config_ext  ! for use out of TRAVIS
!---
 IMPLICIT NONE
 CHARACTER(Len=1), PRIVATE :: adj_appr(6)    ! adjoint approach switcher
!-------
 REAL(wp_), PRIVATE :: r2,q2,gp1,Rfactor
!--- for H.M. subroutines (operation with DKES-data) ---
! CHARACTER(Len=lfn_), PRIVATE    :: name_dk
! INTEGER,   PRIVATE, PARAMETER   :: nfc_dk = 30  ! number of rad. points in DKES-data
! INTEGER,   PRIVATE, PARAMETER   :: mnit = 20    ! upper number of iterat.
! INTEGER,   PRIVATE, PARAMETER   :: ndv  = 200   ! velocity array dimens. (< 500)
! INTEGER,   PRIVATE, PARAMETER   :: ndl  = 200   ! magn. moment array dimens. (< 500)
! REAL(dp_), PRIVATE, PARAMETER   :: vnu  = umax_ ! upper velocity value 
! REAL(dp_), PRIVATE, ALLOCATABLE :: avn (:)      ! velocity array
! REAL(dp_), PRIVATE, ALLOCATABLE :: aftp(:)      ! trap. fraction array
! REAL(dp_), PRIVATE, ALLOCATABLE :: agsf(:)      ! gen. Spitz.func. array
! REAL(dp_), PRIVATE, ALLOCATABLE :: almb(:)      ! magn. moment array
!-------
 REAL(wp_), PRIVATE, ALLOCATABLE :: alam(:)      ! magn. moment array
 REAL(wp_), PRIVATE, ALLOCATABLE :: aun (:)      ! velocity array
 REAL(wp_), PRIVATE, ALLOCATABLE :: sft (:)      ! trap. fraction array (double)
 REAL(wp_), PRIVATE, ALLOCATABLE :: sfn (:)      ! Spitzer function (double)
 REAL(wp_), PRIVATE, ALLOCATABLE :: fe2 (:,:)    ! e/e coll.freq-cy multipl. by u^2
!-------
 REAL(wp_), PRIVATE, PARAMETER :: delta = 1e-4   ! border for recalculation
!-------
 REAL(wp_), PRIVATE :: sfd(1:4)
 INTEGER,   PRIVATE, PARAMETER :: nre = 2        ! order of rel. corr.
 REAL(wp_), PRIVATE, PARAMETER :: vp_mee(0:4,0:4,0:2) = &
            RESHAPE((/0.0,      0.0,      0.0,      0.0,      0.0,      &
                      0.0,      0.184875, 0.484304, 1.06069,  2.26175,  &
                      0.0,      0.484304, 1.41421,  3.38514,  7.77817,  &
                      0.0,      1.06069,  3.38514,  8.73232,  21.4005,  &
                      0.0,      2.26175,  7.77817,  21.4005,  55.5079,  &
                      !                                                 &
                      0.0,     -1.33059,-2.57431, -5.07771, -10.3884,   &
                     -0.846284,-1.46337, -1.4941, -0.799288, 2.57505,   &
                     -1.1601,  -1.4941,   2.25114,  14.159,   50.0534,  &
                     -1.69257, -0.799288, 14.159,   61.4168,  204.389,  &
                     -2.61022,  2.57505,  50.0534,  204.389,  683.756,  &
                      !                                                 &
                      0.0,      2.62498,  0.985392,-5.57449, -27.683,   &
                      0.0,      3.45785,  5.10096,  9.34463,  22.9831,  &
                     -0.652555, 5.10096,  20.5135,  75.8022,  268.944,  &
                     -2.11571,  9.34463,  75.8022,  330.42,   1248.69,  &
                     -5.38358,  22.9831,  268.944,  1248.69,  4876.48/),&
                    (/5,5,3/))
 REAL(wp_), PRIVATE, PARAMETER :: vp_mei(0:4,0:4,0:2) = &
            RESHAPE((/0.0,     0.886227, 1.0,      1.32934,  2.0,       &
                      0.886227,1.0,      1.32934,  2.0,      3.32335,   &
                      1.0,     1.32934,  2.0,      3.32335,  6.0,       &
                      1.32934, 2.0,      3.32335,  6.0,      11.6317,   &
                      2.0,     3.32335,  6.0,      11.6317,  24.0,      &
                      !                                                 &
                      0.0,     0.332335, 1.0,      2.49251,  6.0,       &
                      1.66168, 1.0,      2.49251,  6.0,      14.5397,   &
                      3.0,     2.49251,  6.0,      14.5397,  36.0,      &
                      5.81586, 6.0,      14.5397,  36.0,     91.5999,   &
                      12.0,    14.5397,  36.0,     91.5999,  240.0,     &
                      !                                                 &
                      0.0,    -0.103855, 0.0,      1.09047,  6.0,       &
                      0.726983,0.0,      1.09047,  6.0,      24.5357,   &
                      3.0,     1.09047,  6.0,      24.5357,  90.0,      &
                      9.81427, 6.0,      24.5357,  90.0,     314.875,   &
                      30.0,    24.5357,  90.0,     314.875,  1080.0 /), &
                    (/5,5,3/))
 REAL(wp_), PRIVATE, PARAMETER :: vp_oee(0:4,0:4,0:2) = &
            RESHAPE((/0.0,     0.56419,  0.707107, 1.0073,   1.59099,   &
                      0.56419, 0.707107, 1.0073,   1.59099,  2.73981,   &
                      0.707107,1.0073,   1.59099,  2.73981,  5.08233,   &
                      1.0073,  1.59099,  2.73981,  5.08233,  10.0627,   &
                      1.59099, 2.73981,  5.08233,  10.0627,  21.1138,   &
                      !                                                 &
                      0.0,     1.16832,  1.90035,  3.5758,   7.41357,   &
                      2.17562, 1.90035,  3.5758,   7.41357,  16.4891,   &
                      3.49134, 3.5758,   7.41357,  16.4891,  38.7611,   &
                      6.31562, 7.41357,  16.4891,  38.7611,  95.4472,   &
                      12.4959, 16.4891,  38.7611,  95.4472,  244.803,   &
                      !                                                 &
                      0.0,     2.65931,  4.64177,  9.6032,   22.6941,   &
                      4.8652,  4.64177,  9.6032,   22.6941,  59.1437,   &
                      9.51418, 9.6032,   22.6941,  59.1437,  165.282,   &
                      21.061,  22.6941,  59.1437,  165.282,  485.785,   &
                      50.8982, 59.1437,  165.282,  485.785,  1483.22/), &
                    (/5,5,3/))
 REAL(wp_), PRIVATE, PARAMETER :: vp_g(0:4,0:2) = &
            RESHAPE((/1.32934, 2.0,      3.32335,  6.0,      11.6317,   &
                      2.49251, 0.0,      2.90793,  12.0,     39.2571,   &
                      1.09047, 6.0,      11.45,    30.0,     98.9606/), &
                    (/5,3/))
!--- for SYNCH solver (fully relativistic Green's function) ---
! CHARACTER(Len=lfn_), PRIVATE :: name_gfunc
!########################################################################

 CONTAINS

!#######################################################################
 
! SUBROUTINE Setup_SpitzFunc(Scen)
 SUBROUTINE Setup_SpitzFunc(adj_appr1)
!=======================================================================
! USE plasma_profiles, ONLY : ne_prof,Te_prof,Zeff_prof 
 IMPLICIT NONE
! TYPE(Scenario), INTENT(in) :: Scen
 CHARACTER(Len=1), INTENT(in) :: adj_appr1(6)
! INTEGER   :: ierr,i,n,ibuff
! REAL(wp_) :: rmie,SS,dumm
! LOGICAL   :: renew
!=======================================================================
 adj_appr(:) = adj_appr1(:)
 RETURN
!=======================================================================
! adj_appr(:) = Scen%adj_appr(:)
!=======================================================================
! IF (adj_appr(5) == 's') adj_appr(4) = 'n'  !<- temporal definition
!=======================================================================
! IF (adj_appr(3) == 'K') THEN
! ENDIF
!=======================================================================
! IF (adj_appr(1) /= 'g' .and. adj_appr(5) /= 's') RETURN
!========================================================================
! loading the DKES database
!=======================================================================
! name_dk = Scen%name_dk
!---
! IF (adj_appr(1) == 'g') CALL tc_dkes_l(nfc_dk,name_dk)
!=======================================================================
! arrays of the velocity and the "effective trapped particle" fraction
! (subroutines of H. Maassberg)
!=======================================================================
! IF (ALLOCATED(avn)) DEALLOCATE(avn,aftp,aun,sft)
! IF (.NOT.ALLOCATED(avn)) THEN
! ENDIF
!!---
! CALL gsppr_v(ndv,vnu,avn)
! aun(:) = avn(:)
!!---
! aftp(:) = 0;  sft(:) = 0
!=======================================================================
! arrays needed for the generalized Spitzer function: 
! direct solution of the kinetic equation (subroutines of H. Maassberg)
!=======================================================================
! IF (adj_appr(5) == 's') THEN
!   IF (adj_appr(4) == 'n') THEN
!   ENDIF
! ENDIF
!=======================================================================
! arrays needed for the model with barely-trapped part. included
! (subroutines of H. Maassberg)
!=======================================================================
! IF (adj_appr(3) == 'b') THEN
! ENDIF
!=======================================================================
 RETURN
 END SUBROUTINE Setup_SpitzFunc

!########################################################################
!########################################################################
!########################################################################

! SUBROUTINE GreenFunction(RayInfo,gam,qq,qpar,lamb, X,dXdw,dXdqpar) 
 SUBROUTINE GreenFunction(SS,Te,Zeff,b,bav,b2av,ft, &
                          gam,qq,qpar,lamb, X,dXdw,dXdqpar) 
!=======================================================================
! Returns X, dX/dw and dX/dqpar, i.e. derivatives of the Green function
! (X=G/Fm is the solution of the adjoint kinetic equation for G),
!
! X =-upar/|upar| * H(p) * K(u),
! 
! where fc = 1-ft  - fraction of circulating particles,
!       H(p)       - angle-part of the Green-function,
!       p = upar/u - pitch of the momentum 
!       K(u)       - Spitzer function.
!
! Refs.: M. Rome' et al., Plasma Phys. Contr. Fus. 40 (1998) 511;
!        Y.R. Lin-Liu et al., Phys. Plasmas 10 (2003) 4064.
!
! Inputs:
! RayInfo     - structure with the local plasma and config. parameters
!---
! SS          - flux-surface label
! Te          - temperarure, keV
! Zeff        - Zeff
! b           = B/Bmax (Bmax on the given surface)
! bav         = <b>
! b2av        = <b^2>
! ft          - trapped particles fraction
!---
! qq          = (p/mc)**2
! qpar        - normalized parallel momentum, with q=p/mc
! lamb        = (1-(qpar/q)**2)/b, i.e. normalized magnetic moment
!
! Outputs:
!  X          - Green-function
! dXdw        = dX/dwr, i.e. derivative over the norm. energy
! dXdqpar     = dX/dqpar
!=======================================================================
! July 2009: revised version, which is in agreement with both the
!            SYNCH solver (S. Kasilov) and CQL3D (R. Harvey & R. Prater).
!=======================================================================
 IMPLICIT NONE
! TYPE(RayInfolocal),INTENT(in)  :: RayInfo
 REAL(wp_), INTENT(in)  :: SS,Te,Zeff,b,bav,b2av,ft
 REAL(wp_), INTENT(in)  :: gam,qq,qpar,lamb
 REAL(wp_), INTENT(out) :: X,dXdw,dXdqpar
 REAL(wp_) :: sgn,cs,q,u,lamb1,uu,uper2,upar,mu
 REAL(wp_) :: bet,modB,maxB,minB,b1
 REAL(wp_) :: fc,fteff,fceff,dftdu,rc,drc
 REAL(wp_) :: K,dKdu,dKdw
 REAL(wp_) :: H,dHdlm,H0,dH0dlm,dH0du,H1,dH1dlm,dH1du,dHdu,dHdw
 REAL(wp_) :: uper,dXduper,dXdupar,const
 INTEGER   :: ii,flag
!=======================================================================
  X      = 0
 dXdw    = 0
 dXdqpar = 0
!---
 IF (qq < comp_eps) RETURN
!=======================================================================
! fully relativistic model of S. Kasilov (Phys. Plasmas 3, 1996, p.4115)
!=======================================================================
! IF (adj_appr(3) == 'K') THEN
!   SS   = RayInfo%SS
!   b    = RayInfo%b
!   bav  = RayInfo%bav
!   b2av = RayInfo%b2av
!   bet  = RayInfo%bet
!   upar = qpar/bet
!   uper = sqrt(qq-qpar**2)/bet
!!print*,'NEO2:'
!!   CALL gfunc_NEO2(uper,-upar, X,dXduper,dXdupar)
!!print*,'Synch:'
!!   prop=.true.  ! in this case to recalculate
!!   IF (lamb >= 1) RETURN
!   CALL gfunc_synch(name_gfunc,SS,b,uper,-upar, X,dXduper,dXdupar)
!!---
!   const   = b2av/bav
!    X      = const*X; 
!   dXduper = const*dXduper; 
!   dXdupar =-const*dXdupar
!!---
!   dXdw    = gam /(uper+comp_eps) * dXduper/bet**2
!   dXdqpar =-upar/(uper+comp_eps) * dXduper/bet + dXdupar/bet
!   RETURN
! ENDIF
!=======================================================================
! contribution from the barely trapped electrons (not ready)
!=======================================================================
! IF (lamb >= 1) THEN
!   SELECT CASE(adj_appr(1))
!     CASE('l')
!!      ...
!       RETURN 
!     CASE('g')
!       IF (adj_appr(3) /= 'b') RETURN
!!      ...
!       RETURN 
!     CASE default
!       PRINT*,'GreenFunction: something is wrong, lambda > 1'
!       RETURN
!   END SELECT
! ENDIF
!=======================================================================
! SS   = RayInfo%SS     ! normalized magnetic flux as the surface label
! Te   = RayInfo%Te     ! Te in keV
! Zeff = RayInfo%Zeff   ! >= 1
! mu   = RayInfo%mu     ! mc**2 / Te
! bet  = RayInfo%bet    ! vte / c
! b    = RayInfo%b      ! B / Bmax
! b2av = RayInfo%b2av   ! <b**2>
! ft   = RayInfo%ft     ! trapped particle fraction
 mu   = mc2_/Te
 bet  = sqrt(2/mu)     !
 fc   = 1-ft
!---
 q  = sqrt(qq)  ! p/mc
 u  = q/bet     ! p/pte
 cs = qpar/q    ! pitch
!---
 lamb1 = lamb
 b1 = b
 IF (adj_appr(1) == 'c') b1 = 1
!=======================================================================
! Cohen's model (subroutines "F_Cohen" and "H_Cohen" got from D.Farina)
!=======================================================================
! IF (adj_appr(2) == 'h' .and. adj_appr(6) == 'c') THEN
! ENDIF
!=======================================================================
! pitch-angle part
!=======================================================================
 dHdu = 0
 SELECT CASE(adj_appr(1))
   CASE('c')                  !---------- classical limit -------------!
      H    = abs(cs)                                                   !
     dHdlm =-0.5/abs(cs)                                               !
     fc    = 1                                                         !
   CASE('l')                  !--------- collisionless limit ----------!
     IF (lamb <= 1) THEN                                               !
        H    = b2av/(2*fc)*Conf_p_avrg_integr(SS,lamb1)                !
       dHdlm =-b2av/(2*fc)/Conf_p_avrg       (SS,lamb1)                !
     ELSE                                                              !
        H    = 0                                                       !
       dHdlm = 0                                                       !
     ENDIF                                                             !
!   CASE('g')                  !--------- generalized approach ---------!
!     IF (adj_appr(3) == 'b') THEN !---- with barely trapped part. -----!
!     ENDIF      !------------------------------------------------------!
   CASE default
     PRINT*,'GenSpitzFunc: WARNING: Spitzer function is not defined.'
     RETURN
 END SELECT
 dHdw = gam/(q*bet)*dHdu
!=======================================================================
! velocity part (Spitzer function)
!=======================================================================
 CALL GenSpitzFunc(Te,Zeff,fc,u,q,gam, K,dKdu)
 dKdw = gam/(q*bet) * dKdu
!=======================================================================
! final definitions, (w,lamb) -> (w,qpar)
!=======================================================================
 sgn =-sign(unit,qpar)
!---
  X      = sgn *  H *  K
 dXdw    = sgn * (H * dKdw + 2*gam/b1*(cs/q)**2*dHdlm*K + dHdw*K)
 dXdqpar = sgn * dHdlm*(-2*cs/(b1*q))*K
!=======================================================================
 RETURN
 END SUBROUTINE GreenFunction

!#######################################################################
!#######################################################################
!#######################################################################
 
 SUBROUTINE GenSpitzFunc(Te,Zeff,fc,u,q,gam, K,dKdu)

!=======================================================================
! Author:  N.B.Marushchenko
! June 2005: as start point the subroutine of Ugo Gasparino (198?) 
!            SpitzFunc() is taken and modified.
!            1. adapted to the Fortran-95
!            2. derivative of Spitzer function is added
!            3. separation for 2 brunches is done: 
!               1st is referenced as 'with conservation of the moment',
!               2nd - as 'high speed limit'. 
!               The last one is taken from the Lin-Liu formulation 
!               (Phys.Plasmas 10 (2003) 4064) with K = F*fc.
!               The asymptotical high speed limit (Taguchi-Fisch model)
!               is also included as the reference case.
! Feb. 2008: non-relativ. version is replaced by the weakly relativistic;
!            the method is the the same, but the trial-function is 
!            based on the relativistic formulation.
!            The relativistic corrections for the collisional operator 
!            up to the second order, i.e. (1/mu)**2, are applied.
! Sep. 2008: generalized Spitzer function for arbitrary collisionality
!            is implemented. The model is based on the concept of 
!            the "effective trapped particles fraction". 
!            The different.-integral kinetic equation for the generalized
!            Spitzer function is produced with help of subroutines 
!            ArbColl_TrappFract_Array and ArbColl_SpitzFunc_Array,
!            where the subroutines of H. Maassberg are called).
!========================================================================
! Spitzer function with & w/o trapped particle effects is given by:
!
!    K(x) = x/gamma*(d1*x+d2*x^2+d4*x^3+d4*x^4),
! 
! where x = v/v_th and gamma=1 for non-relativistic version (Ugo),
! or    x = p/p_th for relativistic version (N.M., February 2008).
! Note, that somewhere the function F(x) instead of K(x) is applied,
!
!    F(x) = K(x)/fc.
!
! Numerical inversion of the 5x5 symmetric matrix obtained from the
! generalized Spitzer problem (see paper of Taguchi for the equation
! and paper of Hirshman for the variational approach bringing to the
! matrix to be inverted).
! 
! The numerical method used is an improved elimination scheme
! (Banachiewiczs-Cholesky-Crout method).
! This method is particularly simple for symmetric matrix.
! As a reference see "Mathematical Handbook" by Korn & Korn, p.635-636.
! 
! Refs.:  1. M. Rome' et al., Plasma Phys. Contr. Fus. 40 (1998) 511;
!         2. S.P. Hirshman, Phys. Fluids 23 (1980) 1238
!========================================================================
! INPUTS:
!  SS   - flux-surface label (normalized flux)
!  u    - p/sqrt(2mT)
!  q    - p/mc;
!  gam  - relativistic factor;
!  mu   - mc2/Te
!  Zeff - effective charge;
!  fc   - fraction of circulating particles.
!
! OUTPUTS:
!   K   - Spitzer's function
!  dKdu = dK/du, i.e. its derivative over normalized momentum
!=======================================================================
 IMPLICIT NONE
 REAL(wp_), INTENT(in)  :: Te,Zeff,fc,u,q,gam
 REAL(wp_), INTENT(out) :: K,dKdu
 REAL(wp_) :: mu,gam1,gam2,gam3,w,dwdu
!=======================================================================
  K   = 0
 dKdu = 0
 IF (u < comp_eps) RETURN
!---
 mu = mc2_/max(Te,1d-3)
 SELECT CASE(adj_appr(2))
   CASE('m')  !--------------- momentum conservation ------------------!
     IF (adj_appr(5) == 'v') THEN !---- method: variat. principle -----!
       gam1 = gam                                                      !
       IF (adj_appr(4) == 'n') gam1 = 1                                !
       gam2 = gam1*gam1                                                !
       gam3 = gam1*gam2                                                !
        K   = u/gam1*u*(sfd(1)+u*(sfd(2)+u*(sfd(3)+u*sfd(4))))         !
       dKdu = u/gam3*  (sfd(1)*(1+  gam2)+u*(sfd(2)*(1+2*gam2)+ &      !
                     u*(sfd(3)*(1+3*gam2)+u* sfd(4)*(1+4*gam2))))      !
!     ELSEIF (adj_appr(5) == 's') THEN !-- method: solv. of kin. eqn. --!
     ENDIF                        !---------- end of method -----------!
     !--------------------- end momentum conservation -----------------!
   CASE('h')  !---------------- high-speed-limit ----------------------!
     IF (adj_appr(4) == 'n') THEN     !- non-relativ. asymptotic form -!
       K   =   u**4 *fc/(Zeff+1+4*fc) !- (Taguchi-Fisch model)        -!
      dKdu = 4*u**3 *fc/(Zeff+1+4*fc)                                  !
     ELSEIF (adj_appr(4) == 'r') THEN !- relativistic, Lin-Liu form.  -!
       CALL SpitzFunc_HighSpeedLimit(Zeff,fc,u,q,gam, K,dKdu)          !
     ENDIF                                                             !
   CASE default   !----------------------------------------------------!
     PRINT*,'GenSpitzFunc: WARNING: Spitzer function is not defined.'
     RETURN
 END SELECT
!=======================================================================
!=======================================================================
 RETURN
 END SUBROUTINE GenSpitzFunc

!#######################################################################
!#######################################################################
!#######################################################################

! SUBROUTINE SpitzFuncCoeff(SS,ne,Te,Zeff,fc)
 SUBROUTINE SpitzFuncCoeff(Te,Zeff,fc)
 
!=======================================================================
! Calculates the matrix coefficients required for the subroutine 
! "GenSpitzFunc", where the Spitzer function is defined through the
! variational principle.
! 
! Weakly relativistic (upgraded) version (10.09.2008). 
! Apart of the non-relativistic matrix coefficients, taken from the
! old subroutine of Ugo Gasparino, the relativistic correction written
! as series in 1/mu^n (mu=mc2/T) powers is added. Two orders are taken 
! into account, i.e. n=0,1,2.
!
! In this version, the coefficients "oee", i.e. Omega_ij, are formulated
! for arbitrary collisionality. 
!
! INPUT VARIABLES:
! SS   - flux-surface label
! ne   - density, 1/m^3
! Te   - temperature, keV
! Zeff - effective charge
! fc   - fraction of circulating particles
!
! OUTPUT VARIABLES (defined as a global ones): 
! sfd(1),...,sfd(4) - coefficients of the polynomial expansion of the 
!                 "Spitzer"-function (the same as in the Hirshman paper)
!=======================================================================
 IMPLICIT NONE
! REAL(wp_), INTENT(in) :: SS,ne
 REAL(wp_), INTENT(in) :: Te,Zeff,fc
 INTEGER   :: n,i,j,ij,jj
 REAL(wp_) :: rtc,rtc1,mu,y,tn(1:nre)
 REAL(wp_) :: me(0:4,0:4),mi(0:4,0:4),m(0:4,0:4),g(0:4)
 REAL(wp_) :: om(0:4,0:4)
 REAL(wp_) :: gam11,gam21,gam31,gam41,gam01, &
                    gam22,gam32,gam42,gam02, &
                          gam33,gam43,gam03, &
                                gam44,gam04,gam00
 REAL(wp_) :: alp12,alp13,alp14,alp10, &
                    alp23,alp24,alp20, &
                          alp34,alp30,alp40
 REAL(wp_) :: bet0,bet1,bet2,bet3,bet4,d0
 LOGICAL   :: renew,rel,newTe,newne,newZ,newfc
 REAL(wp_), SAVE :: sfdx(1:4) = 0
 REAL(wp_), SAVE :: ne_old =-1, Te_old =-1, Zeff_old =-1, fc_old =-1
!=======================================================================
 rel   = Te > 1
! newne = abs(ne  -ne_old  ) > delta*ne
 newTe = abs(Te  -Te_old  ) > delta*Te
 newZ  = abs(Zeff-Zeff_old) > delta*Zeff
 newfc = abs(fc  -fc_old  ) > delta*fc
! SELECT CASE(adj_appr(1))
!   CASE ('l','c')
     renew = (newTe .and. rel) .OR. newZ .OR. newfc
!   CASE ('g')
!     renew = newTe .OR. newne .OR. newZ .OR. newfc
! END SELECT
!---
 IF (.not.renew) THEN
   sfd(:) = sfdx(:)
   RETURN
 ENDIF
!=======================================================================
 tn(:) = 0
 IF (adj_appr(4) == 'r') THEN
   IF (nre > 0) THEN
     mu = mc2_/max(Te,1d-3)
     tn(1) = 1/mu
     DO n=2,min(2,nre)
       tn(n) = tn(n-1)/mu
     ENDDO
   ENDIF
 ENDIF
!---
 SELECT CASE(adj_appr(1))
   CASE ('l','c')     !---- both classical & collisionless limits  ----!
     rtc = (1-fc)/fc;   rtc1 = rtc+1                                   !
     !---                                                              !
     DO i=0,4                                                          !
       g(i) = vp_g(i,0)                                                !
       DO n=1,min(2,nre)                                               !
         g(i) = g(i) + tn(n)*vp_g(i,n)                                 !
       ENDDO                                                           !
       !---                                                            !
       DO j=0,4                                                        !
         IF (i == 0 .or. j == 0 .or. j >= i) THEN                      !
             y =      vp_mee(i,j,0) + rtc *vp_oee(i,j,0) + &           !
                                 Zeff*rtc1*vp_mei(i,j,0)               !
           DO n=1,min(2,nre)                                           !
             y = y + (vp_mee(i,j,n) + rtc *vp_oee(i,j,n) + &           !
                                 Zeff*rtc1*vp_mei(i,j,n))*tn(n)        !
           ENDDO                                                       !
           m(i,j) = y                                                  !
         ENDIF                                                         !
       ENDDO                                                           !
     ENDDO                                                             !
     DO i=2,4                                                          !
       DO j=1,i-1                                                      !
         m(i,j) = m(j,i)                                               !
       ENDDO                                                           !
     ENDDO                                                             !
     m(0,0) = 0                                                        !
   CASE ('g')         !----------- arbitrary collisionality -----------!
   CASE default       !------------------------------------------------!
     PRINT*,'Green_Func: WARNING: Adjoint approach is not defined.'
     RETURN
 END SELECT
!=======================================================================
 gam11 = m(1,1)
 gam21 = m(2,1)
 gam31 = m(3,1)
 gam41 = m(4,1)
 gam01 = m(0,1)
!
 alp12 = m(1,2)/m(1,1)
 alp13 = m(1,3)/m(1,1)
 alp14 = m(1,4)/m(1,1)
 alp10 = m(1,0)/m(1,1)
! 
 gam22 = m(2,2)-gam21*alp12
 gam32 = m(3,2)-gam31*alp12
 gam42 = m(4,2)-gam41*alp12
 gam02 = m(0,2)-gam01*alp12
!
 alp23 = gam32/gam22
 alp24 = gam42/gam22
 alp20 = gam02/gam22
!
 gam33 = m(3,3)-gam31*alp13-gam32*alp23
 gam43 = m(4,3)-gam41*alp13-gam42*alp23
 gam03 = m(0,3)-gam01*alp13-gam02*alp23
!
 alp34 = gam43/gam33
 alp30 = gam03/gam33
!
 gam44 = m(4,4)-gam41*alp14-gam42*alp24-gam43*alp34
 gam04 = m(0,4)-gam01*alp14-gam02*alp24-gam03*alp34
!
 alp40 = gam04/gam44
!
 gam00 = m(0,0)-gam01*alp10-gam02*alp20-gam03*alp30-gam04*alp40
!
 bet1 =  g(1)/m(1,1)
 bet2 = (g(2)-gam21*bet1)/gam22
 bet3 = (g(3)-gam31*bet1-gam32*bet2)/gam33
 bet4 = (g(4)-gam41*bet1-gam42*bet2-gam43*bet3)/gam44
 bet0 = (g(0)-gam01*bet1-gam02*bet2-gam03*bet3-gam04*bet4)/gam00
!
 d0     = bet0
 sfd(4) = bet4-alp40*d0
 sfd(3) = bet3-alp30*d0-alp34*sfd(4)
 sfd(2) = bet2-alp20*d0-alp24*sfd(4)-alp23*sfd(3)
 sfd(1) = bet1-alp10*d0-alp14*sfd(4)-alp13*sfd(3)-alp12*sfd(2)
!=======================================================================
 fc_old   = fc
! ne_old   = ne
 Te_old   = Te
 Zeff_old = Zeff
!---
 sfdx(1:4) = sfd(1:4)
!=======================================================================
 RETURN
 END SUBROUTINE SpitzFuncCoeff

!#######################################################################
!#######################################################################
!#######################################################################

 SUBROUTINE SpitzFunc_HighSpeedLimit(Zeff,fc,u,q,gam, K,dKdu)
!=======================================================================
! Calculates the "Spitzer function" in high velocity limit, relativistic 
! formulation: Lin-Liu et al., Phys.Pl. (2003),v10, 4064, Eq.(33).
! 
! Inputs:
! Zeff - effective charge
! fc   - fraction of circulating electrons
! u    - p/(m*vte)
! q    - p/mc
! gam  - relativ. factor
! 
! Outputs:
!  K   - Spitzer function
! dKdu - its derivative
!=======================================================================
 IMPLICIT NONE
 REAL(wp_), INTENT(in)  :: Zeff,fc,u,q,gam
 REAL(wp_), INTENT(out) :: K,dKdu
 INTEGER :: ierr,nfun
 REAL(8) :: gam2,err,flag,Integr
 REAL(8), PARAMETER :: a = 0d0, b = 1d0, rtol = 1d-4, atol = 1d-12
!=======================================================================
 r2 = (1+Zeff)/fc   ! global parameter needed for integrand, HSL_f(t)
!------------------
 IF (u < 1e-2) THEN
    K   =   u**4/(r2+4)
   dKdu = 4*u**3/(r2+4)
   RETURN
 ENDIF
!=======================================================================
 q2  = q*q       ! for the integrand, HSL_f
 gp1 = gam+1     ! ..
!---
 CALL quanc8(HSL_f,zero,unit,atol,rtol,Integr,err,nfun,flag)
!=======================================================================
 gam2 = gam*gam
!---
  K   = u**4 * Integr
 dKdu = (u/gam)**3 * (1-r2*gam2*Integr)
!=======================================================================
 RETURN
 END SUBROUTINE SpitzFunc_HighSpeedLimit

!#######################################################################
!#######################################################################
!#######################################################################

 FUNCTION HSL_f(t) RESULT(f)
!=======================================================================
! Integrand for the high-speed limit approach (Lin-Liu's formulation)
!=======================================================================
 IMPLICIT NONE
 REAL(8), INTENT(in) :: t
 REAL(8) :: f,g
   g = sqrt(1+t*t*q2)
   f = t**(3+r2)/g**3 * (gp1/(g+1))**r2
 END FUNCTION HSL_f

!#######################################################################

 END MODULE green_func_ext

!N Bertelli 14 July 2014: following line commented
!#######################################################################
!!!!! INCLUDE 'quanc8.f90'
!#######################################################################
