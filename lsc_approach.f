
c     Reference: ``Current in wave-driven plasmas,'' by
c     Charles F. F. Karney and Nathaniel J. Fisch,
c     Phys. Fluids {\bf 29} 180-192 (1986) 

c     Review of definitions and conventions of the reference:
c     q         carries sign of electron charge == - e
c     E         is parallel to B, and in the positive direction
c     v_par     is positive in the direction of E and B
c     \Gamma == \ln \Lambda \frac{n q^4}{4\pi \epsilon_0^2 m^2}
c     v_r    == -sign(qE) \sqrt{m\Gamma / \abs{qE} } (a positive number)
c     Dreicer Velocity v_D = - \sqrt{2 + Z} v_r
c     u         v /(\abs{ v_r})
c     u_\par    v_\par / v_rc     u_\perp   v_\perp / \abs{v_r}

c     Following analytical formuar for Ws is for u < 0.5
c     Here:
c     Mu1=u_par/u
c     Mu0=1, Mu2=Mu1**2, Mu4=Mu1**4
c     U0=u, U04=u**4, U06=u**6, U08=u**8,U10=u**10,
c     Z is the effective charge

 
c01  WsloDwn = Mu1 * U04 / (5.+Z)
c02 ^   -    (2. + Z + 3.*Mu2)                                * U06 / 
c03 ^                                     (y 3.*(3.+Z)*(5.+Z) )
c04 ^   + 2.*( (24.+19.*Z+3.*Z**2)*Mu1 + (9.+Z)*Mu3 )         * U08 /
c05 ^                       ( (3.+Z)*(5.+Z)*(7.+3.*Z)*(9.+Z) ) 
c06 ^  -(
c07 ^    (1041.+1864.*Z + 1189.*Z**2 + 316.*Z**3 + 30.*Z**4 ) * Mu0
c08 ^   +( 417.+ 497.*Z +  181.*Z**2 +  21.*Z**3) * 10.       * Mu2
c09 ^   +  (9.+Z)*(13.+3.*Z)                      *  5.       * Mu4
c10 ^   )                                                     * U10 /
c11 ^      ( 5.*(2.+Z)*(3.+Z)*(5.+Z)*(7.+3.*Z)*(9.+Z)*(13.+3.*Z) )
c

c     Keeping the signs straight is a problem.  We say that positive
c     phase velocity tends to drive current in the supportive, or correct
c     direction.  Positive current drive.  Thus, presumably, the Edc
c     is also in the positive direction.  Electrons tend to run away
c     and be accelerated in the opposite direction.  Thus for positive
c     electron velocity which goes with positive phase velocity, the 
c     normalized   u \equiv v/v_r is negative, and the mu for this positive
c     v is negative.
c
c     v  E      u    mu
c     +  +      -     -
c     -  +      +     +
c     +  -      +     +
c     -  -      -     -
c
c     Code is written with the presumtion of E positive and mu minus
c     for positive velocity and mu plus for negative velocity.  If
c     E is negative, the quantitiy muminus becomes +1; muplus becomes -1
c
      REAL*8 FUNCTION RunProb ( uGiven , MuGiven , ZGiven )
c-----calculate runaway probabilty R(u,mu)

c-----from LSC code file fits,F
      IMPLICIT NONE
c-----input
      real*8 
     & uGiven,  ! v /(\abs{ v_r})
     & MuGiven, ! v_par/v
     & ZGiven   ! effective charge
c-----local
      real*8 u,Mu,Z,x

c-----Table 1 from the Karney Fisch article p.190
c     Coefficients for approximation to  runaway probabilty  R(u,mu=1)
c
c     For mu=1 and 1.4 < u < 8    
c     R(U,mu=1)=exp{Sum_i=0,3 {a_i(u-1)**i]/[Sum_i=0,3 {b_i(u-1)**i]}
c     where b_1=1
c     Table 1 gives maximum relative error 1 percent
c
c     For mu=1 and 1 < u < 1.4    
c     R(U,mu=1)=exp{Sum_i=0,3 {a_i(u-1)**i]/[Sum_i=0,3 {b_i(u-1)**i]}
c     where b_1=1
c     Table 1 gives small absolute error but large
c    
c     For u < 1 and all mu R=0    

      REAL*8 
     ^     RpZ01a0,  RpZ01a1,  RpZ01a2,  RpZ01a3,  RpZ01b2,  RpZ01b3 ,
     ^     RpZ02a0,  RpZ02a1,  RpZ02a2,  RpZ02a3,  RpZ02b2,  RpZ02b3 ,
     ^     RpZ05a0,  RpZ05a1,  RpZ05a2,  RpZ05a3,  RpZ05b2,  RpZ05b3 ,
     ^     RpZ10a0,  RpZ10a1,  RpZ10a2,  RpZ10a3,  RpZ10b2,  RpZ10b3
 
      DATA 
     ^     RpZ01a0,  RpZ01a1,  RpZ01a2,  RpZ01a3,  RpZ01b2,  RpZ01b3 /
     ^    -3.68063,  4.23913, -4.55894, -0.39755, -1.22774,  1.41450 /
      DATA
     ^     RpZ02a0,  RpZ02a1,  RpZ02a2,  RpZ02a3,  RpZ02b2,  RpZ02b3 /
     ^    -4.97636,-16.09015,  0.83188,  0.21737,  6.84615, -0.98649 /
      DATA 
     ^     RpZ05a0,  RpZ05a1,  RpZ05a2,  RpZ05a3,  RpZ05b2,  RpZ05b3 /
     ^    -4.27687, -4.33629,  0.30338,  0.05697,  3.21315, -0.47749 /
      DATA
     ^     RpZ10a0,  RpZ10a1,  RpZ10a2,  RpZ10a3,  RpZ10b2,  RpZ10b3 /
     ^    -4.94597, -1.53482,  0.10112,  0.03087,  2.45288, -0.36896 /

      u    = abs(uGiven)
      Mu = + MuGiven
      Z    = ZGiven

c
c     .                                 Electrons of velocity less than 
c     .                                 the runaway velocity simply do not
c     .                                 run away:
      if ( u*u .le. 1. ) then
        RunProb = 0.
        return
      endif

c
c     Backward-running electrons run away easily once u > 1, and it
c     does not depend much on Z.  However, there is no fit given in the 
c     paper, even though a graph is given.  This is MY crude fit to
c     the Fig. 2 of the reference.
      if ( Mu .le. -0.90 ) then
        x = 0.85 * (10./Z)**2 * ( (u*u-1.)/3. )**2
        if (x .gt. 0.85)
     ^                  x = 1.00 - 0.15*dexp(-(x-0.85)**2)
        RunProb = x
        return
      endif

c     No fit given, no graph given for -1 < Mu < 1 , so I set RunProb to 0.0.
      if ( Mu .gt. -0.90 .and. Mu .lt. 0.90 ) then
        RunProb = 0.00
        return
      endif


c     So we are left with u > 1, and Mu = 1. (or at least Mu > 0.90)
c     This uses the fit data from Table I., page 190.
      x = (u - 1.)
      if        ( Z .le. 1.5 ) then
        RunProb = exp (
     ^  ( RpZ01a0 + RpZ01a1*x + RpZ01a2*x**2 + RpZ01a3*x**3 ) / 
     ^  (               1.0*x + RpZ01b2*x**2 + RpZ01b3*x**3 )
     ^                )
        else if ( Z .le. 3.0 ) then
        RunProb = exp (
     ^  ( RpZ02a0 + RpZ02a1*x + RpZ02a2*x**2 + RpZ02a3*x**3 ) / 
     ^  (               1.0*x + RpZ02b2*x**2 + RpZ02b3*x**3 )
     ^                )
        else if ( Z .le. 7.0 ) then
        RunProb = exp (
     ^  ( RpZ05a0 + RpZ05a1*x + RpZ05a2*x**2 + RpZ05a3*x**3 ) / 
     ^  (               1.0*x + RpZ05b2*x**2 + RpZ05b3*x**3 )
     ^                )
        else if ( Z .le. 10. ) then
        RunProb = exp (
     ^  ( RpZ10a0 + RpZ10a1*x + RpZ10a2*x**2 + RpZ10a3*x**3 ) / 
     ^  (               1.0*x + RpZ10b2*x**2 + RpZ10b3*x**3 )
     ^                )
      endif
      return
      END
c
c     -----------------------------------------------------------------|-------
c
      REAL FUNCTION WsloDwn ( uGiven , MuGiven , ZGiven )
c     calulate W_s(u,mu,Z)
c     WsloDwn   is the energy in units of m v_runaway^2 imparted to the
c               electric field by an electron as it slows down.
c-----from LSC code file fits,F
      IMPLICIT NONE
c-----input
      real*8 
     & uGiven,  ! v /(\abs{ v_r})
     & MuGiven, ! v_par/v
     & ZGiven   ! effective charge
c-----local
      real*8 u,Mu,Z,x
     
      REAL*8
     ^      Mu0 , Mu1 , Mu2 , Mu3 , Mu4 , U02, U04, U06, U08, U10

c-----Table 2 from the Karney Fisch article p.190
c
c     Coefficients for apprpximation to W_s(u,mu,z)
c     For mu=1 and 0 < u < 5
c     W_s(U,mu=1)=Sum_i=2,4 {a_i*u**2i]/[Sum_i=0,3{b_i*u**2i]}
c     Maximum relative error = 2 percents

      REAL
     ^     WpZ01a2,  WpZ01a3,  WpZ01a4,  WpZ01b1,  WpZ01b2,  WpZ01b3 ,
     ^     WpZ02a2,  WpZ02a3,  WpZ02a4,  WpZ02b1,  WpZ02b2,  WpZ02b3 ,
     ^     WpZ05a2,  WpZ05a3,  WpZ05a4,  WpZ05b1,  WpZ05b2,  WpZ05b3 ,
     ^     WpZ10a2,  WpZ10a3,  WpZ10a4,  WpZ10b1,  WpZ10b2,  WpZ10b3
      DATA
     ^     WpZ01a2,  WpZ01a3,  WpZ01a4,  WpZ01b1,  WpZ01b2,  WpZ01b3 /
     ^     0.16612, -0.01495,  0.00775,  0.37136,  0.02240,  0.01645 /
      DATA
     ^     WpZ02a2,  WpZ02a3,  WpZ02a4,  WpZ02b1,  WpZ02b2,  WpZ02b3 /
     ^     0.14200, -0.04048,  0.01145,  0.12253,  0.00384,  0.02440 /
      DATA
     ^     WpZ05a2,  WpZ05a3,  WpZ05a4,  WpZ05b1,  WpZ05b2,  WpZ05b3 /
     ^     0.09880, -0.05152,  0.01113, -0.19484,  0.00559,  0.02362 /
      DATA
     ^     WpZ10a2,  WpZ10a3,  WpZ10a4,  WpZ10b1,  WpZ10b2,  WpZ10b3 /
     ^     0.06537, -0.03895,  0.00738, -0.32456,  0.02797,  0.01526 /

c-----table 3 from the Karney Fisch article p.191
c     Coefficients for apprpximation to W_s(u,mu,z)
c     For mu= - 1 and 0 < u < 1
c     W_s(U,mu = -1)=Sum_i=2,5 {a_i*u**2i]}
c     Maximum relative error = 1.5 percents
      REAL
     ^     WmZ01a2,  WmZ01a3,  WmZ01a4,  WmZ01a5 ,
     ^     WmZ02a2,  WmZ02a3,  WmZ02a4,  WmZ02a5 ,
     ^     WmZ05a2,  WmZ05a3,  WmZ05a4,  WmZ05a5 ,
     ^     WmZ10a2,  WmZ10a3,  WmZ10a4,  WmZ10a5
      DATA
     ^     WmZ01a2,  WmZ01a3,  WmZ01a4,  WmZ01a5 /
     ^    -0.16483, -0.13420,  0.15346, -0.24314 /
      DATA
     ^     WmZ02a2,  WmZ02a3,  WmZ02a4,  WmZ02a5 /
     ^    -0.14186, -0.09297,  0.06661, -0.12870 /
      DATA
     ^     WmZ05a2,  WmZ05a3,  WmZ05a4,  WmZ05a5 /
     ^    -0.09975, -0.04781,  0.00606, -0.03545 /
      DATA
     ^     WmZ10a2,  WmZ10a3,  WmZ10a4,  WmZ10a5 /
     ^    -0.06651, -0.02797, -0.00247, -0.00934 /
c
c
      u    = abs(uGiven)
      Mu = + MuGiven
      Z    = ZGiven

c
c
c     If Mu < 0.5, then use the MACSYMA-derived formula of page 191.
      U02 = u*u
      if ( U02 .le. 0.25 ) then

        Mu0 = 1.
        Mu1 = Mu
        Mu2 = Mu*Mu
        Mu3 = Mu*Mu2
        Mu4 = Mu2*Mu2

        U04 = U02*U02
        U06 = U02*U04
        U08 = U04*U04
        U10 = U04*U06
c       --------------------------------------------------------------|
        WsloDwn = Mu1 * U04 / (5.+Z)
     ^   -    (2. + Z + 3.*Mu2)                                * U06 / 
     ^                                     ( 3.*(3.+Z)*(5.+Z) )
     ^   + 2.*( (24.+19.*Z+3.*Z**2)*Mu1 + (9.+Z)*Mu3 )         * U08 /
     ^                       ( (3.+Z)*(5.+Z)*(7.+3.*Z)*(9.+Z) ) 
     ^  -(
     ^    (1041.+1864.*Z + 1189.*Z**2 + 316.*Z**3 + 30.*Z**4 ) * Mu0
     ^   +( 417.+ 497.*Z +  181.*Z**2 +  21.*Z**3) * 10.       * Mu2
     ^   +  (9.+Z)*(13.+3.*Z)                      *  5.       * Mu4
     ^   )                                                     * U10 /
     ^      ( 5.*(2.+Z)*(3.+Z)*(5.+Z)*(7.+3.*Z)*(9.+Z)*(13.+3.*Z) )

c       --------------------------------------------------------------|

        return
      endif

c     So if u > 0.5 but -1 < Mu < 1 , then set WsloDwn to 0.00
c     The paper gives no fits for this region, but Fig. 4 does graph it.
      if ( Mu .gt. -0.90 .and. Mu .lt. 0.90 ) then
        WsloDwn = 0.00
        return
      endif


      x = u*u

c     Mu is +1, so use the Table II. on page 190.
c     This is valid for u < 5.
      if ( Mu .ge. 0.90 .and. u .le. 5.00 ) then
             if ( Z  .le. 1.5  )  then
                WsloDwn =
     ^                 ( WpZ01a2*x**2 + WpZ01a3*x**3 + Wpz01a4*x**4 ) /
     ^            ( 1. + WpZ01b1*x    + WpZ01b2*x**2 + WpZ01b3*x**3 )
        else if ( Z  .le. 3.0  )  then
                WsloDwn =
     ^                 ( WpZ02a2*x**2 + WpZ02a3*x**3 + Wpz02a4*x**4 ) /
     ^            ( 1. + WpZ02b1*x    + WpZ02b2*x**2 + WpZ02b3*x**3 )

        else if ( Z  .le. 7.0  )  then
                WsloDwn =
     ^                 ( WpZ05a2*x**2 + WpZ05a3*x**3 + Wpz05a4*x**4 ) /
     ^            ( 1. + WpZ05b1*x    + WpZ05b2*x**2 + WpZ05b3*x**3 )

        else
                WsloDwn =
     ^                 ( WpZ10a2*x**2 + WpZ10a3*x**3 + Wpz10a4*x**4 ) /
     ^            ( 1. + WpZ10b1*x    + WpZ10b2*x**2 + WpZ10b3*x**3 )
        endif
      return
      endif

c     If Mu = -1 and 0 < u < 1 then use the fit of Table III. on page 191.
      if ( Mu .le.-0.90 .and. u*u .le. 1.00 ) then
             if ( Z  .le. 1.5  )  then
                WsloDwn = WmZ01a2*x**2 
     ^                 + WmZ01a3*x**3 + WmZ01a4*x**4 + WmZ01a5*x**5
        else if ( Z  .le. 3.0  )  then
                WsloDwn = WmZ02a2*x**2 
     ^                 + WmZ02a3*x**3 + WmZ02a4*x**4 + WmZ02a5*x**5
        else if ( Z  .le. 7.0  )  then
                WsloDwn = WmZ05a2*x**2 
     ^                 + WmZ05a3*x**3 + WmZ05a4*x**4 + WmZ05a5*x**5
        else
                WsloDwn = WmZ10a2*x**2 
     ^                 + WmZ10a3*x**3 + WmZ10a4*x**4 + WmZ10a5*x**5
        endif
      return
      endif

c     Unanticipated input parameter, no data; set return to 0.00.
      WsloDwn = 0.00
      return
      END


      SUBROUTINE WsloPrm ( uGiven , MuGiven , ZGiven,  
     ^                     dWsduou, iWhichWay)

c     Provides for a linear interpolation of Zeff in the region 1.0 <Zeff< 10.
c     The interpolation is calculated with nodes at Z = 1, 2, 5, 10 as given
c     in the Karney and Fisch paper.

c     Written by D. Enright, June 1992.

      implicit none
c-----input
      real*8 
     & uGiven,  ! v /(\abs{ v_r})
     & MuGiven, ! v_par/v
     & ZGiven   ! effective charge
c-----output
      INTEGER iWhichWay
      Real*8 dWsduou
      
c-----local
      REAL*8
     ^        Z ,
     ^        dWsduou1, dWsduou2, dWsduou5, dWsduou10
c
      Z    = ZGiven
c
c
      if (Z .le. 1.0) then
         call WsloPrmZ(uGiven, MuGiven, 1.0, dWsduou, iWhichWay)
      else if (Z .le. 2.0) then
         call WsloPrmZ(uGiven, MuGiven, 1.0, dWsduou1, iWhichWay)
         call WsloPrmZ(uGiven, MuGiven, 2.0, dWsduou2, iWhichWay)
         dWsduou = (dWsduou2 - dWsduou1)/(2.0-1.0)*(Z-1.0) + dWsduou1
      else if (Z .le. 5.0) then
         call WsloPrmZ(uGiven, MuGiven, 2.0, dWsduou2, iWhichWay)
         call WsloPrmZ(uGiven, MuGiven, 5.0, dWsduou5, iWhichWay)
         dWsduou = (dWsduou5 - dWsduou2)/(5.0-2.0)*(Z-2.0) + dWsduou2   
      else if (Z .le. 10.0) then
         call WsloPrmZ(uGiven, MuGiven, 5.0, dWsduou5, iWhichWay)
         call WsloPrmZ(uGiven, MuGiven, 10.0, dWsduou10, iWhichWay)
         dWsduou = (dWsduou10 - dWsduou5)/(10.0-5.0)*(Z-5.0) + dWsduou5
      else
         call WsloPrmZ(uGiven, MuGiven, 10.0, dWsduou, iWhichWay)
      endif   
      
      END

c
c     -----------------------------------------------------------------|-------
c
      SUBROUTINE WsloPrmZ( uGiven , MuGiven , ZGiven, 
     ^                     dWsduou, iWhichWay)
c     calclate (d_W_s/d_u)/u
c     Ratio of power coupled form the rf source into electromagnetic energy to
c     the rf power absorbed by electrons; P_el / P_in.

      implicit none
c-----input
      real*8 
     & uGiven,  ! v /(\abs{ v_r})
     & MuGiven, ! v_par/v
     & ZGiven   ! effective charge
c-----output
      INTEGER iWhichWay
      Real*8 dWsduou
c-----locals
      REAL*8
     ^        u , Mu , Z , 
     ^        x
      REAL*8    e01, e02, e03, e04, e1, e2, A0,  A1,  A2

c-----Table 4 from the Karney Fisch article p.191
c     Coefficients for approximation to (d_W_s/d_u)/u at (u,mu=1)
c
c     For mu = 1 and 0 < u < 5  the function (d_W_s/d_u)/u 
c     is approximated by 
c     (d_W_s/d_u)/u at (u,mu=1)={Sum_i=1,3 {a_i*u**2i]/
c                                Sum_i=0,3 {b_i*u**2i]}
c     b_0=1
c     Maximum relative error is 5 percents

      REAL*8
     ^    WPpZ01a1, WPpZ01a2, WPpZ01a3, WPpZ01b1, WPpZ01b2, WPpZ01b3 ,
     ^    WPpZ02a1, WPpZ02a2, WPpZ02a3, WPpZ02b1, WPpZ02b2, WPpZ02b3 ,
     ^    WPpZ05a1, WPpZ05a2, WPpZ05a3, WPpZ05b1, WPpZ05b2, WPpZ05b3 ,
     ^    WPpZ10a1, WPpZ10a2, WPpZ10a3, WPpZ10b1, WPpZ10b2, WPpZ10b3 
      DATA
     ^    WPpZ01a1, WPpZ01a2, WPpZ01a3, WPpZ01b1, WPpZ01b2, WPpZ01b3 /
     ^     0.66445, -0.36032,  0.07328,  0.17769, -0.25452,  0.07278 /
      DATA
     ^    WPpZ02a1, WPpZ02a2, WPpZ02a3, WPpZ02b1, WPpZ02b2, WPpZ02b3 /
     ^     0.56760, -0.38984,  0.08634, -0.04019, -0.24673,  0.08508 /
      DATA
     ^    WPpZ05a1, WPpZ05a2, WPpZ05a3, WPpZ05b1, WPpZ05b2, WPpZ05b3 /
     ^     0.39906, -0.32879,  0.07670, -0.28281, -0.16275,  0.07436 /
      DATA
     ^    WPpZ10a1, WPpZ10a2, WPpZ10a3, WPpZ10b1, WPpZ10b2, WPpZ10b3 /
     ^     0.27028, -0.23261,  0.05272, -0.39140, -0.07526,  0.04981 /

c-----Table 5 from the Karney Fisch article p.192
c     Coefficients for approximation to (d_W_s/d_u)/u at (u,mu= - 1)
c
c     For mu = - 1 and 0 < u < 1  the function (d_W_s/d_u)/u 
c     is approximated by 
c     (d_W_s/d_u)/u at (u,mu=-1)=Sum_i=1,4 {a_i*u**2i]
c     b_0=1
c     Maximum relative error is 3 percents

      REAL*8
     ^    WPmZ01a1, WPmZ01a2, WPmZ01a3, WPmZ01a4 ,
     ^    WPmZ02a1, WPmZ02a2, WPmZ02a3, WPmZ02a4 ,
     ^    WPmZ05a1, WPmZ05a2, WPmZ05a3, WPmZ05a4 ,
     ^    WPmZ10a1, WPmZ10a2, WPmZ10a3, WPmZ10a4 

      DATA
     ^    WPmZ01a1, WPmZ01a2, WPmZ01a3, WPmZ01a4 /
     ^    -0.63673, -1.39960,  3.37662, -4.23684 /
      DATA
     ^    WPmZ02a1, WPmZ02a2, WPmZ02a3, WPmZ02a4 /
     ^    -0.55777, -0.80763,  1.43144, -2.03866 /
      DATA
     ^    WPmZ05a1, WPmZ05a2, WPmZ05a3, WPmZ05a4 /
     ^    -0.39704, -0.33811,  0.23607, -0.51011 / 
      DATA
     ^    WPmZ10a1, WPmZ10a2, WPmZ10a3, WPmZ10a4 /
     ^    -0.26600, -0.17342,  0.01896, -0.13349 /
c
c

c      write(*,*)'input WsloPrmZ, uGiven,MuGiven,ZGiven',
c     & uGiven,MuGiven,ZGiven

      u    = abs(uGiven)
      Mu = + MuGiven
      Z    = ZGiven
      iWhichWay = 0

c      write(*,*)'WsloPrmZ, u,Mu,Z',u,Mu,Z

c
c     .....
c
c     For Mu = 1 and 0 < u < 5 use Table IV. on page 191.
      if ( Mu .ge. 0.90 ) then

             if ( u .gt. 5.00 ) then
c     .                                 If the velocity is too large, limit it
c     .                                 to a value covered by the table, and
c     .                                 report WhichWay the velocity is large.
               u  = 5.00
               iWhichWay = +1
             endif
c
             x = u*u
c
             if ( Z  .le. 1.5  )  then
                dWsduou =
     ^           (   WPpZ01a1*x    + WPpZ01a2*x**2 + WPpZ01a3*x**3 ) /
     ^       ( 1. +  WPpZ01b1*x    + WPpZ01b2*x**2 + WPpZ01b3*x**3 )

        else if ( Z  .le. 3.0  )  then
                dWsduou =
     ^           (   WPpZ02a1*x    + WPpZ02a2*x**2 + WPpZ02a3*x**3 ) /
     ^       ( 1. +  WPpZ02b1*x    + WPpZ02b2*x**2 + WPpZ02b3*x**3 )
  
        else if ( Z  .le. 7.0  )  then
                dWsduou =
     ^           (   WPpZ05a1*x    + WPpZ05a2*x**2 + WPpZ05a3*x**3 ) /
     ^       ( 1. +  WPpZ05b1*x    + WPpZ05b2*x**2 + WPpZ05b3*x**3 )

        else
                dWsduou =
     ^           (   WPpZ10a1*x    + WPpZ10a2*x**2 + WPpZ10a3*x**3 ) /
     ^       ( 1. +  WPpZ10b1*x    + WPpZ10b2*x**2 + WPpZ10b3*x**3 )

        endif
      return
      endif

c     For Mu = -1 and 0 < u < 1   use Table V. on page 192.
      if ( Mu .le.-0.90 .and.  u .le. 1.00 ) then
c
             x = u*u
c
             if ( Z  .le. 1.5  )  then
                dWsduou = WPmZ01a1*x 
     ^                 + WPmZ01a2*x**2 + WPmZ01a3*x**3 + WPmZ01a4*x**4

        else if ( Z  .le. 3.0  )  then
                dWsduou = WPmZ02a1*x 
     ^                 + WPmZ02a2*x**2 + WPmZ02a3*x**3 + WPmZ02a4*x**4

        else if ( Z  .le. 7.0  )  then
                dWsduou = WPmZ05a1*x 
     ^                 + WPmZ05a2*x**2 + WPmZ05a3*x**3 + WPmZ05a4*x**4

        else
                dWsduou = WPmZ10a1*x 
     ^                 + WPmZ10a2*x**2 + WPmZ10a3*x**3 + WPmZ10a4*x**4

        endif
      return
      endif


c     For Mu = -1 and 0 < u < 1   use Table V. on page 192.
      if ( Mu .le.-0.90 .and.  u .gt. 1.00 ) then
c     .                                 If the velocity is too large,
c     .                                 damp the results at high u in such
c     .                                 a way that the value and derivative
c     .                                 are continuous at u=e01 (1), but the
c     .                                 values returned at u > e01 (1) are not
c     .                                 too large.  The table is not valid.
c     .                                 The hope is to induce E field
c     .                                 in TSC.
c

c     e01 etc:   power 01 etc of the u^2 (E-like) at which we expand
c     e1, e2:    first and second power of the expansion parameter u^2-e01
c     A0,01,02:  expansion coeficients around u^2=e01
             x = u*u
             iWhichWay = -1
             e01 =1.
             e02=e01*e01
             e03=e01*e02
             e04=e01*e03
             e1 =(x-e01)
             e2 =(x-e01)**2
c
             if ( Z  .le. 1.5  )  then
                A0 =    WPmZ01a1*e01 +    WPmZ01a2*e02
     ^             +    WPmZ01a3*e03 +    WPmZ01a4*e04
                A1 =    WPmZ01a1     + 2.*WPmZ01a2*e01
     ^             + 3.*WPmZ01a3*e02 + 4.*WPmZ01a4*e03
                A2 =                      WPmZ01a2
     ^             + 3.*WPmZ01a3*e01 + 6.*WPmZ01a4*e02

        else if ( Z  .le. 3.0  )  then
                A0 =    WPmZ02a1*e01 +    WPmZ02a2*e02
     ^             +    WPmZ02a3*e03 +    WPmZ02a4*e04
                A1 =    WPmZ02a1     + 2.*WPmZ02a2*e01
     ^             + 3.*WPmZ02a3*e02 + 4.*WPmZ02a4*e03
                A2 =                      WPmZ02a2
     ^             + 3.*WPmZ02a3*e01 + 6.*WPmZ02a4*e02

        else if ( Z  .le. 7.0  )  then
                A0 =    WPmZ05a1*e01 +    WPmZ05a2*e02
     ^             +    WPmZ05a3*e03 +    WPmZ05a4*e04
                A1 =    WPmZ05a1     + 2.*WPmZ05a2*e01
     ^             + 3.*WPmZ05a3*e02 + 4.*WPmZ05a4*e03
                A2 =                      WPmZ05a2
     ^             + 3.*WPmZ05a3*e01 + 6.*WPmZ05a4*e02

        else
                A0 =    WPmZ10a1*e01 +    WPmZ10a2*e02
     ^             +    WPmZ10a3*e03 +    WPmZ10a4*e04
                A1 =    WPmZ10a1     + 2.*WPmZ10a2*e01
     ^             + 3.*WPmZ10a3*e02 + 4.*WPmZ10a4*e03
                A2 =                      WPmZ10a2
     ^             + 3.*WPmZ10a3*e01 + 6.*WPmZ10a4*e02


        endif
        dWsduou = A0 + A1*e1 + A2*e2
      return
      endif


c     Unanticipated parameter range, set return to 0.00.
      dWsduou = 0.00
      return
      END
c                                                                      |
c     Fits.F ends       -----------------------------------------------|









c-------------------------------------------------------
c       normalization in gr
c     cvac - light velocity in (cm/sec)
c      cvac=2.997930d+10
c




c     ve=(sqrt(Te/me)),    ve in (cm/sec),Te(keV)     !article p.1935
c      ve=1.32d9*dsqrt(temp)	                      !temp keV
c     efficiency  in (A/m**2)/(W/m**3) 
c       [then, converted to (A/cm**2)/(erg/(sec*cm**3))]
c     temperature temp in kev
c     density     dens in 10**19 /m**3
c     cln -Coulomb Ln
c      arg1 = 1.d3/temp*dsqrt(10.d0*den)	!temp KeV, den 10**13/cm**3
c      cln=24.d0-dlog(arg1)
c      eta=cprof*amprof*eta0*rprof 		      !article (13)
c      effKarn=eta*3.84d0*temp/(cln*den)	!(A/m**2)/(W/m**3)
c      write(*,*)'effKarn   (A/m**2)/(W/m**3)', effKarn
c-----------------------------------------------------------
c     efficiency  in (A/cm**2)/(erg/(sec*cm**3))
c     temperature temp in kev
c     density     den in 10**19/m**3
c-----------------------------------------------------------
ctest beg
c     LH wave
c      efficien=2.d0/(5.d0+z_eff)*
c     1         (u1*u1+(7.0d0/4.0d0+9.0d0/4.0d0/(3.d0+z_eff)))
c      efficien=efficien*temp/den*(17.d0/cln)*9.d0*1.d-6
c      write(*,*)'u1',u1  
c     write(*,*)'asimptotic efficiency A/cm^2/(egr/(sec*cm**3))',efficien
c      efficien=8.d0/(5.d0+z_eff)*u1*u1 !test
c      write(*,*)'test efficien ',efficien

c      efficien=efficien*temp/den*(17.d0/cln)*9.d0*1.d-6 
c     write(*,*)'efficiency A/cm^2/(egr/(sec*cm**3))',efficien
c      efficien=efficien*3.84d0*temp/den/cln
c      write(*,*)'asymptotic efficiency(A/m**2)/(W/m**3))',efficien
ctest end         determination the of the current drive sign
c      if (u1.ge.0.d0) then
c         s=1.d0
c      else
c         s=-1.d0
c      endif
c      efficien=efficien*s
cSm050923
c      effKarn=effKarn*s*1.d-5  !  A/cm**2/(erg/(sec*cm**3)
c      effKarn=-effKarn*s*1.d-5  !  A/cm**2/(erg/(sec*cm**3)

      subroutine calc_d_eps_r_d_p_dev_d_eps_r_d_omega(z,r,phi,cnpar,
     &cnper,cnz,cnr,cm,d_eps_r_d_p_dev_d_eps_r_d_omega)
c-----calculate
c     calc_d_eps_d_p_dev_d_eps_d_omegaf=(d_eps_r/d_p)/[d_eps_r/d_omega]
c-----uses------------------------------------
c     cold plasma dielectric tensor for LH wave case
c     with the thermal correctron for hermitian part
c     and imaginary eps_33_i for maxwellian f_e
c
c     D.W.Ignat, E.J.Valeo, S.C. Jardin,
c     Dynanmic modelling of lower hybrid current drive,
c     Nuclear Fusion, Vol 34, No.6(1994) p.837-852
c---------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'lsc_approach.i'

c-----input
      real*8
     &z,r,phi,            !space coordinates
     &cnper,              !perpendicular refructive index  
     &cnpar,              !parallel      refructive index
     &cnz,cnr,  
     &cm                  !n_phi=cm/r  

c-----output
      real*8 d_eps_r_d_p_dev_d_eps_r_d_omega

      real*8
     &S,sum_S,sum_alpha,D,P,alpha,
     &x_e,x_i,y_e,y_i,v_t_e,v_t_i,
     &clight,T_e_kev,rme,T_i_kev,
     &d_eps_r_d_omega,
     & d_eps_r_d_p,
     &d_alpha_d_omega,d_P_d_omega,d_D_d_omega,d_S_d_omega,
     &Ds,cnpers,cnper4,cnper6,cnpars

      integer i

c-----externals
      real*8
     &x,y,b,temperho  

c     write(*,*)'in calc_d_eps_r_d_p_dev_d_eps_r_d_omega'

      clight=2.99792458d10              ! light speed [cm.sec

      rme=9.1094d-28                    ! electron mass [g]

c-----calculate small radius rho inside b,
c     rho will be in common one.i      
      bmod=b(z,r,phi)

      T_e_kev=temperho(rho,1)
      v_t_e=dsqrt(T_e_kev*1.6022d-9/rme) !(sqrt(T_e/m_e) electron thermal 
                                         !velocity [cm/sec]

c-----calculate sums S and alpha nondimensional
c
c     S=1+Xe/Y_e^2-Sum{i=2,nbuklk}(X_i)  (8 article)
c
c     alpha=(3.d0/4.d0)*(Xe/Y_e^4)*(v_te/c)^2+3sum{i=2,nbulk}(x_i*(v_ti/c)^2) 
c
c     aplha_article dimensional=alpha/k_0^2      (9 article)
c     k_0=omega/c
c----------------------------------------------------------
      sum_S=0.d0
      sum_alpha=0.d0

      x_e=x(z,r,phi,1)
      y_e=y(z,r,phi,1)

      do i=2,nbulk !sum by ions
         x_i=x(z,r,phi,i)
         y_i=y(z,r,phi,i)
         T_i_kev=temperho(rho,i) 
         v_t_i=dsqrt(T_i_kev*1.6022d-9/(rme*dmas(i)))  ! sqrt(T_i/m_i) ion
                                                       ! thermal
                                                       ! velocity [cm/sec]
         sum_S=sum_S+x_i                               ! nondimensional
                                                       ! (8 article)

         sum_alpha=sum_alpha+x_i*(v_t_i/clight)**2     ! nondimensional                                                                ! =apha_article*
                                                       ! k_0^2 (9 article)
      enddo !i

      S=1.d0+x_e/y_e**2-sum_S                            ! nondimensional
                                                         ! (8 article)
      alpha=(3.d0/4.d0)*x_e/y_e**4*(v_t_e/clight)**2 +   ! nondimensional
     &      3.d0*sum_alpha                               ! =article_apha*
                                                         ! (omega/clight)**2 
                                                         ! (9 article) 

      D=x_e/y_e                                          ! (6 article)
      Ds=D*D

      P=1.d0-x_e                                         ! (10 article)

      cnpars=cnpar*cnpar
      cnpers=cnper*cnper
      cnper4=cnpers*cnpers
      cnper6=cnper4*cnpers

c-----real part of te dispersion funcrion:
c     eps_r=-alpha*k_perp**6 + k_perp**4*S+
c           +k_perp**2*[(P+S)(k_par**2-k_0**2*S)+k0**2*D**2]+
c           +P*[(k_par**2-k_0**2*S)**2-k0**4*D**2]
c-----------------------------------------------------------------
c     d_eps_r_d_p=d_eps_r_d_p_dimensional/k0^4
c     d_eps_r_d_omega=d_eps_r_d_omega_dimensional*omega/k_0^4
c     Here 
c     d_eps_r_d_omega is nondimensional
c     d_eps_r_d_p is nondimensional
c----------------------------------------------------------------

      d_eps_r_d_p=cnpers*(cnpars-S)+((cnpars-S)**2-Ds)  !nondimensional (5)

      d_S_d_omega=2.d0*sum_S                            !nondimensional (13)
      d_alpha_d_omega=-12.d0*sum_alpha                  !nondimensional (14)
      d_P_d_omega=2.d0*x_e                              !nondimensional (15)
      d_D_d_omega=-x_e/y_e                              !nondimensional (16)

       d_eps_r_d_omega=-cnper6*d_alpha_d_omega
     &   +cnper4*d_S_d_omega
     &   +cnpers*((d_P_d_omega+d_S_d_omega)*(cnpars-S)-
     &          (P+S)*(2.d0*S+d_S_d_omega)+(2.d0*Ds+2.d0*D*d_D_d_omega))
     &   +d_P_d_omega*((cnpars-S)**2-Ds)    
     &   +P*(-2.d0*(cnpars-S)*(2.d0*S+d_S_d_omega)-
     &        (4.d0*Ds+2.d0*D*d_D_d_omega))

c       write(*,*)'cnper6,d_alpha_d_omega',cnper6*d_alpha_d_omega
c       write(*,*)'cnper4,d_S_d_omega',cnper4,d_S_d_omega
c       write(*,*)'d_P_d_omega,cnpars,S,Ds',d_P_d_omega,cnpars,S,Ds
c       write(*,*)'cnpers',cnpers
c       write(*,*)'P,S,Ds,2.d0*D*d_D_d_omega',P,S,Ds,2.d0*D*d_D_d_omega

c       write(*,*)'1 -cnper6*d_alpha_d_omega',-cnper6*d_alpha_d_omega
c       write(*,*)'2 cnper4*d_S_d_omega',cnper4*d_S_d_omega
c       write(*,*)'3 d_P_d_omega*((cnpars-S)**2-Ds)',
c     &            d_P_d_omega*((cnpars-S)**2-Ds)
c       write(*,*)'4 (d_P_d_omega+d_S_d_omega)*(cnpars-S)',
c     &              (d_P_d_omega+d_S_d_omega)*(cnpars-S)
c       write(*,*)'5 (P+S)*(2.d0*S+d_S_d_omega)+
c     &              (2.d0*Ds+2.d0*D*d_D_d_omega)',
c     &              (P+S)*(2.d0*S+d_S_d_omega)+
c     &              (2.d0*Ds+2.d0*D*d_D_d_omega)
c       write(*,*)'6 P*(-2.d0*(cnpars-S)*(2.d0*S+d_S_d_omega)-
c     &        (4.d0*Ds+2.d0*D*d_D_d_omega))',
c     &              P*(-2.d0*(cnpars-S)*(2.d0*S+d_S_d_omega)-
c     &        (4.d0*Ds+2.d0*D*d_D_d_omega))
c       write(*,*)'6_1 P*(-2.d0*(cnpars-S)*(2.d0*S+d_S_d_omega))',
c     &                P*(-2.d0*(cnpars-S)*(2.d0*S+d_S_d_omega))
c       write(*,*)'6_2 (4.d0*Ds+2.d0*D*d_D_d_omega)',
c     &                (4.d0*Ds+2.d0*D*d_D_d_omega)
c-----result nondimensional 
      d_eps_r_d_p_dev_d_eps_r_d_omega=d_eps_r_d_p/d_eps_r_d_omega

c      write(*,*)'d_eps_r_d_p,d_eps_r_d_omega',
c     &           d_eps_r_d_p,d_eps_r_d_omega
c      write(*,*)'d_eps_r_d_p_dev_d_eps_r_d_omega',
c     &           d_eps_r_d_p_dev_d_eps_r_d_omega

c      stop ' d_eps_r_d_p_dev_d_eps_r_d_omega'
      return
      end


      subroutine calc_K_zz_i(x_e,v_par_res,K_zz_i,
     &sigma_lsc,rho)
c-----calculates K_zz_i using integral from delta function
      implicit none
      include 'param.i'
      include 'lsc_approach.i'
 
c-----input
      real*8
     &x_e,            !(omega_pe/omega)**2
     &v_par_res,      !resonance velocity normalized by v_te
     &sigma_lsc,      !width of the delta function
     &rho
      
c-----externals 
      real*8  delta_function_lsc
c-----output
      real*8 K_zz_i

c-----K_zz_i=-pi*x_e*d_f_e_d_v(v_par_res)*v_par_res*dabs(v_par_res)

c-----locals
      integer nv
      real*8 x,pi,h_v,hrho
      integer i_rho !number of radial point !1=<i_rho=< n_psi_TSC
cfor test
c      real*8 rholoc,zeffloc,deriv
c      real*8 zeffrho,d_cl !externals
cend for test

      if(dabs(v_par_res).gt. v_par_lsc_dev_vte_max) then
        K_zz_i=0.d0
        return
      endif

      pi=4.0*datan(1.d0)

      hrho=1.d0/dble(n_psi_TSC)
      i_rho=rho/hrho+1
c      write(*,*)' in K_zz_i rho,hrho,i_rho',rho,hrho,i_rho

      if(i_rho.gt.n_psi_TSC) then
         write(*,*)'WARNING: lsc_approach.f in lcalc_K_zz_i'
         write(*,*)'i_rho>n_psi_TSC'
         write(*,*)'LSC works at rho.le.1, but rho=',rho
         write(*,*)'it was set i_rho=n_psi_TSC'
         i_rho=n_psi_TSC
      endif

cfor  test
c      rholoc=hrho*(i_rho-0.5d0)
c      zeffloc=zeffrho(rholoc)    
cend  for test

      h_v=v_par_mesh_lsc(1)-v_par_mesh_lsc(0)
      K_zz_i=0.d0
      do nv = -nv_lsc_ql_lh,nv_lsc_ql_lh
         x=v_par_res-v_par_mesh_lsc(nv)
         K_zz_i= K_zz_i+ 
     &           v_par_mesh_lsc(nv)*d_fe_dv_lsc(i_rho,nv)*
     &           h_v*delta_function_lsc(x,h_v,sigma_lsc)

cSAP100519
c        if (delta_function_lsc(x,h_v,sigma_lsc).gt.1.d-1) then
c         write(*,*)'in K_zz nv,d_fe_dv_lsc(i_rho,nv),fe_lsc(i_rho,nv)',
c     &                     nv,d_fe_dv_lsc(i_rho,nv),fe_lsc(i_rho,nv)
c         write(*,*)'d_ql_lsc_lh_ar((i_rho,nv)',d_ql_lsc_lh_ar(i_rho,nv)
c         write(*,*)'d_cl(v_par_mesh_lsc(nv),zeffloc)',
c     &              d_cl(v_par_mesh_lsc(nv),zeffloc)
c         write(*,*)'d_ql_lsc_lh_ar(i_rho,nv)/
c     &    d_cl(v_par_mesh_lsc(nv),zeffloc)',
c     &    d_ql_lsc_lh_ar(i_rho,nv)/d_cl(v_par_mesh_lsc(nv),zeffloc)
c         deriv=v_par_mesh_lsc(nv)*d_cl(v_par_mesh_lsc(nv),zeffloc)/
c     &   (d_cl(v_par_mesh_lsc(nv),zeffloc)+d_ql_lsc_lh_ar(i_rho,nv))
c         write(*,*)' deriv',deriv
c        endif

      enddo
c      write(*,*)'K_zz_i,x_e,v_par_res',K_zz_i,x_e,v_par_res
      K_zz_i=-K_zz_i*pi*x_e*dabs(v_par_res)
c      write(*,*)'end K_zz_i',K_zz_i
      return
      end
          
      subroutine calc_fe_lsc_arrays
c------------------------------------------------------------------
c     calculate 
c     array of electron dictribution function and it's derivative
c
cSAP100122
c     fe_lsc(1:n_psi_TSC-1,-nv_lsc_ql_lh : nv_lsc_ql_lh)
c     d_fe_dv_lsc(1:Nn_psi_TSC,-nv_lsc_ql_lh : nv_lsc_ql_lh)
c   
c     fe is normaliazed at unit
c
c     integral[v_par_mesh_ls(k,-nv_lsc_ql_lh),v_par_mesh_ls(k,nv_lsc_ql_lh)}
c             [f_e*d_v_par]=1
c    
c-------------------------------------------------------------------
      implicit none
c      include 'onetwo.i'
      include 'param.i'
      include 'lsc_approach.i'

c-----output
c
cSAP100122
c     fe_lsc(1:n_psi_TSC,-nv_lsc_ql_lh : nv_lsc_ql_lh) 
c     d_fe_dv_lsc(1:n_psi_TSC,-nv_lsc_ql_lh : nv_lsc_ql_lh)
c     these arrays will be in lsc_approch_no_nml.i

c-----externals
      real*8 
     &zeffrho,d_cl
c-----local
      real*8 h_v,h_rho,zeff,v_par_p,v_par_m,rho,norma,coef,v_par
              
      integer n,k

      h_v=v_par_lsc_dev_vte_max/(nv_lsc_ql_lh)

      h_rho=1.d0/n_psi_TSC


c-----calculate integral-----------------------------------
c     integral{0,v_par}[d_cl/(d_cl+d_ql)]
c-----initialization--------------------------------------
      integral_fe_lsc=0.d0
c----------------------------------------------------------   
      do k=1,n_psi_TSC
         rho=h_rho*(k-0.5d0)
         zeff=zeffrho(rho)
    
         do n=1,nv_lsc_ql_lh
            v_par_p=v_par_mesh_lsc(n)
            v_par_m=v_par_mesh_lsc(-n)       

            integral_fe_lsc(k,n)=integral_fe_lsc(k,n-1)+h_v*v_par_p*
     &           d_cl(v_par_p,zeff)/
     &           (d_cl(v_par_p,zeff)+d_ql_lsc_lh_ar(k,n)) 

            integral_fe_lsc(k,-n)=integral_fe_lsc(k,-(n-1))-h_v*v_par_m*
     &           d_cl(v_par_m,zeff)/
     &           (d_cl(v_par_m,zeff)+d_ql_lsc_lh_ar(k,-n))          
             
         enddo
      enddo

c-----calculate array of electron disribution function
      do k=1,n_psi_TSC
         do n=-nv_lsc_ql_lh,nv_lsc_ql_lh           
            fe_lsc(k,n)=dexp(-integral_fe_lsc(k,n))  
         enddo
      enddo

c-----normalization of the electron disribution function
      do k=1,n_psi_TSC 
         norma=0.d0       
         do n=-nv_lsc_ql_lh,nv_lsc_ql_lh  
            norma=norma+fe_lsc(k,n)*h_v  
         enddo
      
         do n=-nv_lsc_ql_lh,nv_lsc_ql_lh
            fe_lsc(k,n)=fe_lsc(k,n)/norma         
         enddo
      enddo


c-----calculate array of derivative df_e_lsc_d_v_par
      do k=1,n_psi_TSC
         do n=-nv_lsc_ql_lh,nv_lsc_ql_lh
            v_par=v_par_mesh_lsc(n)
 
            coef=v_par*d_cl(v_par,zeff)/
     &          (d_cl(v_par,zeff)+d_ql_lsc_lh_ar(k,n))  

            d_fe_dv_lsc(k,n)=fe_lsc(k,n)*(-coef)
         enddo
       enddo

      return
      end

      real*8 function delta_function_lsc(x,h_v,sigma_lsc)
c-----approximation of the delta function in LSC approach
c
c     Integral{delta_function_lsc*(x)*dx}=1
c-------------------------------------------------------
      implicit none 
      include 'param.i'
      include 'lsc_approach.i'
c-----input
      real*8 x, ! argument of delta function
     &h_v,      ! step of velocity mesh
     &sigma_lsc ! parameter to smooth delta_function

c-----locals
      integer j
      real*8 pi

      if(dabs(x).gt. v_par_lsc_dev_vte_max) then
        delta_function_lsc=0.d0
        return
      endif

      goto 10

      j=sigma_lsc/h_v
c      write(*,*)'in delta function sigma_lsc,h_v,j',sigma_lsc,h_v,j
c      j=0 !for test
c      j=10 !for test
      if (dabs(x).gt.(j+0.5)*h_v) then
         delta_function_lsc=0.d0
      else
         delta_function_lsc=1.d0/(h_v*(2*j+1))
      endif
 10   continue
c-----exp(-(x/(2sigma)**2)
      pi=4.d0*datan(1.d0)
      delta_function_lsc=dexp(-(x/(1.d0*sigma_lsc))**2)/
     &                       (1.d0*sigma_lsc*dsqrt(pi))

      return
      end


      real*8 function d_cl(v_d_v_te,z_eff)
c-----calculates normalized Coulomb difusion coefficient
c
c     D_c=beta_z/[1+(v_d_v_te)^2]^3/2
c
c     beta_z=(1+Z)/5

      implicit none
c-----input
      real*8 
     &v_d_v_te,    ! parallel velocity normalized at v_te
     &z_eff

c-----local
      real*8 beta_z,p
      p=1.d0+v_d_v_te**2

      beta_z=(1.d0+z_eff)/5.d0
 
      d_cl=beta_z/(p*dsqrt(p))

      return
      end


      SUBROUTINE d_ql_lh_prof(is,ir)
c----------------------------------------------------------------
c     this subroutine calculates LH quasilinear diffusion coefficients
c     using power along the ray
c     
c-----------------------------------------------------------------
      implicit none
      include 'param.i'     
      include 'one.i' 
      INCLUDE 'rho.i'
      include 'lsc_approach.i'      
      INCLUDE 'onetwo.i'
      INCLUDE 'write.i'
      include 'writencdf.i'
      
c-----input
  
      integer
     & is,   !ray element
     & ir    !ray number

c     following data are in lsc_approach_nml.i
c     integer nv_lsc_ql_lh   !number of velocity mesh points 
c     real*8  v_dev_vte_max  !maximal normilazed velosity of the velocity mesh

c-----output
c      
c     d_ql_lsc_lh_ar(n_psi_TSC,-nv_lsc_ql_lh,nv_lsc_ql_lh
c     v_par_mesh_lsc(-nv_lsc_ql_lh,nv_lsc_ql_lh),

c-----locals      
      real*8
     &mass_e,charge_e,clight,k_1_kev,
     &hrho,h_v,
     &temp_kev,dense,v_te,arg1,cln,tau_n,
     &rholeft,rhoright,rhomax,rhomin,del_s_pol,delrho,
     &r_m,z_m,
     &d_ql_lsc_coef_1,sigma_lsc,del_d_ql,v_res,
     &t,u(6),deru(6),vgr(3),vgr_pol,
     &rhobegin,rhoend,cnpar,cnper,
     &epsilon_0,  ! for test
     &masse_kg,   ! electron mass [kg]
     &charge_e_c, ! electron charge (coulomb)
     &k_si_1_kev,  
     &d_eps_r_d_p_dev_d_eps_r_d_omega,
     &z_l,r_l,phi_l,cnpar_l,cnper_l,
     &cnz_l,cnr_l,cm_l,delta_func_value,
     &dense_m3,v_te_m_sec,coef_watt_m3,
     &power_dens_watt_m3_loc,power_dens_erg_cm3_sec_loc,
     &ratio,delta_pow_loc,hv,delta_pow_is,delta_pow_is_m_1,
     &p_abs,
     &lsc_eps_coef,deltaW,  !for test only 
     &zeff_loc,d_cl_loc, 
     &del_s_pol_bin_test,    !for test only
     &delta_power_prep3d,    !for test only
     &K_zz_i,x_e,r00

      complex*16
     &eps(3,3)      !for teest only   
      integer 
     &i,j,jbinmin,jbinmax,nv,i_geom_optic_loc,
     &k
      integer jjkk
      logical first
              
c-----external
      real*8
     &temperho,densrho,rhov,rhos,delta_function_lsc,b,
     &zeffrho,d_cl,x

      data first / .true./
      save first
      save delta_pow_is_m_1

      mass_e=9.1094d-28    ! electron mass       [g]
      charge_e=4.8032d-10  ! electron charge     [statcoulomb]
      clight=2.99792458d10 !                     [cm/sec]
      k_1_kev=1.6022d-9    ! egrs in 1 KeV       [erg]

c-----for test in SI
      masse_kg=9.1094d-31      ! electron mass   [kg]
      charge_e_c=1.6022d-19    ! electron charge [coulomb]
      k_si_1_kev=1.6022d-16    ! Joules in 1 KeV [Joule]
      epsilon_0=1.d0/(4.d0*pi*9.d9)

      pi=4.d0*datan(1.d0)

      hrho=1.d0/dble(n_psi_TSC)                 !radial TSC mesh step 

      h_v=v_par_lsc_dev_vte_max/(nv_lsc_ql_lh)  ! velosity mesh step 
          
c      write(*,*)'inSUBROUTINE d_ql_lh_prof is,ir',is,ir

c-----coefficient for QL difusion calculations 
c     8*pi^2*(e/m_e)^2/(clight*omega): 
c     [(g*cm^3/sec^2) /g^2 /(cm/sec^2]= [cm^2/g]
c
c     e^2: [g*cm^3/sec^2]
c     
c-----test in SI
      d_ql_lsc_coef_1=2.d0*(pi/epsilon_0)*(charge_e_c/masse_kg)**2/
     &                (clight*1.d-2*2*pi*frqncy*1.d9) ![m^2/kg]

c      write(*,*)'SI d_ql_lsc_coef_1',d_ql_lsc_coef_1

c-----in Gaussian Units
      d_ql_lsc_coef_1=8.d0*pi*pi*(charge_e/mass_e)**2/
     &                (clight*2*pi*frqncy*1.d9) ![cm^2/g] 
                                                !Here frqncy is in [GHZ]
c-----Here:
c     clight [cm/sec] is from normalized Vgr_pol
c     omega=2*pi*frqncy*1.d9 [HZ] is from delta function
c-------------------------------------------------------
c     now proposed that r0x=1 m
      r00=100.d0*r0x

      if (first) THEN
c------- create velocity mesh
         do nv=-nv_lsc_ql_lh ,nv_lsc_ql_lh 
            v_par_mesh_lsc(nv)=h_v*nv !normalized at v_te
                                      !v_te is different 
                                      !at different radial points    
         enddo
c------- create arrays:
c        v_te_bin_center_lsc(i), 
c        tau_n_bin_center_lsc(i), i=1,n_psi_TSC

         hrho=1.d0/n_psi_TSC
         do i=1,n_psi_TSC+1
            rho_bin_lsc(i)=hrho*(i-1)
         enddo
  
         do i=1,n_psi_TSC 

            rho_bin_center_lsc(i)=
     &           0.5d0*(rho_bin_lsc(i)+rho_bin_lsc(i+1))
  
            binvol_lsc(i)=voltot*(rhov(rho_bin_lsc(i+1))**2-
     &                            rhov(rho_bin_lsc(i))**2)*1.d6 ![cm**3]

            binarea_lsc(i)=areatot*(rhos(rho_bin_lsc(i+1))**2-
     &                              rhos(rho_bin_lsc(i))**2)*1.d4 ![cm**2]
	    temp_kev=temperho(rho_bin_center_lsc(i),1)

            dense=densrho(rho_bin_center_lsc(i),1)   ! 10**13 /cm**3

            v_te=dsqrt(temp_kev*1.6022d-9/mass_e)    ! (sqrt(T_e/m_e) electron
                                                     ! thermal velocity [cm/sec]
            v_te_bin_center_lsc(i)=v_te               !cm/sec
           
c-----------cln=Coulomb Ln   
            arg1 = 1.d3/temp_kev*dsqrt(10.d0*dense)  ! temp [KeV],
                                                     ! dense [10**13/cm**3]
            cln=24.d0-dlog(arg1)
c----------------------------------------------------------------------
c           tau_n character time used in Karney code
c------------------------------------------------------------------
c           tau_n=4*pi*epsolin_0**2*T**3/2*m**1/2/(e**4*dens*Cln) [SI]=
c              =T**3/2*m**1/2/(4*pi*e**4*dens*Cln)  [cgs]
c-------------------------------------------------------------------
            tau_n=1.d0/(4.d0*pi)*
     &      (k_1_kev*temp_kev)**1.5d0*dsqrt(mass_e)/
     &      (charge_e**4*dense*1.d13*cln) ! sec

            tau_n_bin_center_lsc(i)= tau_n          
	 enddo
         first = .false.
      endif 

c-----wr_nc and wz_nc are in (cm)

c-----cordinates of the end point of the ray element
      r_m=wr_nc(is,ir)*0.01d0/r0x ! normalization
      z_m=wz_nc(is,ir)*0.01d0/r0x ! normalization
      bmod=b(z_m,r_m,0.d0)
      rhoend=rho

c-----coordinates of the initial point of the ray element
      if (is.gt.1) then
        r_m=wr_nc(is-1,ir)*0.01d0/r0x ! normalization
        z_m=wz_nc(is-1,ir)*0.01d0/r0x ! normalization
        bmod=b(z_m,r_m,0.d0)
        rhobegin=rho
      endif

      cnpar=wnpar_nc(is-1,ir)
      cnper=wnper_nc(is-1,ir) 

c  radius rho is in common/one/,it was calculated in subroutine: b
c  rho is used in functions:tempe,dense,z_eff
c  The magnetic field bmode is in common/one/ .bmode is used
c  in function y(z,r,phi)

c      write(*,*)'ir,is,rhobegin,rhoend',ir,is,rhobegin,rhoend

      if (is.eq.1) then
c--------calculations on the first step                
c        nothing to do
	 goto 130
      endif

c-----------------------------------------------------------
c     index j of TSC mesh point is a number of central bin point
c
c     Calculate minimal and maximal indexes (jbinmin,jbinmax)
c     of radial bins intersected by the ray elements (is-1,is)  
c
c     rhobegin is small radius of ray point  (is-1)
c     rhoend   is small radius of ray point  (is)
c
c---------------------------------------------------------------
      del_s_pol=ws_nc(is,ir)-ws_nc(is-1,ir)
c      goto 100 
      if(rhoend.gt.rhobegin) then
         rhomax=rhoend
         rhomin=rhobegin
      else
	 rhomin=rhoend
	 rhomax=rhobegin
      endif

      if(rhomax.gt.1.d0) rhomax=1.d0
      if(rhomin.gt.1.d0) rhomin=1.d0-1.d-10    

      hrho=1.0d0/n_psi_TSC 
     
      jbinmin=rhomin/hrho+1
      jbinmax=rhomax/hrho+1
     
c       del_s_pol=ws_nc(is,ir)-ws_nc(is-1,ir)
       
      if(jbinmin.gt.n_psi_TSC)jbinmin=n_psi_TSC 
      if(jbinmax.gt.n_psi_TSC)jbinmax=n_psi_TSC 

c      write(*,*)'jbinmin,jbinmax',jbinmin,jbinmax
c--------------------------------------------------------------------
c     del_l_pol (cm)
c-----------------------------------------------------------
            
      if (jbinmin.eq.jbinmax) then       
          del_s_pol_bin_lsc(jbinmin)=del_s_pol          
         goto 30
      endif

      delrho=rhomax-rhomin
      
      del_s_pol_bin_lsc(jbinmin)=del_s_pol*(hrho*jbinmin-rhomin)/delrho
      del_s_pol_bin_lsc(jbinmax)=del_s_pol*(rhomax-hrho*(jbinmax-1))/
     &                           delrho

      if(jbinmax.gt.(jbinmin+1)) then
         do j=(jbinmin+1),(jbinmax-1)
            del_s_pol_bin_lsc(j)=del_s_pol*hrho/delrho
         enddo    
      endif
 100  continue
cSAP100504----to use only one rho bin
      hrho=1.0d0/n_psi_TSC
cAYP aint->idint 
      jbinmin=idint(rhobegin/hrho)+1      ! the number of radial bin
      if(jbinmin.gt.n_psi_TSC) jbinmin=n_psi_TSC 
      jbinmax=jbinmin
 
      del_s_pol_bin_lsc(jbinmin)=del_s_pol
c-------------
cSAP100419 test del_pol_bin-----------------------
c      del_s_pol_bin_test=0.d0
c      do j=jbinmin,jbinmax
c         del_s_pol_bin_test =del_s_pol_bin_test+
c     &                        del_s_pol_bin_lsc(j)                  
c      enddo
c      write(*,*)'is,del_s_pol_bin_test,del_s_pol',
c     &           is,del_s_pol_bin_test,del_s_pol      
cend test ---------------------

 30   continue

      i_geom_optic_loc=i_geom_optic
      i_geom_optic=1                 !to get group velocity in deru 

      t=0.d0

c---- wr_nc and wz_nc are in (cm) 
      u(1)=wz_nc(is-1,ir)*0.01d0/r0x ! normalization
      u(2)=wr_nc(is-1,ir)*0.01d0/r0x ! normalization        
      u(3)=wphi_nc(is-1,ir)
      u(4)=wn_z_nc(is-1,ir)  
      u(5)=wn_r_nc(is-1,ir)     
      u(6)=u(2)*wn_phi_nc(is-1,ir)

      call rside1(t,u,deru)
      i_geom_optic=i_geom_optic_loc

      vgr(1)=deru(1)
      vgr(2)=deru(2)
      !vgr(3)=u(1)*deru(3) !YuP[2020-08-28] BUG? Not important
      vgr(3)=u(2)*deru(3) !YuP[2020-08-28] Correct?  R*(dphi/dt) Not used anyway

      vgr_pol=dsqrt(vgr(1)**2+vgr(2)**2) ! poloidal group velocity
                                         ! normalized to clight    
      z_l= u(1)
      r_l= u(2)
      phi_l= u(3)
      cnpar_l=wnpar_nc(is-1,ir)
      cnper_l=wnper_nc(is-1,ir)
      cnz_l= u(4)
      cnr_l= u(5)
      cm_l= u(6)
          
c------------------------------------------------------------------
c     calculate multiplier: d_eps_r_d_p_dev_d_eps_r_d_omega
c     at 'is-1' ray point
c     It wil be used at QL coefficient calculations
c     d_eps_r_d_p_dev_d_eps_r_d_omega=(d_eps_r/d_p)/[d_eps_r/d_omega]
c------------------------------------------------------------------  
      call  calc_d_eps_r_d_p_dev_d_eps_r_d_omega(z_l,r_l,phi_l,
     &      cnpar_l,cnper_l,cnz_l,cnr_l,cm_l,
     &      d_eps_r_d_p_dev_d_eps_r_d_omega) 

c      write(*,*)'1 rho,rhobegin',rho,rhobegin
c---- for ratio calculations     
      x_e=x(z_l,r_l,phi_l,1)  !for calculation K_xx_i

      delta_pow_is =0.d0 !for test only

c      write(*,*)'in d_ql_lh_prof: is,ir,jbinmin,jbinmax',
c     &                            is,ir,jbinmin,jbinmax  

      do j=jbinmin,jbinmax !j is a number of the radial bin at TSC mesh

         v_res=(clight/v_te_bin_center_lsc(j))/cnpar_l  ! resonance velocity
                                                        ! normalized to v_te(j)
                                                        ! v_te is different
                                                        ! at different
                                                        ! radial points j       
         sigma_lsc=dabs(v_res*wdnpar_nc(is,ir)/cnpar_l)
c-----------------------------------------------------------------    

         call calc_K_zz_i(x_e,v_res,K_zz_i,          !rho 
     &        sigma_lsc,rho_bin_center_lsc(j))       !for j point of TSC mesh                                                         !xe from 'is-1' point
         lsc_eps_coef=d_eps_r_d_p_dev_d_eps_r_d_omega*K_zz_i 
c--------------------------------------------------------------------       
c        d_ql_lsc_lh_ar_loc(-nv_lsc_ql_lh : nv_lsc_ql_lh)
c        is incemental D_ql at radial point "j" at TSC mesh
c        from ray with number "ir" and from ray bin (is-i,is)
c------------------------------------------------------------------ 
c         write(*,*)'in d_ql_lh_prof is,ir,delpwr_nc(is-1,ir)',
c     &                              is,ir,delpwr_nc(is-1,ir)
         d_ql_lsc_lh_ar_loc=0.d0

         do nv=-nv_lsc_ql_lh,nv_lsc_ql_lh

            if(dabs(v_res).gt.v_par_lsc_dev_vte_max) then
               del_d_ql=0.d0
               goto 50
            endif
  
            delta_func_value=delta_function_lsc(
     &                       v_res-v_par_mesh_lsc(nv),h_v,sigma_lsc)

c           write(*,*)'nv,v_par_mesh_lsc(nv),delta_func_value',
c     &                nv,v_par_mesh_lsc(nv),delta_func_value
c           write(*,*)'delpwr_nc(is-1,ir)',delpwr_nc(is-1,ir)
c           write(*,*)'del_s_pol_bin_lsc(j)',del_s_pol_bin_lsc(j)
c           write(*,*)'vgr_pol,binvol_lsc(j)',vgr_pol,binvol_lsc(j)
c           write(*,*)'d_ql_lsc_coef_1',d_ql_lsc_coef_1
c           write(*,*)'d_eps_r_d_p_dev_d_eps_r_d_omega',
c     &                d_eps_r_d_p_dev_d_eps_r_d_omega
c           write(*,*)'v_res',v_res
c           write(*,*)'tau_n_bin_center_lsc(j)',tau_n_bin_center_lsc(j)
c           write(*,*)'v_te_bin_center_lsc(j)',v_te_bin_center_lsc(j)
            del_d_ql=delpwr_nc(is-1,ir)*                   !power in ray chanel
     &      (del_s_pol_bin_lsc(j)/vgr_pol)/binvol_lsc(j)* !deltaTime/deltaVolume
     &      d_ql_lsc_coef_1*                           !8pi^2(e/m)^2/(c*omega)
     &      d_eps_r_d_p_dev_d_eps_r_d_omega*
     &      dabs(v_res)*delta_func_value*      !delta(omega-kv)=|vres|/omega*
                                               !delta(vres-v) ,
                                               !vres, v are normalized to Vte
     &      (tau_n_bin_center_lsc(j)/v_te_bin_center_lsc(j)**2) !normalization

 50         continue
c            write(*,*)'nv,del_d_ql',nv,del_d_ql

            d_ql_lsc_lh_ar_loc(nv)=del_d_ql
          
            if(del_d_ql.lt.0.d0) then
              write(*,*)'lsc_approach d_ql_lh_prof(is,ir)'
              write(*,*)'del_d_ql.lt.0.d0,is,ir,j,nv',is,ir,j,nv
              write(*,*)'delpwr_nc(is-1,ir)',delpwr_nc(is-1,ir)
              write(*,*)'del_s_pol_bin_lsc(j)',del_s_pol_bin_lsc(j)
              write(*,*)'binvol_lsc(j)',binvol_lsc(j)
              write(*,*)'vgr_pol',vgr_pol
              write(*,*)'d_ql_lsc_coef_1',d_ql_lsc_coef_1
            write(*,*)'tau_n_bin_center_lsc(j)',tau_n_bin_center_lsc(j)
              write(*,*)'v_te_bin_center_lsc(j)**2',
     &                   v_te_bin_center_lsc(j)**2
              write(*,*)'v_res',v_res
              write(*,*)'d_eps_r_d_p_dev_d_eps_r_d_omega',
     &                   d_eps_r_d_p_dev_d_eps_r_d_omega
              write(*,*)'delta_func_value',delta_func_value
              write(*,*)'del_d_ql',del_d_ql
              stop 'del_d_ql.lt.0.d0'
            endif

         enddo !nv

c--------calculate power density:
c
c        power_dens_watt_m3_loc and power_dens_erg_cm3_sec_lo
c
c        from incremental d_ql_lsc_lh_ar_loc
c        and using electron distribution d_fe/d_v with old total D_ql

         hv = v_par_mesh_lsc(1)-v_par_mesh_lsc(0)
         power_dens_watt_m3_loc=0.d0
         do nv = -nv_lsc_ql_lh,nv_lsc_ql_lh
            power_dens_watt_m3_loc = power_dens_watt_m3_loc
     &      - v_par_mesh_lsc(nv)*d_ql_lsc_lh_ar_loc(nv)*
     &        d_fe_dv_lsc(j,nv)*hv
         enddo

         dense_m3 = densrho(rho_bin_center_lsc(j),1)*1.d19  !electron density
                                                            ! 1/m**3
         v_te_m_sec=v_te_bin_center_lsc(j)*1.d-2  !thermal velocity in [m/sec]

         coef_watt_m3 = dense_m3*masse_kg*v_te_m_sec**2/
     &                  tau_n_bin_center_lsc(j) !watt/m**3

         power_dens_watt_m3_loc=
     &          power_dens_watt_m3_loc*coef_watt_m3   !watt/m**3

         power_dens_erg_cm3_sec_loc=
     &          power_dens_watt_m3_loc*10.d0   !erg/cm**3/sec

c--------incremental QL power from old total fe from 'is' ray bin in j volume bin [erg]
         delta_pow_loc=power_dens_erg_cm3_sec_loc* binvol_lsc(j)

c--------incremental QL power from old total fe from 'is' ray bin in all volumes  [erg] 
         delta_pow_is=delta_pow_is+delta_pow_loc   
                  
c         write(*,*)'bef ratio calculation: delta_pow_loc into j volume',
c     &              j,delta_pow_loc

c         write(*,*)'(delpwr_nc(is-1,ir)-delpwr_nc(is,ir))*
c     &              del_s_pol_bin_lsc(j)/del_s_pol',
c     &              (delpwr_nc(is-1,ir)-delpwr_nc(is,ir))*
c     &              del_s_pol_bin_lsc(j)/del_s_pol
                
c         write(*,*)'del_s_pol_bin_lsc(j)/del_s_pol',
c     &              del_s_pol_bin_lsc(j)/del_s_pol

c--------calculate for test: delta power at one ray bin like in prep3d_lsc
c        delta_power_prep3d.
c        It should be equal to (delpwr_nc(is-1,ir)-delpwr_nc(is,ir))
c
c        calculate ratio: 
c
c        ratio=[(delpwr_nc(is-1,ir)-delpwr_nc(is,ir))*del_s_pol_bin_lsc(j)/del_s_pol]
c        /delta_pow_loc

CAYP201124 added for error tracing as delta delpwr can be 0
c         write(*,*)'delta_pow_loc/delwpr_nc',
c     &              delta_pow_loc,delpwr_nc(is-1,ir)

CAYP201124 avoid zero powers
       jjkk=0
       if(dabs(delpwr_nc(is-1,ir)).lt.1.d-177) jjkk=jjkk+1
       if(dabs(delta_pow_loc).lt.1.d-177) jjkk=jjkk+1
!YuP       if(jjkk.gt.0) then
!YuP        ratio=1.d0
c         write(*,*)'delta_pow_loc/delwpr_nc',
c     &              delta_pow_loc,delpwr_nc(is-1,ir)
!YuP       else
         if(dabs(delta_pow_loc).gt.1.d-8*dabs(delpwr_nc(is-1,ir))) then
         !YuP[2020-11-25] rearranged if() to avoid divisions.
!YuP         if (dabs(delta_pow_loc/delpwr_nc(is-1,ir)).gt.1.d-8) then
!YuP           if(((delpwr_nc(is-1,ir)-delpwr_nc(is,ir))/
!YuP     &        delpwr_nc(is-1,ir))
!YuP     &.       gt.1.d-8) then
           if( delpwr_nc(is-1,ir)-delpwr_nc(is,ir)
     &         .gt. 1.d-8*delpwr_nc(is-1,ir)) then

c-------------calculation delta_power like in prep3d_lsc for test                            
c             lsc_eps_coef=d_eps_r_d_p_dev_d_eps_r_d_omega*K_zz_i 
c             d_eps_r_d_p_dev_d_eps_r_d_omega from is-i ray point
c             K_zz_i  from j point of TSC mesh, but xe from is-i ray point
c----------------------------------------------------------------------
              delta_power_prep3d=2.d0*lsc_eps_coef*delpwr_nc(is-1,ir)*
     &               2.d0*pi*frqncy*1.d9*
     &               (ws_nc(is,ir)-ws_nc(is-1,ir))/(vgr_pol*clight)
 
c              write(*,*)'delta_power_prep3d',delta_power_prep3d

c              write(*,*)'delpwr_nc(is-1,ir)-delpwr_nc(is,ir)',
c     &                   delpwr_nc(is-1,ir)-delpwr_nc(is,ir)

c---------------------------------------------------------
              ratio=(delpwr_nc(is-1,ir)-delpwr_nc(is,ir))*
     &               del_s_pol_bin_lsc(j)/del_s_pol
     &               /delta_pow_loc
                      
c             write(*,*)'is,j,ratio',is,j,ratio
           else
             ratio=1.d0 
           endif
         else
           ratio=1.d0            
         endif
c         ratio=1.d0
c         write(*,*)'j,ratio',j,ratio
!YuP        endif !jjkk

c--------if D_QL is too big: D_QL*ratio > D_Coulomb*1.d14 then
c           it will be set D_QL=D_Coulomb*1.d14
c        else
c          incremental D_QL*ratio will be added to the total D_QL
c -------endif
         do nv = -nv_lsc_ql_lh,nv_lsc_ql_lh
            zeff_loc=zeffrho(rho_bin_center_lsc(j))
            d_cl_loc=d_cl(v_par_mesh_lsc(nv),zeff_loc)
            if( d_ql_lsc_lh_ar_loc(nv)*ratio.gt.d_cl_loc*1.d14) then
               d_ql_lsc_lh_ar_loc(nv)=d_cl_loc*1.d14
            else
              d_ql_lsc_lh_ar_loc(nv)=d_ql_lsc_lh_ar_loc(nv)*ratio
            endif
            d_ql_lsc_lh_ar(j,nv)=d_ql_lsc_lh_ar(j,nv)+
     &      d_ql_lsc_lh_ar_loc(nv)
         enddo
 
c------- test-----------------------------------------------
c        calculate power density
c        power_dens_watt_m3_loc, power_dens_erg_cm3_sec_loc
c        and power delta_pow_loc
c        from new total D_QL and old (fe and its derivative)
c        in j TCS radial bin
c----------------------------------------------------------- 
         power_dens_watt_m3_loc=0.d0
         do nv = -nv_lsc_ql_lh,nv_lsc_ql_lh
            power_dens_watt_m3_loc=
     &      power_dens_watt_m3_loc-
     &      v_par_mesh_lsc(nv)*d_ql_lsc_lh_ar(j,nv)*
     &      d_fe_dv_lsc(j,nv)*hv
         enddo

         dense_m3 = densrho(rho_bin_center_lsc(j),1)*1.d19  !electron density 1/m**3
         v_te_m_sec=v_te_bin_center_lsc(j)*1.d-2  ![m/sec]
         coef_watt_m3 = dense_m3*masse_kg*v_te_m_sec**2/
     &                  tau_n_bin_center_lsc(j)       !watt/m**3
         power_dens_watt_m3_loc=
     &          power_dens_watt_m3_loc*coef_watt_m3   !watt/m**3

         power_dens_erg_cm3_sec_loc=
     &          power_dens_watt_m3_loc*10.d0   !erg/cm**3/sec


         delta_pow_loc=power_dens_erg_cm3_sec_loc* binvol_lsc(j)
     
c         write(*,*)'from new total D_QL and old d_fe_dv:'
c         write(*,*)'j,delta_pow_loc',j,delta_pow_loc
c end test

      enddo    !j
       
c      write(*,*)'incremental QL power (not using ratio)'
c      write(*,*)'from old total fe from is ray bin'
c      write(*,*) 'in all j volumes: delta_pow_is [erg]', delta_pow_is
         
c-----test 

c-----calculation incremental QL power (after using ratio)
c     from old total fe from is ray bin
c     in all j =ibinmin,jbinmax volumes: delta_pow_is [erg]

      if (is.eq.1) delta_pow_is_m_1=0.d0

      delta_pow_is=0.d0

      do j=jbinmin,jbinmax 
         power_dens_watt_m3_loc=0.d0
         do nv = -nv_lsc_ql_lh,nv_lsc_ql_lh
            power_dens_watt_m3_loc=
     &      power_dens_watt_m3_loc-
     &      v_par_mesh_lsc(nv)*d_ql_lsc_lh_ar(j,nv)*
     &      d_fe_dv_lsc(j,nv)*hv
         enddo !nv
         dense_m3 = densrho(rho_bin_center_lsc(j),1)*1.d19  !electron density
                                                            ! 1/m**3
         v_te_m_sec=v_te_bin_center_lsc(j)*1.d-2  ![m/sec]
         coef_watt_m3 = dense_m3*masse_kg*v_te_m_sec**2/
     &                  tau_n_bin_center_lsc(j) !watt/m**3
         power_dens_watt_m3_loc=
     &          power_dens_watt_m3_loc*coef_watt_m3   !watt/m**3 
         power_dens_erg_cm3_sec_loc=
     &          power_dens_watt_m3_loc*10.d0   !erg/cm**3/sec  
         delta_pow_loc=power_dens_erg_cm3_sec_loc* binvol_lsc(j)
         delta_pow_is=delta_pow_is+delta_pow_loc
      enddo !j
c      write(*,*)'QL power from total D_QL using ratio'
c      write(*,*)'from old total fe from jbinmin=< j =<jbinmax'
c      write(*,*)'for given is ray bin'
c      write(*,*)'in all j volumes: delta_pow_is [erg]', delta_pow_is
         
c-----end test      

 130   continue

      return
      END

      subroutine LH_absorption_LSC(z,r,phi,cnpar,
     &cnper,cnz,cnr,cm,eps,d_f_e_d_v,dnpar,cnprim)
c-----Calculate------------------------------------
c     cold plasma dielectric tensor for LH wave case
c     with the thermal correctron for hermitian part
c     and imaginary eps_33_i for given  
c     electron distribution function f_e and its derivative
c     d_f_e_d_v,
c
c     D.W.Ignat, E.J.Valeo, S.C. Jardin,
c     Dynanmic modelling of lower hybrid current drive,
c     Nuclear Fusion, Vol 34, No.6(1994) p.837-852
c----------------------------------------------------------
      implicit none

      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'onetwo_nml.i'
      include 'lsc_approach.i'
      include 'rho.i'

c-----input
      real*8
     &z,r,phi,            !space coordinates
     &cnper,              !perpendicular refructive index  
     &cnpar,              !parallel      refructive index
     &cnz,cnr,  
     &cm,                  !n_phi=cm/r
     &dnpar                !variation of N_parallel at the ray point
                           !will be used for the calculation of
                           !the width of the delta function: sigma_lsc
c-----output
      complex*16
     &eps(3,3)            !dielectric tensor

      real*8
     &cnprim              !Im(N_perpendicular)

c-----local
      real*8    
     &v_gr_perp,lsc_eps_coef,
     &u(6),deru(6),vgr(3),vdotb,vgperps,bf(3),hrho,h_v,
     &v_te,      !theramal velocity sqrt(T_e/m_e)  [cm/c] 
     &mass_e,    ! electron mass [g]
     &charge_e,  ! electron charge (statcoulomb)
     &clight,    ! [cm/sec]
     &k_1_kev,   ! egrs in 1 KeV      (erg)    
     &sigma_lsc, !the width of the delta function 
     &d_eps_r_d_p_dev_d_eps_r_d_omega,K_zz_i,x_e,v_res,
c-----for test in SI
     & masse_kg,      ! electron mass [kg]
     & charge_e_c,    ! electron charge (coulomb)
     & k_si_1_kev,    ! Joules in 1 KeV      (Joule)
     & epsilon_0,
     &arg1,cln,dense,tau_n,temp_kev

      integer i, i_geom_optic_loc,j,k

      logical first

c------externals  
      external d_f_e_d_v  !calculate derivative df_dv_par 
      real*8 d_f_e_d_v,   !from the electron distribution
     &rhov,rhos,
     &b,temperho,densrho,x

      data first /.true./
      save first

      mass_e=9.1094d-28    ! electron mass [g]
      charge_e=4.8032d-10  ! electron charge (statcoulomb)
      clight=2.99792458d10 ! [cm/sec]
      k_1_kev=1.6022d-9    ! egrs i
c-----for test in SI
      masse_kg=9.1094d-31      ! electron mass [kg]
      charge_e_c=1.6022d-19    ! electron charge (coulomb)
      k_si_1_kev=1.6022d-16    ! Joules in 1 KeV      (Joule)
      epsilon_0=1.d0/(4.d0*pi*9.d9)

c      write(*,*)'in LH_absorption_LSC_new before lsc_eps first',first

c--------------------------------

      if (first) then
         first=.false.
c--------initialization of arrays      
         call ainalloc_lsc_approach_no_nml_i

         hrho=1.d0/n_psi_TSC   !step of TSC radial mesh

c------- create velocity mesh
         h_v=v_par_lsc_dev_vte_max/(nv_lsc_ql_lh)
         do j=-nv_lsc_ql_lh ,nv_lsc_ql_lh            
            v_par_mesh_lsc(j)=h_v*j !normalized at v_te
                                    !v_te is different 
                                    !at different radial points  
c            write(*,*)'set v_par_mesh_lsc(j),j',v_par_mesh_lsc(j),j
         enddo

c------- create arrays:
c        v_te_bin_center_lsc(i), 
c        tau_n_bin_center_lsc(i), i=1,n_psi_TSC         
         do i=1,n_psi_TSC+1
            rho_bin_lsc(i)=hrho*(i-1)
         enddo

         do i=1,n_psi_TSC
            rho_bin_center_lsc(i)=
     &           0.5d0*(rho_bin_lsc(i)+rho_bin_lsc(i+1))  
           
            binvol_lsc(i)=voltot*(rhov(rho_bin_lsc(i+1))**2-
     &                            rhov(rho_bin_lsc(i))**2)*1.d6 ![cm**3]

            binarea_lsc(i)=areatot*(rhos(rho_bin_lsc(i+1))**2-
     &                              rhos(rho_bin_lsc(i))**2)*1.d4 ![cm**2]

	    temp_kev=temperho(rho_bin_center_lsc(i),1)
            dense=densrho(rho_bin_center_lsc(i),1)   ! 10**13 /cm**3
            v_te=dsqrt(temp_kev*1.6022d-9/mass_e)    ! (sqrt(T_e/m_e) electron
                                                     ! thermal velocity [cm/sec]

            v_te_bin_center_lsc(i)=v_te               !cm/sec

c-----------cln -Coulomb Ln   
            arg1 = 1.d3/temp_kev*dsqrt(10.d0*dense) ! temp [KeV],
                                                    ! dense [10**13/cm**3]
            cln=24.d0-dlog(arg1)
c----------------------------------------------------------------------
c           tau_n character time used in Karney code
c------------------------------------------------------------------
c           tau_n=4*pi*epsolin_0**2*T**3/2*m**1/2/(e**4*dens*Cln) [SI]=
c              =T**3/2*m**1/2/(4*pi*e**4*dens*Cln)  [cgs]
c-------------------------------------------------------------------
c            tau_n=4.d0*pi*
c     &      epsilon_0**2*
c     &      (k_si_1_kev*temp_kev)**1.5d0*dsqrt(masse_kg)/
c     &      ((charge_e_c)**4*dense*1.d19*cln) ! sec
c            write(*,*)'0 tau_n',tau_n


            tau_n=1.d0/(4.d0*pi)*
     &      (k_1_kev*temp_kev)**1.5d0*dsqrt(mass_e)/
     &      (charge_e**4*dense*1.d13*cln) ! sec
c            write(*,*)'1 tau_n',tau_n


c            tau_n=1.d0/(4.d0*pi)*
c     &      (1.6022d0*temp_kev)*dsqrt(1.6022d0*temp_kev)*
c     &      dsqrt(9.1094d0)*
c     &      dsqrt(1.d-1)/
c     &      (4.8032d0**4* dense*cln) ! [sec]

            tau_n_bin_center_lsc(i)= tau_n

c            write(*,*)'tau_n',tau_n
         enddo
  
c--------set zero QL coefficients to test LSC approach
c        calculate array of derivatives from maxwellian fe

         do i=1,n_psi_TSC
            do j=-nv_lsc_ql_lh,nv_lsc_ql_lh 
               d_ql_lsc_lh_ar(i,j)=0.0d0 
            enddo
         enddo

         call calc_fe_lsc_arrays !at zero QL coefficient

      endif !first

c-----calculate small radius rho inside b,
c     rho is in common one.i
      
      bmod=b(z,r,phi)

      if (rho.gt.1.d0) then
          write(*,*)'WARNING:'
          write(*,*)'in LH_absorption_LSC_new rho.gt.1 rho=',rho
          write(*,*)'for rho>0 LSC approach can not work'
          write(*,*)'in this case it will be set cnprim_e=0'  
          cnprim=0.d0
          return
      endif

c-----calculate------------------------------------------
c     lsc_eps_coef=[(d_eps_r/d_p)*K_zz,I]/[d_eps_r/d_omega]
c---------------------------------------------------------   
      call  calc_d_eps_r_d_p_dev_d_eps_r_d_omega(z,r,phi,
     &      cnpar,cnper,cnz,cnr,cm,
     &      d_eps_r_d_p_dev_d_eps_r_d_omega) 

      hrho=1.d0/n_psi_TSC     ! radial step
      k=idint(rho/hrho)+1      ! the number of radial bin
      v_te=v_te_bin_center_lsc(k)   !(sqrt(T_e/m_e) electron thermal 
                                     !velocity [cm/sec]
      x_e=x(z,r,phi,1)                                  !for calculation K_xx_i
c      temp_kev=temperho(rho,1)
c      v_te=dsqrt(temp_kev*1.6022d-9/mass_e)             ! (sqrt(T_e/m_e) electron
                                                        ! thermal velocity [cm/sec]
      v_res=(clight/v_te)/cnpar                         ! resonance velocity
                                                        ! normalized to v_te(j)
                                                        ! v_te is different
                                                        ! at different
                                                        ! radial points j       
      sigma_lsc=dabs(v_res*dnpar/cnpar)

      if(dabs(v_res).lt.v_par_lsc_dev_vte_max)then
         call calc_K_zz_i(x_e,v_res,K_zz_i,sigma_lsc,rho) 
      else
         K_zz_i=0.d0
      endif

      lsc_eps_coef=d_eps_r_d_p_dev_d_eps_r_d_omega*K_zz_i
cSAP100519
c      write(*,*)'cnpar,v_res,d_eps_r_d_p_dev_d_eps_r_d_omega,K_zz_i',
c     &cnpar,v_res,d_eps_r_d_p_dev_d_eps_r_d_omega,K_zz_i
c      write(*,*)'lsc_eps_coef',lsc_eps_coef

c-----cnprim=eps_i/(d_eps_r_d_omega*v_gr_perp)

      u(1)=z
      u(2)=r
      u(3)=phi
      u(4)=cnz
      u(5)=cnr
      u(6)=cm
 
      i_geom_optic_loc=i_geom_optic
      i_geom_optic=1 !to get group velocity in deru 
      call rside1(0.d0,u,deru)
      i_geom_optic=i_geom_optic_loc
      
      vgr(1)=deru(1)
      vgr(2)=deru(2)
      vgr(3)=r*deru(3)

c-----------------------------------------------------------
c     vdotb -projection of the group velocity  on the
c            magnetic field multiplited by bmod
c-----------------------------------------------------------
      vdotb=vgr(1)*bz+vgr(2)*br+vgr(3)*bphi
c-----------------------------------------------------------
c     vgperps -perpendicular (to magnetic field)
c              component of the group velocity
c              in the second degree
c-----------------------------------------------------------
      bf(1)=bz
      bf(2)=br
      bf(3)=bphi

      vgperps=0.0d0
      do i=1,3
        vgperps=vgperps+(vgr(i)-vdotb*bf(i)/bmod**2)**2
      enddo  

c-----cnprim=-eps_i/(d_eps_r_d_omega*dsqrt(vgperps))
 
      cnprim=lsc_eps_coef/dsqrt(vgperps)

c      write(*,*)'-lsc_eps_coef,dsqrt(vgperps),cnprim',
c     &           -lsc_eps_coef,dsqrt(vgperps),cnprim

      if(cnprim.lt.0.d0) then
         write(*,*)'WARNING: in LH_absorption_LSC cnprim <0'
         write(*,*)'It was set cnprim=dabs(cnprim)'
         cnprim=dabs(cnprim)
      endif

      return
      end



c        ********************** prep3d_lsc***********************
c        *                      -----                       *
c        *  prep3d_lsc -subroutine to prepare the parameters    *
c        *  for	 output files for FP code  (e.g., CQL3D)    *
c        ****************************************************
c         input parameters: iray -number of the ray from antenna  
c                                 is in common/cone/ 
c         output parameter: iraystop (if power in the ray channel 
c         delpwr(is).lt. delpwrmn*delpwr(1), then iraystop=1). 
c 
c-----------------------------------------------------------------------
      subroutine prep3d_lsc(iray,is,iraystop,i_prof)

      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'write.i'
      include 'cefield.i'
c      include 'cone.i'
      include 'grill.i'
      include 'eps.i'
      include 'fourb.i'
      include 'oxb.i'
      include 'output.i'
      include 'three.i'
cSm070128 for test
      include 'six.i'

      include 'writencdf.i'
      include 'lsc_approach.i'
c-----input 
      integer
     &iray,  !number of the ray
     &is,    !number od ray element
     &i_prof !=0 do not calculate CD efficiency and
             !    power and current radial profiles 
             !    by call p_c_prof
             !=1  calculate CD efficiency and
             !    and power current radial profiles 
             !    by call p_c_prof
             

      real*8 vgr(3),bf(3)
      integer iraystop    
      
      real*8 tempiar(nbulka)

      complex*16 cnx
      
      real*8 optical_depth

c-----for reflection lost 
      integer irefl_old
      real*8  tot_pow_absorb_at_refl

c-------------------------------------------
      real*8 
     &u(6),deru(6),  
     &cnprim,argexp,ckvipl_e,ckvi, ckvipl_cl,ckvipl_i,
     &ckvipold,ckvipol,cm,cnper,cnpar,cnprim_e,cnprim_cl,
     &cnprim_i,cnz,cnr,cnt,cnprim_old,
     &dc,ds,dens_e,epar,enp,eplus,eminus,
     &ex,ey,ez,frqncy_l,z,r,phi,p_flux,psi_s,t,pwexp,one,
     &vgperps,r00,t00,temp_e,vdotb,v_gr_perp,vratio,vgrpls,vgrmods,
     &vgrpol,vgrs,wf,xe,ye,z_eff,
     &rhobegin,rhoend,p,delrho,
     &delta_power,lsc_eps_coef,vres,d_eps_r_d_p_dev_d_eps_r_d_omega,   
     &del_s_pol,hrho,rhomax,rhomin,K_zz_i,sigma_lsc,v_res,
     &z_m,r_m,delta_power_collisional
      integer i,id_old,i_geom_optic_loc,iflref_old,j,kk,jbinmax,jbinmin,
     &k

c-----external for lsc LH absorption
cBH110715      real*8  d_f_e_maxw_d_v_lsc,d_f_e_nonmaxw_d_v_lsc
cBH110715      external d_f_e_maxw_d_v_lsc,d_f_e_nonmaxw_d_v_lsc
      real*8   d_f_e_nonmaxw_d_v_lsc
      external d_f_e_nonmaxw_d_v_lsc
      real*8 b,gamma1,tempe,dense,zeff,x,y,cn
c----------------------------------
      data irefl_old /0/
      data tot_pow_absorb_at_refl/0.d0/
c-------------------------------------------
      data optical_depth/0.d0/
      
      save cnprim_old     
      save  optical_depth
      save irefl_old
      save tot_pow_absorb_at_refl

      save ckvipold

     
      iraystop=0
      pi=4.d0*datan(1.d0)
c----------------------------------------
c     cvac (cm/sec)
      cvac=2.997930D+10
c----------------------------------------
c     cld (cm),frgncy(GHz)
      cld=cvac/(2.d0*pi*frqncy*1.0d+09)
c----------------------------------------
c     now it is proposed that r0x=1 m
      r00=100.d0*r0x
    
c-----------------------------------------
c     z,r (m)
c-----cordinates of the end point of the ray element
      r_m=wr_nc(is,iray)*0.01d0/r0x ! normalization
      z_m=wz_nc(is,iray)*0.01d0/r0x ! normalization
      bmod=b(z_m,r_m,0.d0)
      rhoend=rho

c-----coordinates of the initial point of the ray element
      if (is.eq.1)then
         u(1)=wz_nc(is,iray)/r00
         u(2)=wr_nc(is,iray)/r00
         u(3)=wphi_nc(is,iray)
         u(4)=wn_z_nc(is,iray)
         u(5)=wn_r_nc(is,iray)
         u(6)=wn_phi_nc(is,iray)*u(2)
      else                           !data from (is-1) ray point
         u(1)=wz_nc(is-1,iray)/r00
         u(2)=wr_nc(is-1,iray)/r00
         u(3)=wphi_nc(is-1,iray)
         u(4)=wn_z_nc(is-1,iray)
         u(5)=wn_r_nc(is-1,iray)
         u(6)=wn_phi_nc(is-1,iray)*u(2)
      endif
      ws(is)=ws_nc(is,iray)
      wz(is) =wz_nc(is,iray)
      wr(is)=wr_nc(is,iray)
      wphi(is)=wphi_nc(is,iray)
      wn_z(is)=wn_z_nc(is,iray)
      wn_r(is)=wn_r_nc(is,iray)
      wn_phi(is)=wn_phi_nc(is,iray)

      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6)
c----------------------------------------------------
c     poloidal angle 0< theta_pol [degree] <360

c     bmod,bz,br,bphi (Tl)
      bmod=b(z,r,phi)
      rhobegin=rho
      bf(1)=bz
      bf(2)=br
      bf(3)=bphi
      gam=gamma1(z,r,phi,cnz,cnr,cm)
      ds=dsin(gam)
      dc=dcos(gam)
      cnt=cn(r,cnz,cnr,cm)
      cnpar=cnt*dc
      cnper=cnt*ds

c-------------------------------------------------------------
c      write(*,*)'prep3d_lsc bef LH_absorption_LSC'

c      write(*,*)'is,iray,z,r,phi,cnpar,cnper,cnz,cnr,cm',
c     &           is,iray,z,r,phi,cnpar,cnper,cnz,cnr,cm
c      if (is.gt.1) write(*,*)'wdnpar_nc(is-1,iray)',wdnpar_nc(is-1,iray)
 
      cnprim_cl=0.d0
      cnprim_e=0.d0

c      call absorplh(u,cnpar,cnper,temp_e,dens_e,tempiar
c     1                  ,bz,br,bphi,nbulk,bmod,frqncy,z_eff,
c     1                   cnprim_e,cnprim_i,cnprim_cl)
c      write(*,*)'LH iabsorp_2 cnprim_e', cnprim_e

c      call LH_absorption_LSC(z,r,phi,cnpar,
c     &           cnper,cnz,cnr,cm,reps, d_f_e_nonmaxw_d_v_lsc,
c     &           wdnpar_nc(is,iray),cnprim_e)

c      write(*,*)'prep3d_lsc LH LSC nonmaxw cnprim_e', cnprim_e  
 
c      call LH_absorption_LSC_maxwellian(z,r,phi,cnpar,
c     &                  cnper,cnz,cnr,cm,reps,cnprim_e)

c      write(*,*)'prep3d_lsc LH_LSC  cnprim_e', cnprim_e

      cnprim_i=0.d0
      cnprim_cl=0.d0
    
      if (ion_absorption.eq.'enabled') then
         cnprim=cnprim_e+cnprim_i
      else
         cnprim_i=0.d0 !YuP[2022-04-22] Added, for ion_absorption='disabled'
         !if(nbulk.gt.1) cnprim_s(2:nbulk)=0.d0 !YuP[2022-04-22] Added
         cnprim=cnprim_e
      endif
         
      cnprim=cnprim+cnprim_cl

c     cnprim=(cnprim_e+cnprim_i+cnprim_cl)
c      write(*,*)'cnprim_e,cnprim_i,cnprim_cl,cnprim',
c     1             cnprim_e,cnprim_i,cnprim_cl,cnprim
    
c---------------------------------------------------------
c     ws (cm)
      if (is.eq.1) then               
         i_ox_conversion=0
      end if

cSAP100519
c      write(*,*)'is,ws',is,ws_nc(is,iray)
    
c-----calculate the group velocity
      id_old=id

      i_geom_optic_loc=i_geom_optic
      i_geom_optic=1 !to get group velocity in deru
      t=0.d0 
      call rside1(t,u,deru)  !at is-1 ray point
      i_geom_optic=i_geom_optic_loc

      vgr(1)=deru(1)
      vgr(2)=deru(2)
      vgr(3)=r*deru(3)
      
      vgrmods=vgr(1)**2+vgr(2)**2+vgr(3)**2
      if(vgrmods.gt.1.01d0) then !YuP[2020-08-17] Changed gt.1 to gt.1.01
         !otherwise too much of printout; 
         ! vgroup=1.00005 (slightly larger than 1) might happen in near-vacuum
         write(*,*)
         write(*,*) '*************************************************'
       write(*,*)'prep3d_lsc WARNING: vgroup>1, vgroup=',dsqrt(vgrmods)
         write(*,*) '*************************************************'
         write(*,*)
      endif

c-----------------------------------------------------------
c     vdotb -projection of the group velocity  on the
c            magnetic field multiplited by bmod
c-----------------------------------------------------------
      vdotb=vgr(1)*bz+vgr(2)*br+vgr(3)*bphi
c-----------------------------------------------------------
c     vgperps -perpendicular (to magnetic field)
c              component of the group velocity
c              in the second degree
c-----------------------------------------------------------
      vgperps=0.0d0
      do i=1,3
        vgperps=vgperps+(vgr(i)-vdotb*bf(i)/bmod**2)**2
      enddo

c----------------------------------------------------------
c     collisional damping  !!! to change for is point 
c----------------------------------------------------------
      if( iabsorp_collisional.eq.1) then
         temp_e=tempe(z,r,phi,1)
         dens_e=dense(z,r,phi,1)
         z_eff=zeff(z,r,phi)
         frqncy_l=frqncy
         v_gr_perp= dsqrt(vgperps)
         call absorp_collisional(temp_e,dens_e,frqncy_l,z_eff,
     &   v_gr_perp,coll_mult,
     &   cnprim_cl)
c         write(*,*)'prep3d cnprim,cnprim_cl',cnprim,cnprim_cl
         cnprim=cnprim+cnprim_cl

c         cnprim_e=cnprim_e+cnprim_cl
      endif !iabsorp_collisional=1 
c----------------------------------------------------------
      vgrs=vgr(1)**2+vgr(2)**2+vgr(3)**2
      ckvi=dsqrt(vgperps/vgrmods)*cnprim/cld
      vgrpls=vgr(1)**2+vgr(2)**2
      vgrpol=dsqrt(vgrpls)
      wf=frqncy
c----------------------------------------------------------
c     ckvipol  (1/cm)
      vratio=dsqrt(vgperps/vgrpls)
      ckvipol=vratio*cnprim/cld

      ckvipl_e=vratio*cnprim_e/cld
      ckvipl_i=vratio*cnprim_i/cld
      ckvipl_cl=vratio*cnprim_cl/cld
 
c      write(*,*)'prep3d: ckvipl_e,ckvipl_i,ckvipl_cl',
c     &                   ckvipl_e,ckvipl_i,ckvipl_cl

      seikon(is)=0.d0
c---------------------------------------------------------
c     wr and wz (cm),wphi (radian)
c      wr(is)=r*r00
c      wphi(is)=phi
c      wz(is)=z*r00
c      wnpar(is)=cnpar
c      write(*,*)'prep3d wnpar(is)',wnpar(is)
c      wnper(is)=cnper
c      write(*,*)'prep3d is,cnper,wnper(is)',is,cnper,wnper(is)
c      wmtor(is)=u(6)
c-------------------------------------------------------------------
c     Here delpwr(is) is the power(erg/c) in the ray
c     channel.It is equal powini_iray at antenna.
      if (is.eq.1) then
c        powj(iray) (erg/c) was calculated in cone_ec
c        powinilh(iray) (erg/c) was calculated in grill_lh
c        powini=powj or powinilh        
         p=powini
         write(*,*)'prep3d is=1,powini',powini

         delpwr(is)=p
c         write(*,*)'is=1 delpwr(1)',delpwr(1)
         optical_depth=0.d0
    
         iflref_old=0
         ckvipold=ckvipol
      else
       
         argexp=(-(ckvipol+ckvipold)*(ws_nc(is,iray)-ws_nc(is-1,iray)))
         ckvipold=ckvipol
c         write(*,*)'in prep3d_lsc ckvipol,argexp',ckvipol,argexp
	 if (dabs(argexp).gt.90.d0) then
	    write(*,*)'in prep3d argexp.gt.90',argexp
            delpwr(is)=0.d0
         else 
            pwexp=dexp(argexp)
	    if (pwexp.lt.1.d-50) then
	         write(*,*)'in prep3d pwexp.lt.1.d-50',pwexp
                 delpwr(is)=0.d0
	     endif
	 endif
         optical_depth=optical_depth+dabs(argexp)
c        write(*,*)'optical_depth',optical_depth
                  
c-YuP-130605: changed ray-stopping criterion from argexp>0.d0 to this:
           if(argexp.gt. 1.d-30)then
              WRITE(*,*)'**********************************************'
              WRITE(*,*)'WARNING: in lsc_approach: argexp>0 :' , argexp
              WRITE(*,*)'It would give growing ray power-> stopping ray'
              WRITE(*,*)'ckvipol,ckvipold=',ckvipol,ckvipold
              WRITE(*,*)'**********************************************'
              argexp=0.d0
              iraystop=1
              iray_status_one_ray=10
           endif
c-YuP-130605: From print-out: 
c Even though sometimes ckvipol and ckvipold both are zero,
c yet the value of  argexp= -(ckvipol+ckvipold)*(ws(is)-ws(is-1))
c is not zero (but rather a small number ~ 1e-321).
c Because of this seemingly insignificant rounding error,
c the rays were stopped prematurely.
c It only happens on Hopper/PGI, not on IntelVisualFortran.

c        write(*,*)'delpwr(is-1),argexp',delpwr(is-1),argexp
c         delpwr(is)=delpwr(is-1)*dexp(argexp)
         delta_power_collisional=delpwr(is-1)*(1.d0-dexp(argexp))
c         write(*,*)'perp3d_lsc  power from exp(argexp) '
c         write(*,*)'is.delpwr(is-1),delpwr(is),delpwr(is-1)-delpwr(is)',
c     &              is,delpwr(is-1),delpwr(is),delpwr(is-1)-delpwr(is)
cSAP100222
c--------calculate------------------------------------------
c        lsc_eps_coef=[(d_eps_r/d_p)*K_zz,I]/[d_eps_r/d_omega]
c---------------------------------------------------------
         call  calc_d_eps_r_d_p_dev_d_eps_r_d_omega(z,r,phi,
     &   cnpar,cnper,cnz,cnr,cm,
     &   d_eps_r_d_p_dev_d_eps_r_d_omega)

c------ -calculate jbinmin,jbinmax for given rhobegin,rhoend
c        j is the number of the central bin radial point of TSC mesh

         if(rhoend.gt.rhobegin) then
           rhomax=rhoend
           rhomin=rhobegin
         else
	   rhomin=rhoend
	   rhomax=rhobegin
         endif

         if(rhomax.gt.1.d0) rhomax=1.d0
         if(rhomin.gt.1.d0) rhomin=1.d0-1.d-10    

         delrho=rhomax-rhomin

         hrho=1.0d0/n_psi_TSC 
cSAP100504--------like in the old version
         k=idint(rho/hrho)+1      ! the number of radial bin  
         if(k.gt.n_psi_TSC) k=n_psi_TSC
         v_res=(cvac/v_te_bin_center_lsc(k))/cnpar   
         sigma_lsc=dabs(v_res*wdnpar_nc(is,iray)/cnpar)
         xe=x(z,r,phi,1)
         if(dabs(v_res).lt.v_par_lsc_dev_vte_max)then  
              call calc_K_zz_i(xe,v_res,K_zz_i,      !rho
     &        sigma_lsc,rho_bin_center_lsc(k))       !for j point of TSC mesh
                                                     !xe from 'is-1' poin  
         else
              K_zz_i=0.d0
         endif       
         lsc_eps_coef=d_eps_r_d_p_dev_d_eps_r_d_omega*K_zz_i  

cSAP100519
c         write(*,*)'cnpar,v_res,d_eps_r_d_p_dev_d_eps_r_d_omega,K_zz_i',
c     &              cnpar,v_res,d_eps_r_d_p_dev_d_eps_r_d_omega,K_zz_i
         delta_power = 
     &               2.d0*lsc_eps_coef*delpwr(is-1)*
     &               2.d0*pi*frqncy*1.d9*
     &               (ws_nc(is,iray)-ws_nc(is-1,iray))/(vgrpol*cvac)  
cSAP100519
c          write(*,*)'is,lsc_eps_coef,delpwr(is-1),delta_power',
c     &               is,lsc_eps_coef,delpwr(is-1),delta_power

         goto 10
c------------------------------------
         jbinmin=rhomin/hrho+1
         jbinmax=rhomax/hrho+1

         if(jbinmin.gt.n_psi_TSC)jbinmin=n_psi_TSC 
         if(jbinmax.gt.n_psi_TSC)jbinmax=n_psi_TSC 
     
         del_s_pol=ws_nc(is,iray)-ws_nc(is-1,iray)

         del_s_pol_bin_lsc(jbinmin)=del_s_pol*
     &              (hrho*jbinmin-rhomin)/delrho

         del_s_pol_bin_lsc(jbinmax)=del_s_pol*
     &               (rhomax-hrho*(jbinmax-1))/delrho

         if(jbinmax.gt.(jbinmin+1)) then
           do j=(jbinmin+1),(jbinmax-1)
              del_s_pol_bin_lsc(j)=del_s_pol*hrho/delrho
           enddo    
         endif

c        call calc_jbinmax_jbinmin(rhobegin,rhoend,
c     &  n_psi_TSC,radial_mesh,jbinmax,jbinmin)

c------ -calculate delta_power absorped at (is-1,is) ray bin 
         xe=x(z,r,phi,1)  !for calculation K_xx_i

         delta_power=0.d0

         do j=jbinmin,jbinmax
           
           v_res=(cvac/v_te_bin_center_lsc(j))/cnpar      ! resonance velocity
                                                          ! normalized to v_te
                                                          ! v_te is different
                                                          ! at different
                                                          ! radial points j       
           sigma_lsc=dabs(v_res*wdnpar_nc(is,iray)/cnpar)
c           write(*,*)'prep3d bef K_zz j,xe,v_res,sigma_lsc',
c     &                                j,xe,v_res,sigma_lsc
c----------------------------------------------------------------- 
           if(dabs(v_res).lt.v_par_lsc_dev_vte_max)then  
              call calc_K_zz_i(xe,v_res,K_zz_i,      !rho
     &        sigma_lsc,rho_bin_center_lsc(j))       !for j point of TSC mesh
                                                     !xe from 'is-1' poin  
           else
              K_zz_i=0.d0
           endif                                     !xe from 'is-1' point
           lsc_eps_coef=d_eps_r_d_p_dev_d_eps_r_d_omega*K_zz_i 

           delta_power = delta_power+
     &               2.d0*lsc_eps_coef*delpwr(is-1)*
     &               2.d0*pi*frqncy*1.d9*
     &               del_s_pol_bin_lsc(j)/(vgrpol*cvac) 
     
         enddo !j 

 10      continue

c         write(*,*)'delpwr(is-1),delta_power',delpwr(is-1),delta_power
         delta_power=dmin1(delpwr(is-1)*(1.d0-1.d-15),delta_power)

c         write(*,*)'dmin1(delpwr(is-1)*(1.d0-1.d-15),delta_power)',
c     &              delta_power
c         delpwr(is)=delpwr(is-1)- delta_power  
         delpwr(is)=delpwr(is-1)- delta_power-delta_power_collisional
cSAP100519
c         write(*,*)'prep3d_lsc delpwr(is),delta_power', 
c     &                         delpwr(is),delta_power

c         write(*,*)'delpwr(is-1),delpwr(is),delpwr(is-1)-delpwr(is)',
c     &              delpwr(is-1),delpwr(is),delpwr(is-1)-delpwr(is)
c-------------------------------------------------------------------
c        reflection lost at the plasma edge
c------------------------------------------------------------------
c         write(*,*)'prep3d refl_loss,irefl,irefl_old',
c     &             refl_loss,irefl,irefl_old
         tot_pow_absorb_at_refl=tot_pow_absorb_at_refl+
     &           delpwr(is)*refl_loss*(irefl-irefl_old)
      
c         write(*,*)'prep3d before refl_looss delpwr(i)',delpwr(is)

         delpwr(is)=delpwr(is)*(1.d0-refl_loss*(irefl-irefl_old))
         irefl_old=irefl 
         w_tot_pow_absorb_at_refl=tot_pow_absorb_at_refl
c-------------------------------------------------------------------    
c         write(*,*)'prep3d delpwr(is)',delpwr(is)
c         write(*,*)'delpwr(1)',delpwr(1)

         if((delpwr(is).gt.1.d-200).and.(delpwr(1).gt.1.d-200))then
c           if(i_ox.ne.1) write(*,*)'-dlog(delpwr(is)/delpwr(1))',
c     &             -dlog(delpwr(is)/delpwr(1))
         endif
 
c         if(argexp.gt.0.d0)then
c           write(*,*)'******************************************'
c           write(*,*)'******************************************'
c           write(*,*)'WARNING: in prep3d.f argexp>0' 
c           write(*,*)'It will give the growing ray power'
c          write(*,*)'******************************************'
c         endif
           
c         if(delpwr(is).lt.delpwrmn*delpwr(1))then
c            write(*,*)'***in prep3d delpwr(is).lt.delpwrmn*delpwr(1)**'
c           stop ray_iray calculations
c            iraystop=1
c         endif
      end if !is=1
      
c----------------------------------------------------------------------

c     grill conditions for the waves (LH or FW)
c     wdnpar(is)=wdnpar0(iray)
c     write(*,*)'in prep3d new wdnpar(is)',wdnpar(is)
      
c---------------------------------------------------------
c     vgrpol(cm/sec)=vgrpol*r00/t00
c     vrgpol/c=vgrpol/wf
c---------------------------------------     
      salphac(is)=2.d0*ckvipl_cl    ! collisional damping coefficient
      sdpwr(is)=  2.d0*ckvipl_i     ! ion damping coefficient
      salphal(is)=2.d0*ckvipl_e     ! electron damping coefficient
c      write(*,*)'salphal(is)=2.d0*ckvipl_e',salphal(is)
      if (is.gt.1) then
c         write(*,*)'is,(0.5d0*(ws(is)-ws(is-1))),ws(is),ws(is-1)',
c     &              is,(0.5d0*(ws(is)-ws(is-1))),ws(is),ws(is-1)
c         write(*,*)'delpwr(is-1),salphal(is-1),delpwr(is)',
c     &              delpwr(is-1),salphal(is-1),delpwr(is)
cSAP100528
c         salphal(is)=(delta_power/(0.5d0*(ws(is)-ws(is-1)))-
c     &   delpwr(is-1)*salphal(is-1))/delpwr(is)
c         write(*,*)'0.5d0*(ws(is)-ws(is-1))*(delpwr(is-1)*salphal(is-1)
c     &+delpwr(is)*salphal(is))',
c     &0.5d0*(ws(is)-ws(is-1))*(delpwr(is-1)*salphal(is-1)
c     &+delpwr(is)*salphal(is))
         if((ws(is)-ws(is-1)).lt.1d-20) then
           salphal(is)=salphal(is-1)
         else
           if(delpwr(is).gt.1.d-20) then
             salphal(is)=delta_power/((ws(is)-ws(is-1))*delpwr(is))
c             write(*,*)' salphal(is)',salphal(is)
         
c             write(*,*)'(ws(is)-ws(is-1))*delpwr(is)*salphal(is))',
c     &                  (ws(is)-ws(is-1))*delpwr(is)*salphal(is)
c             write(*,*)'delta_power',delta_power
           endif
        endif
      endif 
c-----------------------------------------------------------
c     data for onetwo
c-----------------------------------------------------------i
      if ((i_prof.eq.1).and.(ionetwo.eq.1)) then        
c        call p_c_prof_lsc(is,rhoold,rho,cnpar,cnper,cex,cey,cez)
c-------calculate: rho,cnpar.cnper at is point
        z=   wz(is)/r00
        r=   wr(is)/r00
        phi= wphi(is)
        cnz= wn_z(is)
        cnr= wn_r(is)
        cm=  wn_phi(is)*r
        bmod=b(z,r,phi)
        gam=gamma1(z,r,phi,cnz,cnr,cm)
        ds=dsin(gam)
        dc=dcos(gam)
        cnt=cn(r,cnz,cnr,cm)
        cnpar=cnt*dc
        cnper=cnt*ds

c        write(*,*)'in prep3d_lsc before p__c_prof is',is  
c        write(*,*)'z,r,phi,cnz,cnr,cm,rho',z,r,phi,cnz,cnr,cm,rho
c        write(*,*)'ds,dc,gam,cnt,cnpar,cnper',ds,dc,gam,cnt,cnpar,cnper

cSAP100528
c        call p_c_prof(is,rhoold,rho,cnpar,cnper,cex,cey,cez)
        call p_c_prof_lsc(is,rhoold,rho,cnpar,cnper,cex,cey,cez)

      endif

c      write(*,*)'in prep3d_lsc after p__c_prof'      
      rhoold=rho

      return
      END

      subroutine calc_jbinmax_jbinmin(rhobegin,rhoend,
     &n_rho,radial_mesh,jbinmax,jbinmin)
c-----calculate
c     minimal and maximal numbers of uniform radial mesh point
c     jbinmax,jbinmin
c     (rhobegin,rhoend) is inside 
c     (rho_mesh(jbinmin-1),rho_mesh(jbinmin))

      implicit none
c-----input
      integer
     &n_rho                           !dimension of radial mesh

      real*8
     &rhobegin,rhoend,
     &radial_mesh(1:n_rho)            !uniform radial mesh
                                      !radial_mesh(1)=0
                                      !radial_mesh(n_rho)=1
c-----output
      integer jbinmax,jbinmin
  
c-----locals
      integer j
      real*8 rhomax,rhomin,hrho

      if(rhoend.gt.rhobegin) then
         rhomax=rhoend
         rhomin=rhobegin
      else
	 rhomin=rhoend
	 rhomax=rhobegin
      endif

      if(rhomax.gt.1.d0) rhomax=1.d0
      if(rhomin.gt.1.d0) rhomin=1.d0-1.d-10    

      hrho=radial_mesh(2)-radial_mesh(1)

      jbinmin=rhomin/hrho+1
      jbinmax=rhomax/hrho+1

      if(jbinmin.gt.n_rho)jbinmin=n_rho  
      if(jbinmax.gt.n_rho)jbinmax=n_rho
          
      return
      end

      subroutine LSC_power_iterations
c-----Do iterations by absorpted power at nonmaxwellian 
c     electron distribution
      implicit none

      include 'param.i'
      include 'grill.i'
      include 'one.i'                 
      include 'lsc_approach.i'
      include 'write.i'
      include 'writencdf.i'
      include 'onetwo.i'   

c-----locals
      integer i_lsc_pow,
     &iray,is,i,nv,iraystop,k,n,
     &icount_iter,     !number of iterations     
     &ir_plot,         !number of radial point to plot D_QL,fe,
                       !integral at the exponent argument of fe
                       !ir_plot is set inside this subroutine FOR SOME CASE.
     &i_profile,       !=0 do not calculate CD efficiency and
                       !    power and current radial profiles 
                       !    by call p_c_prof
                       !=1  calculate CD efficiency and
                       !    power and current radial profiles 
                       !    by call p_c_prof
     &i_power_step,
     &i_rho

      real*8
     &rhobegin,rhoend,cnpar,cnper,
     &diff_d_q,diff_fe,diff_delpwr,coef_power
           
      write(*,*)'in LSC_power_iterations icount_iter_max =',
     &icount_iter_max
     
      icount_iter=0          
      i_profile=0

c-----initialization of arrays       
      d_ql_lsc_lh_ar=0.d0     !set zero for array of D_QL
      call calc_fe_lsc_arrays !fe and its derivative at zero  d_ql_lsc_lh_ar

c---- loop in power to use growing initial power  
      write(*,*)'in  LSC_power_iterations n_power_steps',n_power_steps

      do i_power_step=1,n_power_steps     
         coef_power=i_power_step*1.d0/n_power_steps !sets the part 
                                                    !of the initial ray power
         write(*,*)'i_power_step,coef_power',i_power_step,coef_power
         icount_iter=0
         i_profile=0
 10      continue  !loop back point for power interation.

c         fe_lsc_old=fe_lsc
c         call calc_fe_lsc_arrays 
c--------------------------------------------------------------------
         if ((i_profile.eq.1).and.(ionetwo.eq.1)) call onetwoini
c--------------------------------------------------------------------
         delpwr_nc_old_lsc=delpwr_nc      
         d_ql_lsc_lh_ar_old=d_ql_lsc_lh_ar

         d_ql_lsc_lh_ar=0.d0     !set zero for array of D_QL
c-------------------------------------------------------------------
c        calculate QL coefficients from all rays
c-------------------------------------------------------------------       
         write(*,*)'before call d_ql_lh_prof'
         do iray=1,nrayl
c-----------the loop other all rays
            do is=2,nrayelt_nc(iray)
c--------------the loop other ray elements    
               call d_ql_lh_prof(is,iray)
            enddo
         enddo  
         write(*,*)'after call d_ql_lh_prof'
c----------------------------------------------------------------
c       calculate norm of arrays difference
c       diff_d_q=||d_ql_lsc_lh_ar_old- d_ql_lsc_lh_ar||
c----------------------------------------------------------------
         call norm_difference(1,n_psi_TSC,-nv_lsc_ql_lh,nv_lsc_ql_lh,
     &   d_ql_lsc_lh_ar_old,d_ql_lsc_lh_ar,diff_d_q)
c         write(*,*)'icount_iter,diff_d_q', 
c      &              icount_iter,diff_d_q          

c--------averaging procedure D_QL=0.5(D_QL_old+D_QL)
c        do i=1,n_psi_TSC
c          do nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
c            d_ql_lsc_lh_ar(i,nv)=0.5d0*
c     &      (d_ql_lsc_lh_ar(i,nv)+d_ql_lsc_lh_ar_old(i,nv))
c          enddo
c        enddo
c--------------------------------------------
         fe_lsc_old=fe_lsc
         call calc_fe_lsc_arrays !fe and df_dv at new D_QL 
                                 !fe is normaiazed at unit 
         call norm_difference(1,n_psi_TSC,-nv_lsc_ql_lh,nv_lsc_ql_lh,  
     &   fe_lsc,fe_lsc_old,diff_fe)
 
c----------------------------------------------------------------
c        plot at ir_plot radial point:
c        D_ql(v_parallel)
c----------------------------------------------------------------

         ir_plot=69
         ir_plot=10
         ir_plot=90 !for iter lh case
         ir_plot=0.7*n_psi_TSC
c        ir_plot=35
c        ir_plot=137

c         call plot_lsc_d_ql(ir_plot) !plot D_QL(V_parallel)
         
c         call plot_log_fe(ir_plot)   !plot log_fe(V_parallel)
     
c         call plot_integral_lsc(ir_plot) !plot integral(V_parallel)
                                         !in exponent of fe    
c------------------------------------------------------------------
c        calculate absorbed power along all rays
c-------------------------------------------------------------------
         delpwr_nc_old_lsc=delpwr_nc
         
         call total_power_fe  !for test only, print total fe power egr/sec

         write(*,*)'before loop for prep3d'

         do iray=1,nrayl   
c-----------the loop other all rays
            powini=powinilh(iray)
            powini=coef_power*powinilh(iray)
         
            write(*,*)'iray,icount_iter',iray,icount_iter

            write(*,*)'coef_power,powinilh(iray),powini ',
     &                 coef_power,powinilh(iray),powini 

c------------------------------------------------------------------
c           calculate power along all rays for new fe with new D_QL 
c----------------------------------------------------------------------
            if (nrayelt_nc(iray).gt.0) then
               do is=1,nrayelt_nc(iray)
c-----------------the loop other ray elements
                  call prep3d_lsc(iray,is,iraystop,i_profile)
               enddo
            endif !nrayelt_nc(iray).gt.0

            do is=1,nrayelt      
               delpwr_nc(is,iray)=delpwr(is)
            enddo
cBH140729:  Zero out remainder of delpwr_nc, since may contain
cBH140729:  prior interation data.
            do is=nrayelt+1,nrelta
               delpwr_nc(is,iray)=0.d0
            enddo

            call norm_power_difference(diff_delpwr)
c--------------------------------------------------------------------
c           calculation of power(array spower(NR)) and
c           current(array scurrent(NR)) radial profiles
c           as sum of profiles Sum(i=1,iray)
c--------------------------------------------------------------------
            if ((i_profile.eq.1).and.(ionetwo.eq.1)) then
ctest
c         do i=1,NR-1         
c             write(*,*)'before sonetwo i,power(i)',i,power(i)
c         enddo
cend test
              call sonetwo
            endif

c            write(*,*)'i_power_step,icount_iter',
c     &                 i_power_step,icount_iter
c            write(*,*)'diff_d_q',diff_d_q
c            write(*,*)'diff_fe',diff_fe
c            write(*,*)'diff_delpwr',diff_delpwr
c--------------------------------------------------------------------
         enddo ! iray=1,nrayl
         write(*,*)'after loop for prep3d'
c--------------------------------------------------------------------
c        plot old and new power along all rays 
c        after one 'icount_iter' iteration delpwr(s_poloidal)
c--------------------------------------------------------------------
         call plot_lsc_powers

         call norm_power_difference(diff_delpwr)

         icount_iter=icount_iter+1
         write(*,*)'in loop in power:i_power_step,icount_iter',
     &                               i_power_step,icount_iter
       
         write(*,*)'diff_d_q',diff_d_q
         write(*,*)'diff_fe',diff_fe
         write(*,*)'diff_delpwr',diff_delpwr

         write(*,*)'subroutine LSC_power_iterations before'
c         write(*,*)'before goto 10'      
         write(*,*)'i_profile',i_profile
         call total_power_fe  !for test only

         if(i_profile.eq.1) goto 30

         if(icount_iter.gt.icount_iter_max)then
            write(*,*)'icount_ier.gt.icount_iter_max' 
            i_profile=1   
            goto 10
         endif

         if((i_power_step.eq.n_power_steps).and.
     &      (icount_iter.gt.icount_iter_max))then
            write(*,*)'i_power_step.eq.n_power_steps and'
            write(*,*)'icount_ier.gt.icount_iter_max' 
            i_profile=1   
            goto 10
         endif

         if  (diff_d_q.gt.eps_diff_d_ql) then 
            goto 10
         else
            i_profile=1  !=1  to calculate CD efficiency and
                         !    power and current radial profiles 
                         !    by call p_c_prof
            goto 10
         endif  

 30      continue

         write(*,*)'subroutine LSC_power_iterations before'
         write(*,*)'enddo  i_power_step=1,n_power_steps'
         write(*,*)'i_power_step',i_power_step
         call total_power_fe  !for test only

c--------change delpwr for new power level
         if (i_power_step.lt.n_power_steps) then
            do iray=1,nrayl
               do is=1,nrayelt
                 delpwr_nc(is,iray)=delpwr_nc(is,iray)*
     &           (i_power_step-1.d0)/i_power_step
               enddo
            enddo
         endif

      enddo !i_power_step=1,n_power_steps
c---------------------------------------------------------------------------
c     calculate
c     absorbed LH power density raial profile
c     and RF current density radial profile without Edc
c     using non-Maxwellian electron distribution
c---------------------------------------------------------------------------
      write(*,*)'LSC_power_iterations before'//
     &' calc_RF_power_CD_density_profile'
      call total_power_fe  !for test only
    
c      do iray=1,nrayl
c        write(*,*)'iray,nrayelt_nc(iray)',iray,nrayelt_nc(iray)
c        do is=1,nrayelt_nc(iray)
c           write(*,*)'iray,is,delpwr_nc(is,iray),ws_nc(is,iray)',
c     &                iray,is,delpwr_nc(is,iray),ws_nc(is,iray)
c        enddo
c      enddo

      call MK_GRAPH_delpwr
cfor test
c      call plot_lsc_powers_1
cend for test
c      call calc_RF_power_CD_density_profile

      return
      end

      subroutine total_power_fe
c-----calculate total power from fe
      
      implicit none 
      include 'param.i'
      include 'lsc_approach.i'
      include 'onetwo.i'
      include 'one_no_nml.i'

c-----locals
      integer nv
 
      real*8
     &dense_m3,                ! electron density 1/m**3
     &masse_kg,                ! electron mass [kg]
     &coef_watt_m3,
     &hv,v_te_m_sec,
c-----for test     
     &delta_pow_is,power_dens_watt_m3_loc,power_dens_erg_cm3_sec_loc,
     &delta_pow_loc
      integer j
c-----externals
      real*8 densrho,zeffrho

      masse_kg=9.1094d-31      ! electron mass [kg]

      hv= v_par_mesh_lsc(1)-v_par_mesh_lsc(0)
      delta_pow_is=0.d0 
      do j=1,n_psi_TSC
         power_dens_watt_m3_loc=0.d0
         do nv = -nv_lsc_ql_lh,nv_lsc_ql_lh
            power_dens_watt_m3_loc=
     &      power_dens_watt_m3_loc-
     &      v_par_mesh_lsc(nv)*d_ql_lsc_lh_ar(j,nv)*
     &      d_fe_dv_lsc(j,nv)*hv
         enddo !nv
         dense_m3 = densrho(rho_bin_center_lsc(j),1)*1.d19  !electron density
                                                            ! 1/m**3
         v_te_m_sec=v_te_bin_center_lsc(j)*1.d-2  ![m/sec]
         coef_watt_m3 = dense_m3*masse_kg*v_te_m_sec**2/
     &                  tau_n_bin_center_lsc(j) !watt/m**3
         power_dens_watt_m3_loc=
     &          power_dens_watt_m3_loc*coef_watt_m3   !watt/m**3 
         power_dens_erg_cm3_sec_loc=
     &          power_dens_watt_m3_loc*10.d0   !erg/cm**3/sec  
         delta_pow_loc=power_dens_erg_cm3_sec_loc* binvol_lsc(j)
         delta_pow_is=delta_pow_is+delta_pow_loc
      enddo !j

      write(*,*)'in subroutine: total_power_fe'
      write(*,*)'total_power=delta_pow_is [egr/sec] ',delta_pow_is
 
      return
      end

      subroutine norm_difference(ndim1_low,ndim1_top,
     &ndim2_low,ndim2_top,a,bnd,norm_diff)
      !YuP[2020-08-20] renamed b to bnd, to avoid conflict with function b()
c-----calculates norm of arrays a and b difference
      implicit none
c-----input
      integer ndim1_low,ndim1_top,ndim2_low,ndim2_top   !arrays dimensions

      real*8
     &a(ndim1_low:ndim1_top,ndim2_low:ndim2_top),
     &bnd(ndim1_low:ndim1_top,ndim2_low:ndim2_top)  !2D arrays
   

c-----output
      real*8 norm_diff !norm of difference

c-----locals
      real*8 eps,dif
      integer i,n
    
      norm_diff=0.d0

      eps=1.d-32

      do i=ndim1_low,ndim1_top           	 
         do n=ndim2_low,ndim2_top
            dif=dabs(a(i,n)-bnd(i,n))/dabs(a(i,n)+bnd(i,n)+eps)
c            dif=dabs(a(i,n)-bnd(i,n))
            if(norm_diff.lt.dif) norm_diff=dif
c            norm_diff=norm_diff+
c     &      2.d0*dabs(a(i,n)-bnd(i,n))/dabs(a(i,n)+bnd(i,n)+eps)
         enddo
      enddo

c      norm_diff= norm_diff/
c     &           ((ndim1_top-ndim1_low)*(ndim2_top-ndim2_low))     

      return
      end subroutine norm_difference


      subroutine norm_power_difference(norm_diff)
c-----calculates norm of arrays delpwr_nc and delpwr_nc_old_lsc difference
      implicit none

      include 'param.i'
      include 'writencdf.i'
      include 'lsc_approach.i'


c-----output
      real*8 norm_diff !norm of difference

c-----locals
      real*8 eps,dif
      integer iray,is,n_norm
    
      norm_diff=0.d0
      n_norm=0
      eps=1.d-32

      do iray=1,nrayl           	 
         do is=1,nrayelt_nc(iray)    
            dif=dabs(delpwr_nc_old_lsc(is,iray)-delpwr_nc(is,iray))/
     &      dabs(delpwr_nc_old_lsc(is,iray)+delpwr_nc(is,iray)+eps)
            if(norm_diff.lt.dif) norm_diff=dif

c            norm_diff=norm_diff+
c     &      2.d0*dabs(delpwr_nc_old_lsc(is,iray)-delpwr_nc(is,iray))/
c     &      dabs(delpwr_nc_old_lsc(is,iray)+delpwr_nc(is,iray)+eps)
c            n_norm=n_norm + 1
         enddo
      enddo

c      norm_diff=norm_diff/n_norm

      return
      end

      real*8 function d_f_e_nonmaxw_d_v_lsc(v_par,rho)
c-----nondimensional
c     derivative d_f_dv_par
c     from 1D maxwellian distribution
c     versus v_par 
c     v_par is normalized to v_t parallel velocity
c
c     the normalization condition is integral(f_par*dv_par)=1
c
c     It uses pre-calculated  array 
c     d_fe_dv_lsc(1:n_psi_TSC,-nv_lsc_ql_lh,nv_lsc_ql_lh)
c------------------------------------------------------
      implicit none
c      include 'onetwo.i'
      include 'param.i'
      include 'lsc_approach.i'
c-----input
      real*8 
     &v_par, ! parallel velocity normalized to v_t:  m*v_t**2=T 
     &rho    ! small radius. 

c-----local 
      real*8 hrho,coef,h_v,rhol
      integer k,n,n_left,n_right

      hrho=1.d0/n_psi_TSC     ! radial step
      h_v=v_par_lsc_dev_vte_max/(nv_lsc_ql_lh)  ! velocity mesh step 

      rhol=rho

      if(rho.gt.1.d0) then
         d_f_e_nonmaxw_d_v_lsc=0.d0
         goto 10
      endif

      k=idint(rhol/hrho)+1                    ! the number of radial bin
      n=idint(v_par/h_v)                      ! the number of velocity 
                                              ! mesh bin
      if(v_par.ge.0.d0) then
        n_left=n
        n_right=n+1
      else
        n_left=n-1
        n_right=n
      endif

      if(n_left.lt.-nv_lsc_ql_lh) n_left=-nv_lsc_ql_lh
      if(n_right.gt.nv_lsc_ql_lh) n_right=nv_lsc_ql_lh

      if (n_left.eq.n_right) then 
         coef=0.d0
      else
         coef=(v_par-v_par_mesh_lsc(n_left))/
     &     (v_par_mesh_lsc(n_right)-v_par_mesh_lsc(n_left))
      endif

      d_f_e_nonmaxw_d_v_lsc=d_fe_dv_lsc(k,n_left)+
     &      coef*(d_fe_dv_lsc(k,n_right)-d_fe_dv_lsc(k,n_left))

 10   continue

      return
      end

      subroutine calc_RF_power_CD_density_profile
c-----calculte absorbed LH power density
c     and CD_density (without E_DC)
c     CD_density_small_E  
c     radial profiles
c     at TSC uniform radial grid at center bin points

      implicit none 
      include 'param.i'
      include 'lsc_approach.i'
      include 'onetwo.i'
      include 'one_no_nml.i'
c-----locals
      integer i_rho,nv
 
      real*8
     &dense_m3,                !electron density 1/m**3
     &masse_kg,                ! electron mass [kg]
     &charge_e_c,              ! electron charge (coulomb) 
     &masse_g,                 ! electron mass [g]  
     &charge_e,                ! electron charge (statcoulomb)
     &coef_watt_m3,
     &coef_a_m2,
c     &rho,
     &hv,v_te_m_sec,
     &zeff,p,p1,mu,
     &total_power,             !watt 
     &total_current,
     &total_current_small_EDC,           !A
     &total_current_lsc,
     &total_current_East_Karney,
c-----for test     
     &delta_pow_is,power_dens_watt_m3_loc,power_dens_erg_cm3_sec_loc,
     &delta_pow_loc
      integer j
c-----endtest 
c-----externals
      real*8 densrho,zeffrho

      mu=-1.d0 ! =-1 for co-operative current drive
               ! =+1 for current drive into an opposing field
c-----SI
      masse_kg=9.1094d-31      ! electron mass [kg]
      charge_e_c=1.6022d-19    ! electron charge (Coulomb)

c-----CGS
      masse_g=9.1094d-28      ! electron mass [g]     
      charge_e=4.8032d-10     ! electron charge (statcoulomb)

      hv= v_par_mesh_lsc(1)-v_par_mesh_lsc(0)

c-----------------------------------------------------------------
c     arrays of electron dictribution function and it's derivative
c
c     fe_lsc(1:n_psi_TSC,-nv_lsc_ql_lh : nv_lsc_ql_lh)
c     d_fe_dv_lsc(1:n_psi_TSC,-nv_lsc_ql_lh : nv_lsc_ql_lh)
c
c     were calculated before in subroutine calc_fe_lsc_arrays
c--------------------------------------------------------------------
      write(*,*)'in Subroutine calc_RF_power_CD_density_profile'

      do i_rho=1,n_psi_TSC
         rho=rho_bin_center_lsc(i_rho)
         zeff=zeffrho(rho)

         p1=(1.d0+zeff/2.d0+3.d0*mu**2/2.d0)/(3.d0+zeff)

         power_dens_watt_m3_TSC_1D(i_rho)=0.d0
         CD_dens_no_Edc_a_m2_TSC_1D(i_rho)=0.d0  
         CD_dens_small_Edc_a_m2_TSC_1D(i_rho)=0.d0

         do nv = -nv_lsc_ql_lh,nv_lsc_ql_lh 
            p=-v_par_mesh_lsc(nv)*d_ql_lsc_lh_ar(i_rho,nv)*
     &      d_fe_dv_lsc(i_rho,nv)*hv 
            power_dens_watt_m3_TSC_1D(i_rho)=
     &      power_dens_watt_m3_TSC_1D(i_rho)+p


            if((p.lt.0.d0).or.(d_ql_lsc_lh_ar(i_rho,nv).lt.0.d0))then
              write(*,*)'i_rho,nv,p',i_rho,nv,p
              write(*,*)'v_par_mesh_lsc(nv)',v_par_mesh_lsc(nv)
              write(*,*)'d_fe_dv_lsc(i_rho,nv)',d_fe_dv_lsc(i_rho,nv)
           write(*,*)'d_ql_lsc_lh_ar(i_rho,nv)',d_ql_lsc_lh_ar(i_rho,nv)
            endif

            if(nv.ge.1)then
              CD_dens_no_Edc_a_m2_TSC_1D(i_rho)=
     &        CD_dens_no_Edc_a_m2_TSC_1D(i_rho)-
     &        v_par_mesh_lsc(nv)*(fe_lsc(i_rho,nv)-fe_lsc(i_rho,-nv))*hv
            endif

            if (nv.ge.1) then
               mu=1.d0
            else
               mu=-1.d0
            endif 
            
            CD_dens_small_Edc_a_m2_TSC_1D(i_rho)=
     &      CD_dens_small_Edc_a_m2_TSC_1D(i_rho)+
     &      4.d0*v_par_mesh_lsc(nv)**3/(5.d0+zeff)*mu*hv*
     &      d_ql_lsc_lh_ar(i_rho,nv)*d_fe_dv_lsc(i_rho,nv)

         enddo !nv
      
         rho=rho_bin_lsc(i_rho)

c--------transform to SI units
c        power densuty to [watt/m^3]
c        CD density to [a/m**3]

         dense_m3 = densrho(rho,1)*1.d19              !electron density 1/m**3
         v_te_m_sec=v_te_bin_center_lsc(i_rho)*1.d-2  ![m/sec]

         coef_watt_m3 = dense_m3*masse_kg*v_te_m_sec**2/
     &                  tau_n_bin_center_lsc(i_rho) !watt/m**3

         power_dens_watt_m3_TSC_1D(i_rho) = 
     &   power_dens_watt_m3_TSC_1D(i_rho)*coef_watt_m3   !watt/m**3


         coef_a_m2=dense_m3*charge_e_c*v_te_m_sec     !A/m**2

         CD_dens_no_Edc_a_m2_TSC_1D(i_rho)=
     &   CD_dens_no_Edc_a_m2_TSC_1D(i_rho)*coef_a_m2  !A/m**2

         CD_dens_small_Edc_a_m2_TSC_1D(i_rho)=
     &   CD_dens_small_Edc_a_m2_TSC_1D(i_rho)*coef_a_m2  !A/m**2

         call lsc_LH_cur_density(i_rho,j_rf_TSC_1D(i_rho),
     &                           d_j_rf_d_Edc_TSC_1D(i_rho))

         write(*,*)'i_rho,rho_bin_center_lsc(i_rho)',
     &              i_rho,rho_bin_center_lsc(i_rho)

         write(*,*)'power_dens_watt_m3_TSC_1D(i_rho)',
     &              power_dens_watt_m3_TSC_1D(i_rho)

c         write(*,*)'CD_dens_no_Edc_a_m2_TSC_1D(i_rho)',
c     &              CD_dens_no_Edc_a_m2_TSC_1D(i_rho)

         write(*,*)'CD_dens_small_Edc_a_m2_TSC_1D(i_rho)',
     &              CD_dens_small_Edc_a_m2_TSC_1D(i_rho)

         write(*,*)'j_rf_TSC_1D(i_rho)',
     &              j_rf_TSC_1D(i_rho)
         write(*,*)'d_j_rf_d_Edc_TSC_1D(i_rho)',
     &              d_j_rf_d_Edc_TSC_1D(i_rho)

c         write(*,*)'s_cur_den_parallel(i_rho)*1.d4 [a/m**2]',
c     &              s_cur_den_parallel(i_rho)*1.d4
c         write(*,*)'powden_e(i_rho) watt/m3',
c     &              powden_e(i_rho)*1.d-1
        
      enddo !i_rho=1,n_psi_TSC

c-----------------------------------------------------------------
c     write output data to LSC_for_TSC.dat file

      call write_lsc_output
c-----------------------------------------------------------------
      do i_rho =1,n_psi_TSC
           write(*,*)' i_rho',i_rho
           write(*,*)' power_dens_watt_m3_TSC_1D(i_rho)',
     &                 power_dens_watt_m3_TSC_1D(i_rho)
           write(*,*)'CD_dens_small_Edc_a_m2_TSC_1D(i_rho)',
     &                CD_dens_small_Edc_a_m2_TSC_1D(i_rho)
      enddo
       
c------------------------------------------------------------
c     calculate total power and curren from densities

      total_power=0.d0
      total_current_small_EDC=0.d0
      total_current_lsc=0.d0
      do i_rho=1,n_psi_TSC
         total_power=total_power+
     &               power_dens_watt_m3_TSC_1D(i_rho)*
     &               binvol_lsc(i_rho)*1.d-6
         total_current_small_EDC=total_current_small_EDC+
     &                 CD_dens_small_Edc_a_m2_TSC_1D(i_rho)*
     &                 binarea_lsc(i_rho)*1.d-4
         total_current_lsc=total_current_lsc+
     &                 j_rf_TSC_1D(i_rho)*
     &                 binarea_lsc(i_rho)*1.d-4

      enddo 

      write(*,*)'total_power [watt] = ',total_power
      write(*,*)'total_current_small_EDC [A] =',total_current_small_EDC
      write(*,*)'total_current_lsc [A] =',total_current_lsc
      write(*,*)'total_current_small_EDC*powtott*1.d-7/total_power =',
     &           total_current_small_EDC*powtott*1.d-7/total_power

      total_power=0.d0
      total_current_East_Karney=0.d0

      do i_rho=1,NR-1
         total_power=total_power+
     &               powden_e(i_rho)*
     &               binvol(i_rho)*1.d-6*1.d-1
         total_current_East_Karney=total_current_East_Karney+
     &                 s_cur_den_parallel(i_rho)*
     &                 binarea(i_rho)        
      enddo 

      write(*,*)'2 total_power [watt] = ',total_power
      write(*,*)'s_cur_den_parallel total_current [A] =',
     &total_current_East_Karney

      do i_rho=1,NR-1
      write(*,*)'i_rho,powden_e(i_rho)',i_rho,powden_e(i_rho)
      enddo

      total_current=0.d0
      do i_rho=1,NR-1
        
         total_current=total_current+
     &                 currden(i_rho)*
     &                 binarea(i_rho)
      enddo 
      write(*,*)'currden total_current [A] =',total_current
c----test------------------------------------------------------------
      delta_pow_is=0.d0 
      do j=1,n_psi_TSC
         power_dens_watt_m3_loc=0.d0
         do nv = -nv_lsc_ql_lh,nv_lsc_ql_lh
            power_dens_watt_m3_loc=
     &      power_dens_watt_m3_loc-
     &      v_par_mesh_lsc(nv)*d_ql_lsc_lh_ar(j,nv)*
     &      d_fe_dv_lsc(j,nv)*hv
         enddo !nv
         dense_m3 = densrho(rho_bin_center_lsc(j),1)*1.d19  !electron density
                                                            ! 1/m**3
         v_te_m_sec=v_te_bin_center_lsc(j)*1.d-2  ![m/sec]
         coef_watt_m3 = dense_m3*masse_kg*v_te_m_sec**2/
     &                  tau_n_bin_center_lsc(j) !watt/m**3
         power_dens_watt_m3_loc=
     &          power_dens_watt_m3_loc*coef_watt_m3   !watt/m**3 
         power_dens_erg_cm3_sec_loc=
     &          power_dens_watt_m3_loc*10.d0   !erg/cm**3/sec  
         delta_pow_loc=power_dens_erg_cm3_sec_loc* binvol_lsc(j)
         delta_pow_is=delta_pow_is+delta_pow_loc
      enddo !j
 
      write(*,*)'total_power=delta_pow_is',delta_pow_is

c----------------------------------------
      return     
      end

      subroutine lsc_LH_cur_density(i_rho,LH_cd_density,
     &d_LH_cur_dens_d_E_DC )
c-------------------------------------------------
c     calculate parallel to the magnetic field
c     LH current drive density LH_cd_density [A/m**2]
c     and its derivative by E_DC:d_LH_cur_dens_d_E_DC [A/(V*m)]
c     using LSC approach
c     at i_rho_TSC number of central radial bin point at TSC mesh
c-------------------------------------------------

      implicit none

      include 'param.i'
      include 'lsc_approach.i'
      include 'onetwo.i'

c-----input
      integer i_rho ! number  of the point at
                    ! central of radial bins mesh   
c-----output
      real*8
     &LH_cd_density,       ! LH current drive density [A/m**2]
     &d_LH_cur_dens_d_E_DC ! derivative d_LH_cur_dens/d_E_DC [A/(V*m)]
c-----local
      real*8
     &v_r_m_s,        ! V_r [m/sec]
     &v_te_m_s,       ! thermal velocity [m/sec]
     &v_r_d_vt,       !normalized V_r/V_te
     &rho,            ! normalized small radius 
     &h_v,            ! normaized parallel velocity mesh step
     &v_par,          ! normaized parallel velocity
     &masse_kg,       ! electron mass [kg]
     &charge_e_c,     !  electron charge [Coulomb]   
     &tau_n,          ! character time = v_te**3/Gamma
     &dens_e_m3,      ! electron density 1/m**3
     &E_DC_V_m,       !electric field [V/m]
     &delta_E,        ! nondimensional
     &E_DC_V_m_minus, !E_DC_V_m*(1-delta_E)
     &E_DC_V_m_plus,  !E_DC_V_m*(1+delta_E)
     &v_r_m_s_minus,  !V_r_minus from E_DC_V_m_minus [m/sec]
     &v_r_m_s_plus,   !V_r_plus from E_DC_V_m_plus  [m/sec]
     &v_r_d_vt_minus, !normalized V_r_minus/V_te
     &v_r_d_vt_plus,  !normalized V_r_plus/V_te
     &LH_cd_density_minus, !LH current drive density [A/m**2] for E_DC_V_m_minus
     &LH_cd_density_plus,  !LH current drive density [A/m**2] for E_DC_V_m_plus
c-----input arguments for WsloPrmZ
     & uGiven ,  !v_par/v_r
     & MuGiven , !u_par/u
                 !u_par=v_par/v_r
                 !u=v/|v_r|
     & ZGiven,
     & uGiven_minus, !v_par/v_r_minus
     & uGiven_plus,  !v_par/v_r_plus
c-----output from  WsloPrmZ
     &dWsduou,          !=(dW_s/d_u)/u
     &dWsduou_minus,    !=(dW_s/d_u)/u for E_DC_minus
     &dWsduou_plus      !=(dW_s/d_u)/u for E_DC_plus
      integer nv,    
     &iWhichWay  !output from  WsloPrmZ

c-----externals
      real*8 d_f_e_nonmaxw_d_v_lsc,zeffrho,densrho,temperho

c--------------------------------------------------------------------
c     V_r = -sign(e*E_DC)*sqrt(m_e*Gamma/|eE_DC|)
c 
c     for electron e<0
c
c     tau_n=v_te^3/Gamma
c       
c     So
c     V_r = sign(E_DC)*v_te*sqrt(m_e*v_te/(tau_n*|eE_DC|))
c
c     for positive E_DC: V_r >0


      delta_E=1.d-5

      rho=rho_bin_center_lsc(i_rho)

c      write(*,*)'in lsc_LH_cur_density i_rho,rho',i_rho,rho

      E_DC_V_m=EdcTSC_1D(i_rho)   ! E_DC [V/m]

      v_te_m_s = v_te_bin_center_lsc(i_rho)*1.d-2 !thermal velocity [m/sec]
      tau_n=tau_n_bin_center_lsc(i_rho)
 
c      write(*,*)'in lsc_LH_cur_density i_rho,rho,EdcTSC_1D(i_rho)',
c     & i_rho,rho,EdcTSC_1D(i_rho)

      masse_kg=9.1094d-31      ! electron mass [kg]
      charge_e_c=1.6022d-19    ! electron charge (coulomb)
  
      v_r_m_s=dsign(1.d0,E_DC_V_m)*v_te_m_s*
     &dsqrt(masse_kg*v_te_m_s/(tau_n*dabs(charge_e_c*E_DC_V_m ))) ! [m/sec]
     
      v_r_d_vt = v_r_m_s/v_te_m_s

c      write(*,*)'in lsc_LH_cur_densit v_r_m_s,v_te_m_s,v_r_d_vt',
c     &                                v_r_m_s,v_te_m_s,v_r_d_vt 

      LH_cd_density=0.d0
c-----for derivative
      E_DC_V_m_minus=E_DC_V_m*(1-delta_E)
      E_DC_V_m_plus=E_DC_V_m*(1+delta_E)
c      v_r_m_s_minus=dsign(1.d0,E_DC_V_m_minus)*v_te_m_s*
c     &dsqrt(masse_kg*v_te_m_s/(tau_n*dabs(charge_e_c*E_DC_V_m_minus))) ! [m/sec]
c      write(*,*)'1 v_r_m_s_minus',v_r_m_s_minus
      v_r_m_s_minus=v_r_m_s/dsqrt(1-delta_E)
c      write(*,*)'2 v_r_m_s_minus ',v_r_m_s_minus

c      v_r_m_s_plus=dsign(1.d0,E_DC_V_m_plus)*v_te_m_s*
c     &dsqrt(masse_kg*v_te_m_s/(tau_n*dabs(charge_e_c*E_DC_V_m_plus))) ! [m/sec]
c      write(*,*)'1 v_r_m_s_plus',v_r_m_s_plus
      v_r_m_s_plus=v_r_m_s/dsqrt(1+delta_E)
c      write(*,*)'2 v_r_m_s_plus ',v_r_m_s_plus

      v_r_d_vt_minus = v_r_m_s_minus/v_te_m_s
      v_r_d_vt_plus = v_r_m_s_plus/v_te_m_s

      LH_cd_density_minus=0.d0
      LH_cd_density_plus=0.d0
c--------------------------------------------------------
      dens_e_m3=densrho(rho,1)*1.d19 !electron density [1/m**3]

      ZGiven=zeffrho(rho)     
       
c      write(*,*)'in lsc_LH_cur_densit dens_e_m3,ZGiven',dens_e_m3,ZGiven

      h_v=v_par_mesh_lsc(1)-v_par_mesh_lsc(0)

      do nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
c---------------------------------------------------------------
c        calculate derivative from W_s function
c        dWsduou=(dW_s/d_u)/u
c----------------------------------------------------------------
         v_par=v_par_mesh_lsc(nv)
         uGiven=v_par_mesh_lsc(nv)/v_r_d_vt

c         write(*,*)'nv,v_par,v_r_d_vt,uGiven',nv,v_par,v_r_d_vt,uGiven

         if (uGiven .ge.0.d0)  then
            MuGiven=1.d0
         else
            MuGiven=-1.d0
         endif
c-------------------------------------------------------------------
c        Provides for a linear interpolation of Zeff
c        for dW/du
c        in the region 1.0 <Zeff< 10.
         call WsloPrmZ( uGiven , MuGiven , ZGiven,
     &                  dWsduou, iWhichWay) 

c         write(*,*)'dWsduou, iWhichWay',dWsduou, iWhichWay    
         LH_cd_density = LH_cd_density + 
     &   d_ql_lsc_lh_ar(i_rho,nv)*
     &   d_f_e_nonmaxw_d_v_lsc(v_par,rho)*
     &   dWsduou*uGiven*h_v 
c--------for derivative------------------------------

         uGiven_minus=v_par_mesh_lsc(nv)/v_r_d_vt_minus
         uGiven_plus=v_par_mesh_lsc(nv)/v_r_d_vt_plus

         call WsloPrmZ( uGiven_minus , MuGiven , ZGiven,
     &                  dWsduou_minus, iWhichWay) 
         call WsloPrmZ( uGiven_plus , MuGiven , ZGiven,
     &                  dWsduou_plus, iWhichWay) 

         LH_cd_density_minus = LH_cd_density_minus  + 
     &   d_ql_lsc_lh_ar(i_rho,nv)*
     &   d_f_e_nonmaxw_d_v_lsc(v_par,rho)*
     &   dWsduou_minus*uGiven_minus*h_v 

         LH_cd_density_plus = LH_cd_density_plus  + 
     &   d_ql_lsc_lh_ar(i_rho,nv)*
     &   d_f_e_nonmaxw_d_v_lsc(v_par,rho)*
     &   dWsduou_plus*uGiven_plus*h_v 
c------------------------------------------------


c         write(*,*)'d_ql_lsc_lh_ar(i_rho,nv)',d_ql_lsc_lh_ar(i_rho,nv)
c         write(*,*)'d_f_e_nonmaxw_d_v_lsc(v_par,rho)',
c     &              d_f_e_nonmaxw_d_v_lsc(v_par,rho)
c         write(*,*)'d_fe_dv_lsc(i_rho,nv)',d_fe_dv_lsc(i_rho,nv)
c         write(*,*)'LH_cd_density', LH_cd_density 
      enddo

      LH_cd_density=LH_cd_density*charge_e_c*dens_e_m3*
     &v_te_m_s*v_r_d_vt**3        ![A/m**2]

c-----for derivative-------------------------------------                     
      LH_cd_density_minus=LH_cd_density_minus*charge_e_c*dens_e_m3*
     &v_te_m_s*v_r_d_vt_minus**3           ![A/m**2]

      LH_cd_density_plus=LH_cd_density_plus*charge_e_c*dens_e_m3*
     &v_te_m_s*v_r_d_vt_plus**3            ![A/m**2]

      d_LH_cur_dens_d_E_DC=(LH_cd_density_plus-LH_cd_density_minus)/
     &                     (E_DC_V_m*2.d0*delta_E)

      return
      end


      SUBROUTINE p_c_prof_lsc(is,rhobegin,rhoend,cnpar,cnper,
     &cefldx,cefldy,cefldz)
c----------------------------------------------------------------
c     this subroutine calculates absorted power profiles
c     calculates RF current  profiles	 (from deposited power)
c-----------------------------------------------------------------
      !implicit double precision (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
      INCLUDE 'three.i'
      INCLUDE 'onetwo.i'
      INCLUDE 'write.i'
      INCLUDE 'gr.i'
      INCLUDE 'rho.i'
     
c-----input
      real*8 rhobegin,rhoend,cnpar,cnper
      integer is
      complex*16 cefldx,cefldy,cefldz !polarization

c-----locals
      real*8 delpow_s,ppow_s
      dimension delpow_s(nbulka) !Added for indiv ion contrib[BH041009]
      dimension ppow_s(nbulka)  !Added for indiv ion contrib[BH041009]

      real*8 z_r,r_r,ymar,xmar,z_effr,tempr,denr,cnparr,yer,effic_r

      real*8 r_m,z_m,temp,ye,u1,den,z_eff,hro,rholeft,rhoright,
     &rhomax,rhomin,eff_rho_max,eff_rho_min,hrho,delpower,delpow_i,
     &delpow_e,delpow_cl,del,r0,z0,r0m,z0m,psiloc,zfacgeom,aspct,
     &delcurr,rho0,poloidlen,rho0_pol,delrho,
     &ppow,pcur,ppow_e,ppow_i,ppow_cl,
     &r_bmin_right,r_bmin_left,z_bmin_right,z_bmin_left,theta,psi

      integer i,kk,jbinmin,jbinmax,j,k
cSAP070831
      real*8 psi_loc,cos_theta_pol,sin_theta_pol,rho_loc,theta_pol
      integer n_theta_pol,ir
c-----external
      real*8  tempe,y,u_res,dense,zeff,psi_rho,b,efficien,rhov,rhos,
     &qsafety_psi,bmin_psi,bmax_psi,dvol_dpsi,rho_lrho,
     &b_average

      pi=4.d0*datan(1.d0)

c     wr and wz are in (cm) 
      r_m=wr(is)*0.01d0/r0x ! normalization
      z_m=wz(is)*0.01d0/r0x ! normalization

c       write(*,*)'p_c_prof NR,is,ieffic',NR,is,ieffic

c  radius rho is in common/one/,it was calculated in subroutine: b
c  rho is used in functions:tempe,dense,z_eff
c  The magnetic field bmode is in common/one/ .bmode is used
c  in function y(z,r,phi)


c begin if is=1      

c      write(*,*)'prep3d.f p_c_prof_ is,rhobegin,rhoend,rho,cnpar,cnper',
c     &                             is,rhobegin,rhoend,rho,cnpar,cnper

      if(is.eq.1) then

c         write(*,*)'prep3d.f p_c_prof is=1 rho=',rho
cSAP080731
        if (rho.ge.1.d0) then
c---------------------------------------------------------------
c          zero CD efficiency outside the plasma
c----------------------------------------------------------------
           eff(is)=0.d0
           goto 100
        endif
  


c        calculation of efficiency on the first step
         temp=tempe(z_m,r_m,wphi(is),1)
         ye=y(z_m,r_m,wphi(is),1)
         u1=u_res(jwave,cnpar,temp,ye)
         den=dense(z_m,r_m,wphi(is),1)
         z_eff=zeff(z_m,r_m,wphi(is))
         if (ieffic.eq.1) then
c -------------------------------------------------------------------
c     calculation of current drive efficiency using asymptotic formula
c--------------------------------------------------------------------
          eff(is)=efficien(z_eff,u1,jwave,temp,den) 

c         write(*,*)'prep3d asymptotic eff',eff(is)
         endif

         if (ieffic.eq.2) then
c -------------------------------------------------------------------
c     calculation of current drive efficiency using Ehst-Karney formula
c--------------------------------------------------------------------
c          write(*,*)'prep3d bef Karn z_m,r_m,r0x,z_eff,temp,den,jwave',
c     &    z_m,r_m,r0x,z_eff,temp,den,jwave

cfor adj
c          rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2)
c          cos_theta_pol=(r_m-xma)/rho_loc
c          sin_theta_pol=(z_m-yma)/rho_loc
c          if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
c          if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
c          if (sin_theta_pol.ge.0.d0) then
c             theta_pol=+dacos(cos_theta_pol)
c          else  
c             theta_pol=-dacos(cos_theta_pol)
c          endif

c          if (theta_pol.lt.0.d0) then
c             theta_pol=theta_pol+2.d0*pi
c          endif

c          if (theta_pol.gt.2.d0*pi) then
c             n_theta_pol=theta_pol/(2.d0*pi)
c             theta_pol=theta_pol-2.d0*pi*n_theta_pol
c          endif
         
  
c          psi_loc=psi_rho(rho)

c          write(*,*)'prep3d.f in p_c_prof ,cefldy,cefldz',cefldy,cefldz
c          write(*,*)'prep3d.f in p_c_prof cnpar,cnper',cnpar,cnper
c          write(*,*)'prep3d.f in p_c_prof psi_loc,theta_pol',
c     &    psi_loc,theta_pol

c          call CD_adj_LH_efficiency(cnpar,cnper,
c     &    psi_loc,theta_pol,cefldy,cefldz,
c     &    eff(is))
c          write(*,*)'prep3d ADJ is,cnpar,eff(is)',is,cnpar,eff(is)
c-------------------------------------------------------------------
c           write(*,*)'prep3d.f  in p_c_prof before'//
c     &               'asimptotic_power_CD_density_LSC'
 
c           call asimptotic_power_CD_density_LSC(
c     &     z_m,r_m,0.d0,cnpar,cnper)
c     &     CD_dens,power_dens,CD_efficiency)
c           write(*,*)'prep3d.f  in p_c_prof after'//
c     &               'asimptotic_power_CD_density_LSC'

c           write(*,*)'prep3d z_eff',z_eff
           call efKarney(z_m*100.d0*r0x,r_m*100.0d0*r0x,r0x,
     +                 z_eff,temp,den,jwave,
     1                 cnpar,eff(is))
c           write(*,*)'prep3d Karney is,cnpar,eff(is)',is,cnpar,eff(is)
c----------------------------------------------------
c          ef_Karney using Bonoli subroutine:
c          For the FW NSTX case the Bonoli subroutine gave efficiency
c          which coincided exactly with the efficiency obtained 
c          calculated bu subroutine call efKarney
c
c-------------------------------------------------------
c          call efKarney_Bonoli(z_m*100.d0*r0x,r_m*100.0d0*r0x,r0x,
c     +                 z_eff,temp,den,jwave,
c     1                 cnpar,eff(is))
c          write(*,*)'prep3d Karney_Bonoli is,cnpar,eff(is)',
c     &                is,cnpar,eff(is)
c------------------------------------------------------------
	 endif   ! ieffic.eq.2

         if (ieffic.eq.3) then
c -------------------------------------------------------------------
c          calculation of current drive efficiency using curba
c--------------------------------------------------------------------
           z_r=(wz(is))
           r_r=(wr(is))
           ymar=(yma*100.d0*r0x) !cm
           xmar=(xma*100.d0*r0x) !cm
           z_effr=(z_eff)
           tempr=(temp)
           denr=(den)
           cnparr=(cnpar)
           yer=(ye)

c          write(*,*)'in p_c_prof before effcurb is=1 z_r,r_r,ymar,xmar',
c     .    z_r,r_r,ymar,xmar
c          write(*,*)'z_effr,tempr,denr,cnparr,yer',
c     .    z_effr,tempr,denr,cnparr,yer
c          write(*,*)'in p_c_prof before effcurb'

      call effcurb(z_r,r_r,ymar,xmar,r0x,z_effr,tempr,denr,jwave,cnparr,
     +             yer,effic_r)
           write(*,*)'in p_c_prof after effcurb efffic_r',effic_r
           eff(is)=(effic_r)
         endif     ! ieffic.eq.3

         if (ieffic.eq.4) then
c -------------------------------------------------------------------
c          calculation of current drive efficiency using Lin_liu
c          TorGA_curgap
c--------------------------------------------------------------------
           z_r=(wz(is))
           r_r=(wr(is))
           ymar=(yma*100.d0*r0x) !cm
           xmar=(xma*100.d0*r0x) !cm
           z_effr=(z_eff)
           tempr=(temp)
           denr=(den)
           cnparr=(cnpar)
           yer=(ye)

           call eff_Lin_Liu(z_r,r_r,ymar,xmar,r0x,z_effr,tempr,denr,
     &                 jwave,
     &                 cnparr,
     &                 cnper,ioxm,ye,
     &                 cefldx,cefldy,cefldz,
     &                 effic_r)

           write(*,*)'in p_c_prof after eff_Lin_Liu efffic_r',effic_r

           eff(is)=(effic_r)
         endif    ! ieffic.eq.4

         if (ieffic.eq.5) then
c -------------------------------------------------------------------
c        calculation of current drive efficiency using adj function 
c        for n_harmonic=0 case (LH wave)
c--------------------------------------------------------------------
          rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2)
          
          if(rho_loc.eq.0.d0)then ! YuP [06,2015] To avoid 1/rho_loc=Inf
          cos_theta_pol=1.d0
          sin_theta_pol=0.d0
          else
          cos_theta_pol=(r_m-xma)/rho_loc
          sin_theta_pol=(z_m-yma)/rho_loc
          endif
          
          if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
          if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
          if (sin_theta_pol.ge.0.d0) then
             theta_pol=+dacos(cos_theta_pol)
          else  
             theta_pol=-dacos(cos_theta_pol)
          endif

          if (theta_pol.lt.0.d0) then
             theta_pol=theta_pol+2.d0*pi
          endif

          if (theta_pol.gt.2.d0*pi) then
             n_theta_pol=theta_pol/(2.d0*pi)
             theta_pol=theta_pol-2.d0*pi*n_theta_pol
          endif
           
          psi_loc=psi_rho(rho)

c          write(*,*)'prep3d.f in p_c_prof ,cefldy,cefldz',cefldy,cefldz
c          write(*,*)'prep3d.f in p_c_prof cnpar,cnper',cnpar,cnper
c          write(*,*)'prep3d.f in p_c_prof psi_loc,theta_pol',
c     &    psi_loc,theta_pol

          call CD_adj_LH_efficiency(cnpar,cnper,
     &    psi_loc,theta_pol,cefldy,cefldz,
     &    eff(is))
          write(*,*)'prep3d ADJ is,cnpar,eff(is)',is,cnpar,eff(is)
         endif !ieffic.eq.5


         if (ieffic.eq.6) then
c -------------------------------------------------------------------
c        calculation of current drive efficiency using adj function 
c        for all harmoniucs general case
c--------------------------------------------------------------------
          rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2)
          
          if(rho_loc.eq.0.d0)then ! YuP [06,2015] To avoid 1/rho_loc=Inf
          cos_theta_pol=1.d0
          sin_theta_pol=0.d0
          else
          cos_theta_pol=(r_m-xma)/rho_loc
          sin_theta_pol=(z_m-yma)/rho_loc
          endif
      
          if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
          if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
          if (sin_theta_pol.ge.0.d0) then
             theta_pol=+dacos(cos_theta_pol)
          else  
             theta_pol=-dacos(cos_theta_pol)
          endif

          if (theta_pol.lt.0.d0) then
             theta_pol=theta_pol+2.d0*pi
          endif

          if (theta_pol.gt.2.d0*pi) then
             n_theta_pol=theta_pol/(2.d0*pi)
             theta_pol=theta_pol-2.d0*pi*n_theta_pol
          endif
           
          psi_loc=psi_rho(rho)

c          write(*,*)'prep3d.f in p_c_prof cefldx,cefldy,cefldz',
c     &                                    cefldx,cefldy,cefldz
c          write(*,*)'prep3d.f in p_c_prof cnpar,cnper',cnpar,cnper
          write(*,*)'prep3d.f in p_c_prof psi_loc,theta_pol',
     &    psi_loc,theta_pol

      
          call CD_adj_efficiency(cnpar,cnper,
     &    psi_loc,theta_pol,cefldx,cefldy,cefldz,
     &    eff(is))

          write(*,*)'prep3d ADJ is,cnpar,eff(is)',is,cnpar,eff(is)
         endif !ieffic.eq.6
         
          
c--------------------------------------------------------------------
c        toroidal and poloidal current drives efficiencies for is=1
 100     continue
         bmod=b(z_m,r_m,0.d0)        
c--------------------------------------------------------------------
         allpower=0.0d0
         allpw_e=0.0d0
         allpw_i=0.0d0
         allpw_cl=0.0d0

         allcur=0.0d0
        
c------- initialization arrays
c        for power and current
         do i=1,NR
	    power(i)=0.0d0
	    current(i)=0.0d0           
	    power_e(i)=0.0d0
	    power_i(i)=0.0d0
	    power_cl(i)=0.0d0
            cur_den_parallel(i)=0.d0 
	 enddo

c--------binvol(NR) calculations
         theta=0.d0 
         hro=1.d0/dble(NR-1)
         do i=1,NR-1 !for onetwo.i in [cm**3]
            rholeft=hro*(i-1)
            rhoright=rholeft+hro
            binvol(i)=voltot*(rhov(rhoright)**2-rhov(rholeft)**2)*1.d6

c            write(*,*)'NR,i,binvol(i)',NR,i,binvol(i)

           binarea(i)=areatot*(rhos(rhoright)**2-rhos(rholeft)**2)*1.d4

            psi=psi_rho(rholeft)
            call zr_psith(psi,theta,z_bmin_left,r_bmin_left)
            psi=psi_rho(rhoright)
            call zr_psith(psi,theta,z_bmin_left,r_bmin_right)

            binarea_pol(i)=pi*(r_bmin_right-r_bmin_left)*
     &                        (r_bmin_right+r_bmin_left)*1.d4

         enddo
 
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               do i=1,NR
                  power_s(i,kk)=0.0d0
               enddo
            enddo
         endif

	 goto 30
      endif
c end if is=1

c      write(*,*)' p_c_prof: rhobegin,rhoend',rhobegin,rhoend

      if(rhoend.gt.rhobegin) then
         rhomax=rhoend
         rhomin=rhobegin
      else
	 rhomin=rhoend
	 rhomax=rhobegin
      endif
cSAP080831
      if(rhomax.gt.1.d0) rhomax=1.d0
      if(rhomin.gt.1.d0) rhomin=1.d0-1.d-10    

      hrho=1.0d0/(NR-1)
      jbinmin=1
      jbinmax=NR-1

c      write(*,*)'NR,hrho,rhomin,rhomax',NR,hrho,rhomin,rhomax

      do j=1,NR-1
         if(rhomin.lt.(hrho*j)) then
           jbinmin=j
	    goto 10
	 endif
      enddo
10    continue
      do j=jbinmin,NR-1
         if(rhomax.lt.(hrho*j)) then
            jbinmax=j
	    goto 20
	 endif
      enddo
20    continue

      if (jbinmin.eq.NR) jbinmin=NR-1

c      write(*,*)'p_c_prof is,jbinmax&min',is,jbinmax,jbinmin
     
c-----------------------------------------------------------
c     here delpower and allpower are in (erg/sec)
c------------------------------------------------------------
      delpower=delpwr(is-1)-delpwr(is)

cSmirnov970105 beg
cBH001017      delpow_e=0.5d0*(ws(is)-ws(is-1))*
cBH001017     1(delpwr(is-1)*sdpwr(is-1)+delpwr(is)*sdpwr(is))
cBH001017      delpow_i=0.5d0*(ws(is)-ws(is-1))*
cBH001017     1(delpwr(is-1)*salphal(is-1)+delpwr(is)*salphal(is))
cHarvey991017 beg
      delpow_i=0.5d0*(ws(is)-ws(is-1))*
     1(delpwr(is-1)*sdpwr(is-1)+delpwr(is)*sdpwr(is))

      if (iabsorp.eq.3) then
         do kk=2,nbulk
            delpow_s(kk)=0.5d0*(ws(is)-ws(is-1))*
     1        (delpwr(is-1)*salphas(is-1,kk)+delpwr(is)*salphas(is,kk))
         enddo
      endif
cSAP100528      
c      delpow_e=0.5d0*(ws(is)-ws(is-1))*
c     1(delpwr(is-1)*salphal(is-1)+delpwr(is)*salphal(isc))  
      delpow_e=(ws(is)-ws(is-1))*
     1(delpwr(is)*salphal(is))  

c      write(*,*)'p_c_prof is,delpwr(is-1),delpwr(is),delpower,delpow_e'
c     &,is,delpwr(is-1),delpwr(is),delpower,delpow_e

cHarvey991017 end
      delpow_cl=0.5d0*(ws(is)-ws(is-1))*
     1(delpwr(is-1)*salphac(is-1)+delpwr(is)*salphac(is))

      del=delpow_i+delpow_e+delpow_cl
cHarvey970111 beg
cSmirnov050215      if (del.ne.0.d0) then
cSmirnov050215      Change 0.d0 to a small number
      if (del.gt.1.d-100) then
         delpow_e=delpow_e*delpower/del
         delpow_i=delpow_i*delpower/del
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               delpow_s(kk)=delpow_s(kk)*delpower/del
            enddo
         endif
         delpow_cl=delpow_cl*delpower/del

      endif
cSmirnov970105 end
      allpower=allpower+delpower
cSAP091021
c      write(*,*)'delpower,delpow_i,delpow_e',delpower,delpow_i,delpow_e
c      write(*,*)'delpow_s(k),k=2,nbulk',(delpow_s(k),k=2,nbulk)
c      write(*,*)'allpower,delpower',allpower,delpower
cSmirnov970106 beg
      allpw_e=allpw_e+delpow_e
      allpw_i=allpw_i+delpow_i
      allpw_cl=allpw_cl+delpow_cl
      
cSmirnov970106 end
c     write(*,*)' allpower,allpw_e,allpw_i,allpw_cl',
c    1 allpower,allpw_e,allpw_i,allpw_cl
c     write(*,*)'delpower,delpow_e,delpow_cl',
c    1delpower,delpow_e,delpow_cl
c-----------------------------------------------------------
c     r0 (cm)
      r0=0.5d0*(wr(is)+wr(is-1))

cSAP080831
c      write(*,*)'p_c_prof is>1  rho',rho
      if (rho.ge.1.d0) then
c---------------------------------------------------------------
c        zero CD efficiency outside the plasma
c----------------------------------------------------------------
         eff(is)=0.d0
         goto 110
      endif
c---- calculation of the efficiency 
      temp=tempe(z_m,r_m,wphi(is),1)
      ye=y(z_m,r_m,wphi(is),1)
      u1=u_res(jwave,cnpar,temp,ye)
      den=dense(z_m,r_m,wphi(is),1)
      z_eff=zeff(z_m,r_m,wphi(is))
   
      if (ieffic.eq.1) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using asymptotic formula
c--------------------------------------------------------------------
        eff(is)=efficien(z_eff,u1,jwave,temp,den)
c       write(*,*)'prep3d asymptotic eff',eff(is)
      endif

      if (ieffic.eq.2) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using Ehst-Karney formula
c--------------------------------------------------------------------
c      write(*,*)'prep3d bef Karn z_m,r_m,r0x,z_eff,temp,den,jwave',
c     &    z_m,r_m,r0x,z_eff,temp,den,jwave

cfor adj
c          rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2)
c          cos_theta_pol=(r_m-xma)/rho_loc
c          sin_theta_pol=(z_m-yma)/rho_loc
c          if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
c          if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
c          if (sin_theta_pol.ge.0.d0) then
c             theta_pol=+dacos(cos_theta_pol)
c          else  
c             theta_pol=-dacos(cos_theta_pol)
c          endif

c          if (theta_pol.lt.0.d0) then
c             theta_pol=theta_pol+2.d0*pi
c          endif

c          if (theta_pol.gt.2.d0*pi) then
c             n_theta_pol=theta_pol/(2.d0*pi)
c             theta_pol=theta_pol-2.d0*pi*n_theta_pol
c          endif
         
c          psi_loc=psi_rho(rho)
c          write(*,*)'prep3d.f in p_c_prof cefldy,cefldz',cefldy,cefldz
c          write(*,*)'prep3d.f in p_c_prof cnpar,cnper',cnpar,cnper
c          write(*,*)'prep3d.f in p_c_prof psi_loc,theta_pol',
c     &    psi_loc,theta_pol

c          call CD_adj_LH_efficiency(cnpar,cnper,
c     &    psi_loc,theta_pol,cefldy,cefldz,
c     &    eff(is))
c          write(*,*)'prep3d ADJ is,cnpar,eff(is)',is,cnpar,eff(is)
c-------------------------------------------------------------------
c      write(*,*)'prep3d z_eff',z_eff   
c           write(*,*)'prep3d.f  in p_c_prof before'//
c     &               'asimptotic_power_CD_density_LSC'
 
c           call asimptotic_power_CD_density_LSC(
c     &     z_m,r_m,0.d0,cnpar,cnper)
c     &     CD_dens,power_dens,CD_efficiency)
c           write(*,*)'prep3d.f  in p_c_prof after'//
c     &               'asimptotic_power_CD_density_LSC'

      call efKarney(z_m*100.0d0*r0x,r_m*100.0d0*r0x,r0x,
     & z_eff,temp,den,jwave,cnpar,
     &              eff(is))
c      write(*,*)'prep3d Karney is,cnpar,eff(is)',is,cnpar,eff(is)
c----------------------------------------------------
c          ef_Karney using Bonoli subroutine:
c          For the FW NSTX case the Bonoli subroutine gave efficiency
c          which coinsided exectly wiht the efficiency obtained 
c          calculated bu subroutine call efKarney
c
c-------------------------------------------------------
c      call efKarney_Bonoli(z_m*100.0d0*r0x,r_m*100.0d0*r0x,r0x,
c     & z_eff,temp,den,jwave,cnpar,
c     &              eff(is))
c      write(*,*)'prep3d Karney_Bonoli is,cnpar,eff(is)',is,cnpar,eff(is)

      endif  !ieffic.eq.2

      if (ieffic.eq.3) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using curba
c-------------------------------------------------------------------
        z_r=(wz(is))
        r_r=(wr(is))
        ymar=(yma*100.d0*r0x)
        xmar=(xma*100.d0*r0x)
        z_effr=(z_eff)
        tempr=(temp)
        denr=(den)
        cnparr=(cnpar)
        yer=real(ye)


c         write(*,*)'in p_c_prof before effcurbis.ne.1 z_r,r_r,
c     .   ymar,xmar,r0x,jwave',
c     .   z_r,r_r,ymar,xmar,r0x,jwave
c         write(*,*)'z_effr,tempr,denr,cnparr,yer',
c     .   z_effr,tempr,denr,cnparr,yer
c         write(*,*)'in p_c_prof before effcurb' 

ctest
        u1=u_res(jwave,cnparr,tempr,yer)
        write(*,*)'jwave,cnpar,tempr,yer,u1',jwave,cnparr,tempr,yer,u1
        effic_r=efficien(z_effr,u1,jwave,tempr,denr) ! here: just for comparison
        write(*,*)'asimptotic: z_effr,denr,effic_r',z_effr,denr,effic_r
cendtest
        call effcurb(z_r,r_r,ymar,xmar,r0x,z_effr,tempr,denr,jwave,
     +  cnparr,yer,effic_r) ! actual effic_r

        write(*,*)'in p_c_prof after effcurb is,effic_r',is,effic_r
 
        eff(is)=effic_r
      endif !ieffic.eq.3

      if (ieffic.eq.4) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using eff_Lin_Liu
c-------------------------------------------------------------------
        z_r=(wz(is))
        r_r=(wr(is))
        ymar=(yma*100.d0*r0x)
        xmar=(xma*100.d0*r0x)
        z_effr=(z_eff)
        tempr=(temp)
        denr=(den)
        cnparr=(cnpar)
        yer=real(ye)


c         write(*,*)'in p_c_prof before effcurbis.ne.1 z_r,r_r,
c     .   ymar,xmar,r0x,jwave',
c     .   z_r,r_r,ymar,xmar,r0x,jwave
c         write(*,*)'z_effr,tempr,denr,cnparr,yer',
c     .   z_effr,tempr,denr,cnparr,yer
c         write(*,*)'in p_c_prof before effcurb'

        call effcurb(z_r,r_r,ymar,xmar,r0x,z_effr,tempr,denr,jwave,
     +  cnparr,yer,effic_r)

        write(*,*)'in p_c_prof after effcurb is,effic_r',is,effic_r

           call eff_Lin_Liu(z_r,r_r,ymar,xmar,r0x,z_effr,tempr,denr,
     &                 jwave,
     &                 cnparr,
     &                 cnper,ioxm,ye,
     &                 cefldx,cefldy,cefldz,
     &                 effic_r)
        write(*,*)'in p_c_prof after eff_Lin_Liu efffic_r',effic_r
   
        eff(is)=effic_r
      endif !4

      if (ieffic.eq.5) then
c -------------------------------------------------------------------
c        calculation of current drive efficiency using adj function 
c        for n_harmonic=0 case (LH wave)
c--------------------------------------------------------------------
          rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2)
          
          if(rho_loc.eq.0.d0)then ! YuP [06,2015] To avoid 1/rho_loc=Inf
          cos_theta_pol=1.d0
          sin_theta_pol=0.d0
          else
          cos_theta_pol=(r_m-xma)/rho_loc
          sin_theta_pol=(z_m-yma)/rho_loc
          endif

          if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
          if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
          if (sin_theta_pol.ge.0.d0) then
             theta_pol=+dacos(cos_theta_pol)
          else  
             theta_pol=-dacos(cos_theta_pol)
          endif

          if (theta_pol.lt.0.d0) then
             theta_pol=theta_pol+2.d0*pi
          endif

          if (theta_pol.gt.2.d0*pi) then
             n_theta_pol=theta_pol/(2.d0*pi)
             theta_pol=theta_pol-2.d0*pi*n_theta_pol
          endif
           
          psi_loc=psi_rho(rho)

c          write(*,*)'prep3d.f in p_c_prof ,cefldy,cefldz',cefldy,cefldz
c          write(*,*)'prep3d.f in p_c_prof cnpar,cnper',cnpar,cnper
          write(*,*)'prep3d.f in p_c_prof psi_loc,theta_pol',
     &    psi_loc,theta_pol

          call CD_adj_LH_efficiency(cnpar,cnper,
     &    psi_loc,theta_pol,cefldy,cefldz,
     &    eff(is))
          write(*,*)'prep3d ADJ is,cnpar,eff(is)',is,cnpar,eff(is)


c         if (is.ge.65) call CD_adj_LH_efficiency_test(cnpar,cnper,
c     &      psi_loc,theta_pol,cefldy,cefldz,
c     &      eff(is))

      endif ! ieffic.eq.5


      if (ieffic.eq.6) then
c -------------------------------------------------------------------
c        calculation of current drive efficiency using adj function 
c        for n_harmonic=0 case (LH wave)
c--------------------------------------------------------------------
          rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2)

          if(rho_loc.eq.0.d0)then ! YuP [06,2015] To avoid 1/rho_loc=Inf
          cos_theta_pol=1.d0
          sin_theta_pol=0.d0
          else
          cos_theta_pol=(r_m-xma)/rho_loc
          sin_theta_pol=(z_m-yma)/rho_loc
          endif

          if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
          if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
          if (sin_theta_pol.ge.0.d0) then
             theta_pol=+dacos(cos_theta_pol)
          else  
             theta_pol=-dacos(cos_theta_pol)
          endif

          if (theta_pol.lt.0.d0) then
             theta_pol=theta_pol+2.d0*pi
          endif

          if (theta_pol.gt.2.d0*pi) then
             n_theta_pol=theta_pol/(2.d0*pi)
             theta_pol=theta_pol-2.d0*pi*n_theta_pol
          endif
           
          psi_loc=psi_rho(rho)

c          write(*,*)'prep3d.f in p_c_prof ,cefldx,cefldy,cefldz',
c     &                                     cefldx,cefldy,cefldz
c          write(*,*)'prep3d.f in p_c_prof cnpar,cnper',cnpar,cnper
          write(*,*)'prep3d.f in p_c_prof psi_loc,theta_pol',
     &    psi_loc,theta_pol

ctest
        u1=u_res(jwave,cnpar,temp,ye)
        write(*,*)'jwave,cnpar,temp,ye,u1',jwave,cnpar,temp,ye,u1
        effic_r=efficien(z_eff,u1,jwave,temp,den) ! here: just for comparison
        write(*,*)'asimptotic: z_eff,den,effic_r',z_eff,den,effic_r
cendtest
          call CD_adj_efficiency(cnpar,cnper,
     &    psi_loc,theta_pol,cefldx,cefldy,cefldz,
     &    eff(is))

c          if ((is.ge.262).and.(is.le.290)) then 
c            call CD_adj_efficiency_test(cnpar,cnper,
c    &       psi_loc,theta_pol,cefldx,cefldy,cefldz,
c    &       eff(is))
c          endif

          write(*,*)'prep3d ADJ is,cnpar,eff(is)',is,cnpar,eff(is)
      endif ! ieffic.eq.6
 110  continue

c      write(*,*)'p_c_prof after 110' 

c--------------------------------------------------------------------
c     delpower(erg/sec),delcurr(Ampere),r0(cm)
c-----------------------------------------------------------
cSmirnov970105 beg

c-----calculate parallel CD using the geometric factor 1/(2pi*r0) 
c      delcurr=delpow_e*0.5d0*(eff(is-1)+eff(is))/(2*pi*r0)
c      write(*,*)'1/(2*pi*r0)',1.d0/(2.d0*pi*r0)

c------geometric factor 1/(r*pi*R_mag_axis)
c      zfacgeom = 1.0d0 / (2.0d0 * pi * xma)/100.d0
c      write(*,*)'1/(2*pi*xma)/100.d0',1.d0/(2.d0*pi*xma)/100.d0

c-----calculate toroidal CD
      z0=0.5d0*(wz(is)+wz(is-1))
      r0m=r0*1.d-2
      z0m=z0*1.d-2      
      bmod=b(z0m,r0m,0.d0)
c      write(*,*)'z0m,r0m,bmod,rho',z0m,r0m,bmod,rho
      psiloc=psi_rho(rho)
c-----geometric factor ~ 1/b_averaged 
      zfacgeom = -2.d0 * pi * qsafety_psi(psiloc)/
     &     (dvol_dpsi(psiloc)*dpsimax*b_average(psiloc))/100.d0

c      write(*,*)'qsafety_psi(psiloc),b_average(psiloc)',
c     &           qsafety_psi(psiloc),b_average(psiloc)

c      write(*,*)'dvol_dpsi(psiloc),dpsimax',dvol_dpsi(psiloc),dpsimax
c      write(*,*)'zfacgeom',zfacgeom

c-----geometric factor ~ 1/b_ min/(1+aspct) old from Toray
c      zfacgeom = -2.d0 * pi * qsafety_psi(psiloc)/
c     &     (dvol_dpsi(psiloc)*dpsimax*bmin_psi(psiloc))/100.d0
c      aspct=(bmax_psi(psiloc)-bmin_psi(psiloc))/
c     &(bmax_psi(psiloc)+bmin_psi(psiloc))
c      zfacgeom = zfacgeom / (1.d0 + aspct)
c      write(*,*)'old zfacgeom',zfacgeom

      delcurr=delpow_e*0.5d0*(eff(is-1)+eff(is))*zfacgeom !toroidal current
                                                          !created by delpow_e
c      write(*,*)'prep3d zfacgeom', zfacgeom 

c      write(*,*)' prep3d p_c_prof delpow_e,eff(is-1),eff(is),delcurr',
c     +delpow_e,eff(is-1),eff(is),delcurr

cSAP080902 to plot delta power and delta current along the ray
      delpow_e_ar(is)=delpow_e
      delcur_par_ar(is)=delpow_e*0.5d0*(eff(is-1)+eff(is))
      
c-----toroidal and poloidal CD from old genray version      
      rho0=0.5*(spsi(is)+spsi(is-1)) ! the small radius
      rho0_pol=rho_lrho(rho0)
      poloidlen=rho0_pol*totlength*100.d0     ! poloidal length cm

c      write(*,*)'p_c_prof rho0_pol,totlength,,r0,poloidlen',
c     +rho0_pol,totlength,r0,poloidlen
c      write(*,*)'p_c_prof r0,2*pi*r0,eff(is)',r0,2*pi*r0
c      write(*,*)'p_c_prof eff(is)',
c     *eff(is)
c      write(*,*)'bphi,bpol,bmod',bphi,dsqrt(bz**2+br**2),bmod
c      write(*,*)'p_c_prof rho0,rho0_pol,poloidlen(cm)',
c     * rho0,rho0_pol,poloidlen
            
cSmirnov970105 end
      allcur=allcur+delcurr                     !total toroidal current
      
c      write(*,*)'allcur',allcur
c      write(*,*)'delcurr',delcurr

c      write(*,*)'p_c_prof rhoend,rhobegin',rhoend,rhobegin

      if(rhoend.gt.rhobegin) then
          eff_rho_max=eff(is)
          eff_rho_min=eff(is-1)
      else
          eff_rho_max=eff(is-1)
          eff_rho_min=eff(is)
      endif

c      write(*,*)'jbinmin,jbinmax',jbinmin,jbinmax     
       
      if (jbinmin.eq.jbinmax) then

c         write(*,*)'jbinmon=jbinmax,delpower,power(jbinmin)',
c     .              jbinmin,delpower,power(jbinmin)

         power(jbinmin)=power(jbinmin)+delpower

cSAP080731
c        write(*,*)'power(jbinmin)',power(jbinmin)
 
         cur_den_parallel(jbinmin)=cur_den_parallel(jbinmin)+
     &   delpow_e*0.5d0*(eff_rho_min+eff_rho_max)/binvol(jbinmin)

cSAP080731
c        write(*,*)'delpow_e,eff_rho_min,eff_rho_max,binvol(jbinmin)',
c     &             delpow_e,eff_rho_min,eff_rho_max,binvol(jbinmin)

c        write(*,*)'cur_den_parallel(jbinmin)',cur_den_parallel(jbinmin)

         current(jbinmin)=current(jbinmin)+delcurr   !toroidal current from bin

c         write(*,*)'cur_den_parallel(jbinmin)*binarea(jbinmin)',
c     &   cur_den_parallel(jbinmin)*binarea(jbinmin)

c         write(*,*)'current(jbinmin)',current(jbinmin)
c         write(*,*)'zfacgeom',zfacgeom
c         write(*,*)'0.5d0*(eff(is-1)+eff(is))',
c     &   0.5d0*(eff(is-1)+eff(is))
c         write(*,*)'is,eff(is-1),eff(is)',is,eff(is-1),eff(is)
c         write(*,*)'binarea(jbinmin)/binvol(jbinmin)',
c     &   binarea(jbinmin)/binvol(jbinmin)
c         write(*,*)'0.5d0*(eff_rho_min+eff_rho_max)',
c     &   0.5d0*(eff_rho_min+eff_rho_max)
c         write(*,*)'is,eff(is-1),eff(is)',is,eff(is-1),eff(is)
c         write(*,*)'eff_rho_min,eff_rho_max',eff_rho_min,eff_rho_max
      
         power_e(jbinmin)=power_e(jbinmin)+delpow_e
         power_i(jbinmin)=power_i(jbinmin)+delpow_i
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               power_s(jbinmin,kk)=power_s(jbinmin,kk)+delpow_s(kk)
            enddo
         endif
       
         power_cl(jbinmin)=power_cl(jbinmin)+delpow_cl     

         goto 30
      endif

c      write(*,*)'prep3d  after goto30 delpower rhomax,rhomin',
c     . delpower,rhomax,rhomin    
       
c      delrho=dmax1((rhomax-rhomin),hrho) 
cSm040426
      delrho=rhomax-rhomin

      ppow=delpower/delrho 

      pcur=delcurr/delrho
      
cSmirnov970106 beg
      ppow_e=delpow_e/delrho
      ppow_i=delpow_i/delrho
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            ppow_s(kk)=delpow_s(kk)/delrho
         enddo
      endif
      
      ppow_cl=delpow_cl/delrho
cSmirnov970106 end

c------------------------------------------------------------
c     power (erg/sec), current(A)
c------------------------------------------------------------

c      write(*,*)'jbinmin,power(jbinmin),ppow,(hrho*jbinmin-rhomin)',
c     &jbinmin,power(jbinmin),ppow,(hrho*jbinmin-rhomin)

      power(jbinmin)=power(jbinmin)+ppow*(hrho*jbinmin-rhomin)     

c      write(*,*)'power(jbinmin)',power(jbinmin)

cSmirnov970106 beg
      power_e(jbinmin)=power_e(jbinmin)+ppow_e*(hrho*jbinmin-rhomin)
      power_i(jbinmin)=power_i(jbinmin)+ppow_i*(hrho*jbinmin-rhomin)
      power_cl(jbinmin)=power_cl(jbinmin)+ppow_cl*(hrho*jbinmin-rhomin)
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            power_s(jbinmin,kk)=
     &           power_s(jbinmin,kk)+ppow_s(kk)*(hrho*jbinmin-rhomin)
         enddo
      endif

cSmirnov970106 end
c      write(*,*)'p_c_prof delpow_e',delpow_e

      cur_den_parallel(jbinmin)=cur_den_parallel(jbinmin)+
     &(delpow_e/delrho)*(hrho*jbinmin-rhomin)/binvol(jbinmin)*
     &0.5d0*(eff_rho_min+
     &      (eff_rho_min+(eff_rho_max-eff_rho_min)*
     &                   (hrho*jbinmin-rhomin)/(delrho)))

      current(jbinmin)=
     1 current(jbinmin)+pcur*(hrho*jbinmin-rhomin)  !toroidal current from bin

c     write(*,*)'p_c_prof jbinmin'
c     write(*,*)'current(jbinmin)',current(jbinmin)
c     write(*,*)'cur_den_parallel(jbinmin)*binarea(jbinmin)',
c    &cur_den_parallel(jbinmin)*binarea(jbinmin)

   
c      write(*,*)'jbinmax,power(jbinmax),ppow,(hrho*jbinmax-rhomin)',
c     &jbinmax,power(jbinmax),ppow,(hrho*jbinmax-rhomin)
 
      power(jbinmax)=power(jbinmax)+ppow*(rhomax-hrho*(jbinmax-1))

c     write(*,*)' power(jbinmax)', power(jbinmax)
cSmirnov970106 beg
      power_e(jbinmax)=power_e(jbinmax)+ppow_e*(rhomax-hrho*(jbinmax-1))
      power_i(jbinmax)=power_i(jbinmax)+ppow_i*(rhomax-hrho*(jbinmax-1))
      power_cl(jbinmax)=power_cl(jbinmax)+
     1                  ppow_cl*(rhomax-hrho*(jbinmax-1))
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            power_s(jbinmax,kk)=
     &          power_s(jbinmax,kk)+ppow_s(kk)*(rhomax-hrho*(jbinmax-1))
         enddo
      endif
cSmirnov970106 end

      cur_den_parallel(jbinmax)=cur_den_parallel(jbinmax)+
     &(delpow_e/delrho)*(rhomax-hrho*(jbinmax-1))/binvol(jbinmax)*
     &0.5d0*(eff_rho_max+
     &   (eff_rho_min+(eff_rho_max-eff_rho_min)*
     &   (rhomax-hrho*(jbinmax-1))/(delrho)))
     
      current(jbinmax)=
     1	  current(jbinmax)+pcur*(rhomax-hrho*(jbinmax-1)) !toroidal current from bin

c     write(*,*)'p_c_prof jbinmax'
c     write(*,*)'current(jbinmax)',current(jbinmax)
c      write(*,*)'cur_den_parallel(jbinmax)*binarea(jbinmax)',
c     &cur_den_parallel(jbinmax)*binarea(jbinmax)

c      write(*,*)'jbinmin,jbinmax',jbinmin,jbinmax
c     write(*,*)'power(jbinmin),power(jbinmax)',
c    &power(jbinmin),power(jbinmax)

c      write(*,*)'power_e(jbinmin)',
c     1           power_e(jbinmin)
c      write(*,*)'power_e(jbinmax)',
c     1           power_e(jbinmax)


      if(jbinmax.gt.(jbinmin+1)) then
         do j=(jbinmin+1),(jbinmax-1)
c            write(*,*)'j,power(j),ppow,hrho',j,power(j),ppow,hrho
            power(j)=power(j)+ppow*hrho
c            write(*,*)'power(j)',power(j)
cSmirnov970106 beg
            power_e(j)=power_e(j)+ppow_e*hrho
            power_i(j)=power_i(j)+ppow_i*hrho
            power_cl(j)=power_cl(j)+ppow_cl*hrho
            if (iabsorp.eq.3) then
               do kk=2,nbulk
                  power_s(j,kk)=power_s(j,kk)+ppow_s(kk)*hrho
               enddo
            endif

cSmirnov970106 end
            cur_den_parallel(j)=cur_den_parallel(j)+
     &      (delpow_e/delrho)*hrho/binvol(j)*
     &      (eff_rho_min+(eff_rho_max-eff_rho_min)*
     &      (hrho*(j-0.5d0)-rhomin)/delrho)
     
            current(j)=current(j)+pcur*hrho             !toroidal current from bins
         
c            write(*,*)'j,current(j),cur_den_parallel(j)*binarea(j)',
c     &      j,current(j),cur_den_parallel(j)*binarea(j)

c            write(*,*)'j,power(j)',j,power(j)

c            write(*,*)'j,power_e(j)',
c     .                 j,power_e(j)

         enddo
      endif

30    continue
      
      
c*******test
c       if (is.gt.1)then
c          write(*,*)'is,jbinmin,jbinmax.is',is,jbinmin,jbinmax
c          do ir=jbinmin,jbinmax
c            write(*,*)'ir,power(ir)',ir,power(ir)
c          enddo
c       endif
        
c      write(*,*)'allpower',allpower
c       write(*,*)'p_c_prof allcur', allcur
c      do ir=1,NR
c        write(*,*)'ir,power(ir),current(ir),temparr(ir),
c     1   spower(ir),scurrent(ir)',
c     2   ir,power(ir),current(ir),temparr(ir),
c     1   spower(ir),scurrent(ir)
c      enddo
c       do ir=1,NR
c        write(*,*)'ir,power(ir)',ir,power(ir)
c       enddo
c       spower_tot=0.d0
c       power_tot=0.d0
c       do ir=1,NR-1
c         spower_tot=spower_tot+spower(ir)
c         power_tot=power_tot+power(ir)
c       enddo
c       write(*,*)'delpwr(1)-delpwr(is)',delpwr(1)-delpwr(is)
c       write(*,*)'allpower,spower_tot,power_tot',
c     &            allpower,spower_tot,power_tot
c       write(*,*)'allpower,power_tot',
c     &            allpower,power_tot
c       write(*,*)'in p_c_prof before 99 return'
c       do ir=1,NR-1
c       write(*,*)'ir,cur_den_parallel(ir)',ir,cur_den_parallel(ir)
c       enddo
     
c99    return
      END


      subroutine write_lsc_output
c-------------------------------------------------------------------
c     Write data from LSC for TSC to LSC_for_TSC.dat
c     Radial profiles at uniform radial mesh i_rho=1,n_psi_TSC
c     rho_bin_center_lsc(i_rho)
c     power density:     power_dens_watt_m3_TSC_1D(i_rho) [watt/m**3]
c     current densityj:  j_rf_TSC_1D(i_rho)               [A/m**2]
c     derivatives:       d_j_rf_d_Edc_TSC_1D(i_rho))      [A/(V*m)]
c------------------------------------------------------------------
      implicit none
      include 'param.i'    

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
        
      include 'lsc_approach.i'

c-----locals
      integer i_unit,i_rho,kode

      if(myrank.ne.0) return
      
      i_unit=1
c      open(i_unit,file = 'LSC_for_TSC.dat',status='new',iostat=kode)
      open(i_unit,file = 'LSC_for_TSC.dat',iostat=kode)

      write(1,2)
      do i_rho=1,n_psi_TSC
         write(1,1)rho_bin_center_lsc(i_rho),
     &             power_dens_watt_m3_TSC_1D(i_rho),
     &             j_rf_TSC_1D(i_rho),
     &             d_j_rf_d_Edc_TSC_1D(i_rho)   
      enddo

 1    format(4(1pe15.8))
 2    format('rho pow dens [watt/m**3] cur dens [A/m**2] dJ_RF_d_Edc [A/(
     &V*m)]')
      close(1)

      return
      end
