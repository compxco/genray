c August 19, 2004
c This is a Weiss method relativistic, EC dispersion relation
c evaluator by Eric Nelson-Melby.
c Weiss method: Isaac Weiss, JCP 61, p. 403-416 (1985).

      BLOCK DATA eric_disp
      INTEGER PIN
      PARAMETER (PIN=20)  !SET PIN ALSO THE SAME IN Disp_Nelson_Melby
      DOUBLE PRECISION xlag(PIN), wlag(PIN)
      common /quadlag/ xlag,wlag
      DATA xlag /0.0705398896919887533666890D0, 
     + 0.372126818001611443794241D0, 
     + 0.91658210248327356466772D0, 1.707306531028343880688D0, 
     + 2.74919925530943212965D0, 
     + 4.0489253138508869224D0, 5.6151749708616165141D0,
     + 7.459017453671063310D0,  9.59439286958109677D0,
     + 12.03880254696431631D0, 14.8142934426307400D0, 
     + 17.9488955205193760D0, 21.4787882402850110D0,
     + 25.451702793186906D0, 29.932554631700612D0,
     + 35.0134342404790000D0, 40.8330570567285711D0,
     + 47.6199940473465021D0, 55.81079575006389889D0,
     + 66.524416525615753819D0/
      DATA wlag /0.16874680185111386214922D0,
     + 0.2912543620060682817168D0, 
     + 0.266686102867001288550D0, 0.166002453269506840031D0, 
     + 0.07482606466879237054D0, 
     + 0.02496441730928322107D0, 0.00620255084457223685D0,
     + 0.00114496238647690824D0, 
     + 0.00015574177302781197D0, 0.000015401440865224916D0, 
     + 1.0864863665179824D-6, 
     + 5.330120909556715D-8, 1.757981179050582D-9, 
     + 3.725502402512321D-11, 
     + 4.767529251578191D-13, 3.372844243362438D-15, 
     + 1.1550143395003988D-17, 1.5395221405823436D-20, 
     + 5.2864427255691578D-24, 1.656456612499023296D-28/
      END BLOCK DATA

      complex*16 function det(eps,nll,np)
c-----calculate determinant of the Maxwell system

      implicit none

c-----input
      complex*16 eps(3,3), !dielectric tensor
     &np                   !complex perpendicular refractive index
      real*8 nll           !parallel refractive index

c-----locals  
      complex*16 Kxx, Kxy, Kxz, Kyx, Kyy, Kyz, Kzx, Kzy, Kzz
      real*8 nlls
      complex*16 nps

      Kxx = eps(1,1)
      Kyy = eps(2,2)
      Kzz = eps(3,3)
      Kxy = eps(1,2)
      Kxz = eps(1,3)
      Kyx = eps(2,1)
      Kyz = eps(2,3)
      Kzx = eps(3,1)
      Kzy = eps(3,2)

c-----complex dispersion function
      nlls=nll**2
      nps=np**2

      det =(Kxx-nlls) * (Kyy-nlls-nps) * (Kzz-nps)
     .+  Kxy * Kyz * (Kzx+np*nll) 
     .+ (Kxz+np*nll) * Kyx * Kzy
     .- (Kzx+np*nll) * (Kyy-nlls-nps) * (Kxz+np*nll)
     .-  Kzy * Kyz * (Kxx-nlls)
     .- (Kzz - nps) * Kyx * Kxy

      return
      end   !check the formula for det .Here K is a dielectric tensor.
  
     
      subroutine eps_weiss_to_stix(eps_weiss,eps_stix)
c-----transform the dielectric tensor from Weiss to Stix coordinates
c     Weiss coordinates: wave vector K is in y-z plane 
c     Stix coordinates:  wave vector K is in x-z plane
c     Here z axis is along B field

      implicit none

c-----input
      complex*16 eps_weiss(3,3) !the dielectric tensor in Weiss system    

c-----output
      complex*16 eps_stix(3,3) !the dielectric tensor in Weiss system
 
      eps_stix(1,1)= eps_weiss(2,2) 
      eps_stix(2,2)= eps_weiss(1,1) 
      eps_stix(3,3)= eps_weiss(3,3)
      eps_stix(1,2)= eps_weiss(1,2) 
      eps_stix(2,1)=-eps_weiss(1,2) 
      eps_stix(1,3)=-eps_weiss(2,3) 
      eps_stix(3,1)=-eps_weiss(2,3) 
      eps_stix(2,3)= eps_weiss(1,3) 
      eps_stix(3,2)=-eps_weiss(1,3) 

      return
      end

C**** NOTE: to change resolution of integration grid, change
C**** paramters PIN and/or POUT both here and in iwfdispfuneps, below

      SUBROUTINE Disp_Nelson_Melby(T_e,N_parallel,X_e,Y_e,
     +  N_perp,eps,D)

      IMPLICIT NONE
      DOUBLE PRECISION T_e,N_parallel,X_e,Y_e
cBH141221      DOUBLE COMPLEX   N_perp,eps(1:3,1:3),D
      DOUBLE COMPLEX   N_perp,eps(3,3),D

      INTEGER NROOT,components,MAXNROOT,MAXCOMP,points,msg
      PARAMETER (MAXNROOT=15,MAXCOMP=7)
      INTEGER follow(MAXNROOT),lowtemp(MAXCOMP),omce,nopp
C**** MAXCOMP -- maximum number of Maxwellian components of different energy
C*** in the distribution function, COMP is the actual number to be used
      DOUBLE PRECISION alpha(MAXCOMP),f,n1,te(MAXCOMP),ompece(MAXCOMP)
      DOUBLE PRECISION maxwell(MAXCOMP),ompesqom(MAXCOMP),psi(MAXCOMP)
      DOUBLE PRECISION p0(MAXCOMP)
      COMPLEX*16 dispfun
      COMPLEX*16 x(MAXNROOT),nrootguess(MAXNROOT)
COLD NAG      DOUBLE PRECISION S18ACF, S18ADF
      DOUBLE PRECISION PI, MEC2
COLD NAG      EXTERNAL S18ACF, S18ADF
      INTEGER PIN,POUT,PTEMP,order,OUTCHUNK, POUTB
      PARAMETER (PIN=20, POUTB=8, POUT=256)  !POUT is most important for good accuracy in imaginary parts
C*** POUTB is the base Gauss-Legendre numbers (8 for now)
C***** n is really for 2*POUTB but +/- xi are both present with same weight
C*** POUT can now be any multiple of 2*POUTB
C**** NOTE: need to also declare these parameters in dispfun function
C**** number of points for inner integral (gauss-laguerre) PIN
C**** number of points for outer integral (gauss-legendre) POUT
      DOUBLE PRECISION xlag(PIN), wlag(PIN), xleg(POUT), wleg(POUT)
C** xlegb and wlegb are the base values of the Gauss-Legendre integration,
C** and then larger chunks can be made from transformations of these
      DOUBLE PRECISION xlegb(POUTB), wlegb(POUTB)
      DOUBLE PRECISION ALEG,BLEG,varmin,varmax,varstep,fmin,temin
      DOUBLE PRECISION A1,B1,C1,S1,T1
      DOUBLE PRECISION DBESK0, DBESK1

      common /plasma/ alpha,f,n1,psi,maxwell,ompesqom,p0,components,
     + lowtemp,nopp
      common /quadlag/ xlag,wlag
      common /quadleg/ xleg,wleg,OUTCHUNK
      common /debug/ msg
COLD NAG      EXTERNAL         D01BAX, D01BBF,D01BAZ
C*** D01BAX, D01BBF, D01BAZ are for the Gaussian quadrature
cBH141221      EXTERNAL dispfun
      EXTERNAL DBESK0, DBESK1

      integer IFAIL,varyc,k,iostat,j

      DATA PI /3.1415926535897932384626433832795D0/
      DATA MEC2 /510.99906D0/   !me*c^2 (in keV) electron rest mass
C*** Gauss-Laguerre quadrature points obtained from Mathematica (n=20)
C*** see BLOCK DATA above
C*** Gauss-Legendre quadrature points obtained from Mathematica and
C***  Abramowitz and Stegun for a few extra digits on the high weights
C***  n= 2*POUTB points, here only POUTB points because +/- xi are
C***  both present with the same weight
      DATA xlegb /0.09501250983763744018531934D0,
     + 0.281603550779258913230461D0, 0.45801677765722738634242D0,
     + 0.6178762444026437484467D0, 0.7554044083550030338951D0,
     + 0.865631202387831743880D0, 0.944575023073232576078D0,
     + 0.989400934991649932596D0/
      DATA wlegb /0.1894506104550684962853967D0,
     + 0.182603415044923588866764D0, 0.16915651939500253818931D0,
     + 0.1495959888165767320815D0,  0.124628971255533872052D0,
     + 0.09515851168249278481D0, 0.062253523938647892863D0,
     + 0.027152459411754094852D0/

C**** faster to use K0 and K1 to get K2 (necessary for Maxwell factor)
C*** and use K2=K0+2/z*K1
C*** use NAG routines K0 real argument S18ACF, K1 real argument S18ADF
C**** FOR SLATEC: use DBESK0 for double precision real argument K0,
C**** DBESK1 for double precision real argument K1

C**** See iwroots.inp for inputting variables
C*** (use iwroots_interactive to type in from keyboard)

C********* MAKE CHOICES OF WHICH TO VARY, VARY DIFFERENT PARTS
C***** THEN MAKE OUTPUT FILES
C components=1
C p0=0, omce=0 
C msg=0 

      components=1
C*** Need to set up te, ompesqom, p0, f, n1
      ompesqom(1)=X_e  ! omega_pe^2/omega^2
      te(1)=T_e   ! T_e in keV
      p0(1)=0.0D0
      n1=N_parallel
      if (abs(n1).ge.1.0) then
         print *,'WARNING: Magnitude of N_parallel > 1, results '
         print *,' are probably wrong. (The Weiss method with'
         print *,' transformation is designed for N_parallel<1)'
      endif
      f=1.0D0/Y_e    ! Y_e is omega_ce/omega, f is omega/omega_ce

      omce=0   ! for compatibility with older input files  (default omce=0)
      msg=0    ! msg=1 gives more debugging notes

      if (f .eq. 0) then
        print *,'ERROR: results not defined for f=omega/omega_ce=0'
        print *,' Change f to a non-zero value'
        stop
      endif

C****NOTE: in fortran, constants like this should have D to be counted as double
      Aleg=0.0D0
      Bleg=1.0D0
      IFAIL=0
COLD NAG      CALL D01BBF(D01BAX,Aleg,Bleg,0,PIN,wlag,xlag,IFAIL)
C*** NOW Gauss-Laguerre weights and abscissae are store for the standard
C*** 0 to infinity with weight exp(-x) in xlag (see DATA above)
C get Gauss-Laguerre weights and abscissae from 0 to infinity for
C  weight exp(-x). These can be adjusted for each pbar point in the
C  outer integration.
C      print *,'Gauss-Laguerre from ',Aleg,' to Infinity. Weight ',
C     + 'function exp(-',Bleg,'x) ',PIN,' points.'
C      print *,' xlag: ',xlag
C      print *,' wlag: ',wlag
C      read *

      if (msg.ge.1) then
        print *,'Gauss-Legendre (outer integral) with ',POUT,' points'
        print *,'Gauss-Laguerre (inner integral) with ',PIN,' points'
      endif
       DO j=1,components
        IFAIL=0
        alpha(j)=MEC2/te(j)
COLD NAG        maxwell(j)=S18ACF(alpha(j),IFAIL)+2.0D0*S18ADF(alpha(j),
COLD NAG     + IFAIL)/alpha(j)
C*** For SLATEC routines, argument should be less than 705.117
        if (alpha(j).gt.705.0D0) then
          maxwell(j)=0.0D0
        else
          maxwell(j)=DBESK0(alpha(j))+2.0D0*DBESK1(alpha(j))
     + /alpha(j)   ! use SLATEC Modified Bessel Functions of 2nd kind
        endif

C**** if low temperature (large alpha) then do a slightly different way
        if (maxwell(j).eq. 0.0D0 .or. te(j) .lt. 0.9) then
          lowtemp(j)=1
C*** use asymptotic expansion of K_2(alpha) for large alpha, and push
C*** the exp(alpha) into the integral and series
          maxwell(j)=dsqrt(PI/(2*alpha(j)))*(1.0D0+15.0D0/(8*alpha(j))+
     + 105.0D0/(128*alpha(j)*alpha(j))-315.0D0/(1024*alpha(j)*alpha(j)*
     + alpha(j)))
          if (msg.ge.1) print *,'Component ',j,' using low-',
     + 'temperature formulation. Exp(alpha)*Besselk(2,alpha)=',
     +  maxwell(j)
        else
          if (msg.ge.1) print *,'Component ',j,' Besselk(2,',
     + alpha(j),')=',maxwell(j)
          lowtemp(j)=0
        endif
        maxwell(j)=alpha(j)*alpha(j)/(8.0D0*maxwell(j))  !% normalization factor from Maxwellian
       ENDDO
CC%%% should be 3.822035003610377e13  for alpha=511/20

C************** Set up abscissae and weights for Gaussian quadrature
C---- First for outer integral, modification only depends on n1
C%% for Gauss-Legendre integration of outer integral, use POUT
C%% lower limit -n1, upper limit 1
C***   decide if the transformed integral from -n1 to 1 can be used or not
C***%% if alpha > 709*f/5 = 142 * f then use regular -1 to 1 way.
C*** or  f*te < 3.6
C*** find minimum f that will be used 
       fmin=f
C*** find minimum te that will be used 
       temin=te(1)

C*** if minimum temperature * f < 3.6 then use regular -1 to 1 integral.           
       if (abs(n1).gt.0.75D0 .or. temin*fmin.lt.3.6D0) then
         Aleg=-1.0D0
         nopp=1      !set no pprime flag
       else
         Aleg=-n1 
         nopp=0      !do not set no pprime flag
       endif
       Bleg=1.0D0
C%%
C%% try to find a Gauss method with more concentration in the center
C%% instead of around the edges.
C%%
C%%  get POUT point Gauss-Legendre abscissae and coefficients for Aleg to Bleg
COLD NAG: see file iwroots.f for getting Gauss-Legendre abscissae and coeff.
CCC from the NAG routines D01BBF

C**** make sure POUT is a multiple of 2*POUTB the fixed base value
       if (MOD(POUT,2*POUTB).NE.0) then
          print *,'ERROR: POUT must be a multiple of 2*POUTB'
          print *,'2*POUTB=',2*POUTB
          print *,' (For Gauss-Legendre integration) Change parameters'
          print *,' and compile again.'
          stop
       endif
C*** calculate how many chunks of 2*POUTB are needed, and construct the
C*** proper Gauss-Legendre abscissae and weights
       OUTCHUNK=POUT/(2*POUTB)
       C1=(Bleg-Aleg)/OUTCHUNK  ! break region into OUTCHUNK chunks
       A1=Aleg
       B1=Aleg+C1
       DO j=1,OUTCHUNK
         if (j.eq.OUTCHUNK) B1=Bleg  ! to avoid cumulative roundoff err
         S1=(B1-A1)/2.0D0
         T1=(B1+A1)/2.0D0
C*** order of points doesn't matter as long as proper weights accompany
C*** proper abscissae (addition is commutative)
         DO k=1,POUTB
C*** weights and abscissae for negative xi and positive xi
           wleg((j-1)*2*POUTB+k)=S1*wlegb(k)
           wleg((j-1)*2*POUTB+POUTB+k)=S1*wlegb(k)
           xleg((j-1)*2*POUTB+k)=-S1*xlegb(k)+T1
           xleg((j-1)*2*POUTB+POUTB+k)=S1*xlegb(k)+T1
         ENDDO
C*** now move to the next chunk
         A1=A1+C1
         B1=B1+C1
       ENDDO
C*** this kind of error will almost always happen, fixed by setting
C**** B1=Bleg at the last iteration of the loop
C*** if (A1.ne.Bleg) print *,'WARNING: A1.ne.Bleg, precision err'

C got Gauss-Legendre weights and abscissae for points from
C  -n1 to 1 (or possibly -1 to 1)
        if (msg.ge.1) then
          print *,'Gauss-Legendre from ',Aleg,' to ',Bleg,' ',POUT,
     + ' points.'
          print *,' xleg: ',XLEG
          print *,' wleg: ',WLEG
        endif

C  Need to have set up the Gaussian weights and abscissae first.
C
C Now try root-finding with muller.f
C
CC    print *,'1) omega/omega_ce 2) n|| 3) Te 4) omega_pe^2/omega^2

C**** function:      COMPLEX*16 FUNCTION DISPFUN(n2)  !!,PIN,POUT)

C*** Now output eps and D     
      D=DISPFUN(N_perp,eps)
      return

      END !subroutine Disp_Nelson_Melby

C------ dispersion function, root solver, bessel function series expansion

C--------- iwfdispfun.f (see iwfroots.f) dispersion function based on
C  fully relativistic results from Isaac Weiss, Journal of Computational
C  physics, 1985 ----------------------------------------
C (***** DISTINGUISHED FROM IWDISPFUN.F by having freely accessible routines
C  rather than NAG routines. hence iwfdispfun.f and iwfroots.f)
C 8 March 2003 -- added capability to have shifted Maxwellian
C  (from a Lorentz transformation, see Bornatici et al, PPCF 28, pp 629 (1985)
C
C  12 March 2003, fixed bug that crept in because
C  fortran doesn't check if a variable is declared or not, it just
C  stupidly assumes you want it to be zero! problem was that n1sqb was
C  sometimes used in the code as n1bsq (in p2 and p2fac), which was zero.
C  This only showed up when n1 was non-zero (and p0 nonzero) in the input file.
C   Now it should be fixed.
C
C  20 June 2003, found the problem area that was sometimes causing infinities
C  and NaN in the evaluation of the tensor elements, especially at low frequencies,
C  high temperatures or large n||. Set some limits on the combination of these
C  factors to use the regular -1 to 1 range in the outer integral instead of
C  trying to transform with -n1 to 1 and then using p2fac, pprime, etc.
C  this might make the evaluation longer, but will avoid many of the problems.
C
C  *** NOTE: testing has shown that numerical difficulties are encountered 
C   with low temperature (<about 1 or 2 keV) and/or  high parallel wave number
C  (larger than about 0.75 or less than -0.75). Also, near the fundamental
C  cyclotron frequency, the O-mode needs more resolution than POUT=128 to resolve
C  it well sometimes.
C
C *** NOTE: at exactly f=1.0, there is quite a bit of numerical difficulty.
C
C  15 April 2004 -- rewrote some of the routines to use free SLATEC routines
C   instead of NAG routines (for Bessel functions, and for Gaussian integration)
C  29 April 2004 -- when fail, set nperp to 1 (vacuum value) rather than to 0
C   (which just causes more problems)
C
C  17 May 2004 -- use bess_iw_complex instead of bess_rmg_complexfun for
C    cases where SLATEC routines fail.
C  18 May 2004 -- save f0(n) during series for the first dielectric element,
C   to re-use in other 5, to save time (uses more memory, however)
C
C
C NOTE: error returns from zbesj and zbesy
C ierr=0, normal.                    Computation completed.
C ierr=1, Input error.                  No computation.
C ierr=2, overflow.(Im z too large)     No computation.
C ierr=3, Precision warning          Computation completed.
C ierr=4, worse Precision error         No computation.
C ierr=5, Algorithmic error             No computation.

CCCC See if can be speeded up by saving bessel function evaluations
CCCC in infinite series parts for inte11, inte22, etc.
      COMPLEX*16 FUNCTION DISPFUN(n2,eps)  !!,PIN,POUT)
CCC**** n2 is complex, it is nperp to test
      INTEGER PIN,POUT,OUTCHUNK,IFAIL,components,MAXCOMP,MAXSERIES,msg
C*** use msg=1 for a moderate amount of messages, use msg=2 for 
C*** a lot of messages for real debugging, and msg=3 for same as 2 but
C*** a few less messages, especially in the inner integral, to go through faster
      PARAMETER (PIN=20, POUT=256, MAXCOMP=7,MAXSERIES=2000)
C*** PIN, right now fixed at 20. POUT can be as large as you like, a multiple
C**** of 16 (2*POUTB)
      INTEGER pbarind,bind,lowtemp(MAXCOMP),nopp
C*** if the flag nopp=1 then that means do not do the outer integral transformation
C*** that would save on Bessel function evaluations because it would cause other problems
      DOUBLE PRECISION alpha(MAXCOMP),f,n1,psi(MAXCOMP),p0(MAXCOMP)
      DOUBLE PRECISION maxwell(MAXCOMP),ompesqom(MAXCOMP)
C for a component in motion (shifted Maxwellian), use effective
C fb, n1b, n1sqb, n2b -- Doppler-shifted from f,n1,n1sq,n2
      DOUBLE PRECISION p,n1sq,p2fac,outweight,S2,S3,n1b,n1sqb,fb
      DOUBLE PRECISION f1n1p,f1n1p2,l,expfac,Ap,Bp,Bp2,xieval,T1
      DOUBLE PRECISION pperp,PI,X01AAF,C1,S1,psq,f1n1psq,f1n1p2sq
      DOUBLE PRECISION xlag(PIN), wlag(PIN), xleg(POUT), wleg(POUT)
      DOUBLE PRECISION m11,m33a,m33b,m13a,m13b,gcur,gtest,nmin,db
      DOUBLE PRECISION gam0,v0,onen1v0,n1shifta,n1shiftb,avoidint
      DOUBLE PRECISION SERIESTOL,MEC2
cBH141221      COMPLEX*16 eps(1:3,1:3)
      COMPLEX*16 eps(3,3)
      COMPLEX*16 n2, icomp, rho, Z1, n2b
      COMPLEX*16 f11a,f11b,f12a,f12b,f13a,f13b,f22a,f22b,f23a,f23b
      COMPLEX*16 f33a,f33b,inte11,inte12,inte13,inte22,inte23,inte33
      COMPLEX*16 fg11a,fg11b,fg12a,fg12b,fg13a,fg13b,fg22a,fg22b
      COMPLEX*16 fg23a,fg23b,fg33a,fg33b
      COMPLEX*16 i11,i12,i13,i22,i23,i33
      COMPLEX*16 sig11,sig12,sig13,sig22,sig23,sig33
      COMPLEX*16 fn11a(PIN),fn11b(PIN),fn12a(PIN),fn12b(PIN)
      COMPLEX*16 fn13a(PIN),fn13b(PIN),fn22a(PIN),fn22b(PIN)
      COMPLEX*16 fn23a(PIN),fn23b(PIN),fn33a(PIN),fn33b(PIN)
      INTEGER nsave(PIN),nsvc,icnt,doneseries,recount,rmg,ccnt,scaled
      COMPLEX*16 bessjsav(MAXSERIES)
      LOGICAL bjsavscaled(MAXSERIES)
      INTEGER f0snum,nsav2ind,nsav2(MAXSERIES),scaled2

C**** the ZS are the real and imaginary parts of nperp, and the FVEC are the
C**** real and imaginary parts of the dispersion function.
CC**** N should only =2
      COMPLEX*16 CY(3),CWRK(3),CJ(3)
      DOUBLE PRECISION CJR(3), CJI(3), CWRKR(3), CWRKI(3)
      INTEGER cjindex,icntstart
*     .. External Subroutines ..
COLD NAG:      EXTERNAL      S17DEF,S17DCF,X01AAF,bess_rmg_complex
      EXTERNAL zbesj, zbesy, bess_iw_complex
*     .. Executable Statements ..

      common /plasma/ alpha,f,n1,psi,maxwell,ompesqom,p0,components,
     + lowtemp,nopp
      common /quadlag/ xlag,wlag
      common /quadleg/ xleg,wleg,OUTCHUNK
      common /debug/ msg
C***** NOTE: xleg is decreasing, while xleg is increasing (as returned by NAG routines)
      DATA PI /3.14159265358979323846264338D0/
      DATA MEC2 /510.99906D0/   !me*c^2 (in keV) electron rest mass
      DATA db /1.D-4/
      DATA avoidint /1.D-9/   ! to avoid l=integer, add avoidint to l
      DATA SERIESTOL /1.D-12/  ! relative tolerance for series addition to integral
C*** usually SERIESTOL 1.d-6 is good enough

CC      n2=DCMPLX(ZS(1),ZS(2))
      n1sq=n1*n1
      icomp=DCMPLX(0.0D0,1.0D0)

c      print *,'Inside dispfun'
C      print *,'Gauss-Laguerre from ',Aleg,' to Infinity. Weight ',
C     + 'function exp(-',Bleg,'x) ',PIN,' points.'
C      print *,' xlag: ',xlag
C      print *,' wlag: ',wlag

C% alpha=ratio of mc^2 to electron temperature (here in keV)
C%
C% f= omega/omega_ce = frequency / cyclotron frequency
C% n1 = n_parallel
C% n2 = n_perpendicular
      if (msg.ge.1) PRINT *,'nperp (complex) ',n2

CC%% set up the functions for the various dielectric tensor elements

CC% b is omegastar, p is pbar 

C%%% Now do integration by Gauss-Laguerre techniques
C the abscissae and weights
C should already be loaded into global variables from
C the calling procedure.
CC for now, using already modified versions from the NAG call
CCC But, IF these are for standard values from 0 to infinity and exp(-x).
CCC for modified values A to infinity and exp(-B*x) for weight, then
CCC they can be modified as follows:
CCC weights should be multiplied by exp(-A*B)/B, abscissae x -> x/B+A

C%%
C%% Now do large loop to evaluate inner integral at pbar = all the
C%% xioutereval points.

      sig11=0.0D0  ! reset variables for sigma components (outside component loop, once only)
      sig12=0.0D0
      sig13=0.0D0
      sig22=0.0D0
      sig23=0.0D0
      sig33=0.0D0
C%disp('in iweissdisp');
C%keyboard
      DO 3 ccnt=1,components   ! go through (up to 7) (possibly shifted) Maxwellian components of the distribution function
      i11=0.0D0  ! reset variables for outer integral
      i12=0.0D0
      i13=0.0D0
      i22=0.0D0
      i23=0.0D0
      i33=0.0D0
      if (ompesqom(ccnt).gt.0.0D0) then  ! only bother if there is non-zero plasma density for this species
C*** Now if this component has relative motion (with respect to the background
C*** ions) use shifted n2, n1, etc.
      if (p0(ccnt).ne.0.0D0) then
         gam0=DSQRT(1.0D0+p0(ccnt)*p0(ccnt))
         v0=p0(ccnt)/gam0
         onen1v0=1.0D0-n1*v0
         n1b=(n1-v0)/onen1v0
C problem if n1b>1 !! (luckily, that will never happen for n1<1
C  because of the speed limit c on the electrons
CC**** need to transform xleg and wleg for shifted n1b
C*** on outer integral, but don't leave altered for other components
C see E.Nelson-Melby notebook#2, pg. 60 -- this has to do with 
C Gauss-Legendre weights and evaluation points for different limits
C (i.e. going from -n1b to 1 vs. -n1 to 1)
         n1shifta=(1.0D0+n1b)/(1.0D0+n1)
C if POUT.le.64 then 1 chunk, if between 64 and 128, 2 chunks, else 8 chunks
CCC old way for 1 chunk: n1shiftb=(1.0D0-n1b)/2.0D0-n1shifta*(1.0D0-n1)/2.0D0
C if POUT>64, there are two chunks of Gauss-Legendre points
CCC old way for 2 chunks: n1shiftb=(1.0D0-3*n1b)/4.0D0-n1shifta*(1.0D0-3*n1)/4.0D0
CCC           n1shiftc=(3.0D0-n1b)/4.0D0-n1shifta*(3.0D0-n1)/4.0D0
CCC n1shiftc = n1shiftb (try it!)
         n1shiftb=n1shifta*(n1-(1.0D0+n1)/(2*OUTCHUNK))+
     +  (1.0D0+n1b)/(2*OUTCHUNK)-n1b
CCC this is the new way for arbitrary number of chunks
         n1sqb=n1b*n1b
         n2b=n2/(gam0*onen1v0)
C**** relativistic doppler shifts change n1,n2,etc. if p0 is non zero
         if (n2b.ne.n2b) then
            print *,'Error: n2b.ne.n2b'
            print *,'n2=',n2,' p0=',p0(ccnt),' gam0=',gam0,
     +  'onen1v0=',onen1v0
            print *,'Setting n2=1'
cSm            n2=COMPLEX(1.0D0,0.0D0)
            n2=dcmplx(1.0D0,0.0D0)
            dispfun=dcmplx(0.0D0,0.0D0)
            return
         endif
         fb=f*gam0*onen1v0
      else
         n1shifta=1.0D0
         n1shiftb=0.0D0
         gam0=1.0D0
         v0=0.0D0
         onen1v0=1.0D0
         n1b=n1
         n1sqb=n1sq
         n2b=n2
         if (n2b.ne.n2b) then
            print *,'Error: n2b.ne.n2b'
            print *,'n2=',n2,' p0=',p0(ccnt)
            print *,'Setting n2=1'
cSm            n2=COMPLEX(1.0D0,0.0D0)
            n2=dcmplx(1.0D0,0.0D0)
            dispfun=dcmplx(0.0D0,0.0D0 )
            return
         endif
         fb=f
      end if
      if (msg.ge.2) then
        print *,'Component ',ccnt,' n1=',n1,' n1b=',n1b
        print *,'p0=',p0(ccnt),' v0=',v0,' gam0=',gam0
        print *,'onen1v0=',onen1v0,' n1shifta=',n1shifta
        print *,'n1shiftb=',n1shiftb
        print *,'n2=',n2,' n2b=',n2b
        print *,'f=',f,' fb=',fb
      endif
      DO 5 pbarind=1,POUT
        inte11=0.0D0  ! reset variables for inner integral
        inte12=0.0D0
        inte13=0.0D0
        inte22=0.0D0
        inte23=0.0D0
        inte33=0.0D0
        nsvc=0
        p=n1shifta*xleg(pbarind)+n1shiftb  ! Gauss-Legendre point
        C1=n1b*p  ! temp. holder
        if (nopp.eq.0) then
          p2=-(p*(n1sqb+1)+2*n1b)/(1+n1sqb+2*C1)  !%% this is pprime, for other side of -n1b in outer integral
          p2fac=(n1sqb-1)/(1+n1sqb+2*C1)    !%% factor accompanying G(pprime)
          p2fac=p2fac*p2fac
        else
          p2fac=0.0D0
        endif
C%%% -------- ALSO COMPARE WITH METHOD OF BREAKING UP INTO MANY
C%%%% INTEGRALS FROM N-1/2 to N+1/2, and add until relative error is
C%%%% small, like with infinite sum contribution.
C%%% then do rest in a Gauss-Laguerre way.
C%%%
C%%%%%%-----------------------------
        f1n1p=fb*(1+C1)
        f1n1psq=f1n1p*f1n1p
        if (nopp.eq.0) then
          f1n1p2=fb*(1+n1b*p2)
          f1n1p2sq=f1n1p2*f1n1p2
        endif
        psq=p*p
C** if low temperature, then push exp(alpha) through and redefine
C** variable of integration
        l=f1n1p/dsqrt(1-psq)
        if (l .eq. dnint(l)) then
           print *,'WARNING: l=integer, slightly adjusting to avoid ',
     + 'complicated evaluation.'
          l=l+avoidint
        end if
        nmin=dint(l)+1   ! int is basically like floor for positive numbers
        if (msg.ge.2) print *,'l=',l,' nmin=',nmin

C%
C% Compare to method of breaking integral up into chunks
C% around n-1/2 to n+1/2, and doing each with Gauss-Legendre
C%% techniques, until 

C%%% for our case: A=l, B=alpha/f*(1+n1b*p)  (In int from A to infinity of f(x) exp(-B*x)
        Ap=l
        Bp=alpha(ccnt)/f1n1p
        if (nopp.eq.0) then
          Bp2=alpha(ccnt)/f1n1p2
        else
          Bp2=0.0D0
        endif

C%-------------
        if (msg.ge.2) print *,'p=',p,' Ap=l=',Ap,' Bp=',Bp
        if (lowtemp(ccnt).eq.1) then
          expfac=exp(-Ap*Bp+alpha(ccnt))/Bp    ! part of weight in front of Gauss-Laguerre sum
          if (msg.ge.2) print *,'low temp, arg of exp:',
     +  -Ap*Bp+alpha(ccnt)
        else
          expfac=exp(-Ap*Bp)/Bp    ! part of weight in front of Gauss-Laguerre sum
          if (msg.ge.2) print *,'normal temp, arg of exp:',
     +  -Ap*Bp
        endif
        if (msg.ge.2) print *,'expfac=',expfac
        outweight=maxwell(ccnt)*wleg(pbarind)*n1shifta 
C outweight is the maxwellian and Gauss-Legendre weight
C (which is possibly adjusted due to Lorentz-transformed n1, if this
C  component is in relative motion to the background ions)
C******  expfac is part of the Gaussian-Laguerre weight  to go along
C*** with the standard part of the weight wlag(bind)
C***  outweight : maxwell is the alpha^2/K_2(alpha) maxwellian factor and
C**** wleg(pbarind) is the Guassian-Legendre weight (includes the modification
C**** for outer integral from -n1 to 1), only calculated once outside
C*** bind loop
C        print *,'p=',p,' expfac: ',expfac
C        print *,' p2fac: ',p2fac
        recount=0
        if (expfac.ne.0) then
          gcur=-1.0D0
          DO 20 bind=1,PIN
            xieval=xlag(bind)/Bp+Ap  !Gauss-Laguerre point
C%%% l < b < nmin+0.5, then g(b)=nmin
C%%% nmin+0.5 < b < infinity, g(b)=round(b)
            if (xieval .lt. nmin+0.5D0) then
              gtest=nmin
              if (msg.ge.2) print *,'gtest=nmin=',gtest
            else
              gtest=dnint(xieval)
              if (msg.ge.2) print *,'gtest=dnint(xieval)=',gtest
            end if
            if (gcur .ne. gtest) then
CCC need to calculate fo(g) -- same as below but for integer xieval
              gcur=gtest
              m11=(1.0D0-psq)*gcur*gcur/f1n1psq-1.0D0
              pperp=dsqrt(m11)
              if (msg.ge.2) print *,'(gcur.ne.gtest)m11(pperpsq)=',m11
     +   ,' pperp=',pperp
              m11=-m11  ! used as -pperp^2 in the matrix elements
CC%%% pperp is the same whether p is p or pprime
              rho=n2b*fb*pperp
              IFAIL=-1
C compute besselj functions  --- see bessjcomplex.f,
C calculate only J_nu and get J_-nu from the relation
C J_-nu = -1^nu J_nu  for integral nu, since gcur is always integral
CCC***** SPECIAL CASE IF gcur is less than 1, but this will never occur
COLD NAG              CALL S17DEF(-1+gcur,rho,3,'U',CJ,NZ,IFAIL)
C  Here is the bessel function call for the positive omegastars
C  the flag 1 is to use unscaled bessel functions, 3 is to return 3
C  in a row (-1+gcur, gcur, and 1+gcur)
cSm              CALL zbesj(realpart(rho),imagpart(rho),-1+gcur,1,3,

cSAP091104
c              write(*,*)'in FUNCTION DISPFUN(n2,eps) bef zbesj
c     &rho,dreal(rho),dimag(rho),-1+gcur',
c     &rho,dreal(rho),dimag(rho),-1+gcur

              CALL zbesj(dreal(rho),dimag(rho),-1+gcur,1,3,
     +          CJR, CJI, nz,IFAIL)
              if (IFAIL .ne. 0 .and. IFAIL .ne. 3) then
                print *,'WARNING: BesselJ failed for gcur=',gcur,
     + ' rho=',rho,' nperp=',n2b,
     + ' $$$ During integer part of inner integral $$$ ',
     + ' Stopping root finder for this root, and setting nperp to 1.'
                n2=dcmplx(1.0D0,0.0D0)
                dispfun=dcmplx(0.0D0,0.0D0)
                return
              endif
              DO I=1,3
                 CJ(I)=DCMPLX(CJR(I),CJI(I))
              END DO
C  Now CJ(1) = J_-1+g(rho), CJ(2) = J_g(rho), CJ(3)=J_1+g(rho)
CCC Now calculate J_-nu from J_-nu = -1 ^ nu J_nu
              if (dmod(-1+gcur,2.0D0) .eq. 0) then
                C1=1.0D0
              else
                C1=-1.0D0
              end if
              DO I=1,3
                CJ(I)=DCMPLX(CJR(I),CJI(I))
                CWRK(I)=CJ(I)*C1
                C1=-C1
              ENDDO
C  Now CWRK(1) = J_1-g(rho), CWRK(2)=J_-g(rho), CWRK(3)=J_-1-g(rho)
CCC /export/soft/nag/examples/source   
CCCC f77 bessjcomplex.f -lnag -o bessjcomplex.exe
C              print *,'rho(g): ',rho
C              print *,'nu1(g): ',-1+gcur,' J(3series) ',CJ
C              print *,'-nu1(g): ',1-gcur,' J(3series) ',CWRK
C**** Now bessel functions are calculated. Calculate integrand for G(r)
C**** work from the lowest level upward
CC              m11=-pperp*pperp  ! m12 is just -m11 (NOTE: m11 is assigned above)
              m33a=p*gcur/f1n1p    ! needs to be squared after use in next statement
              m13a=m33a*pperp   !m23 is the same as m13
              m33a=m33a*m33a
C              print *,'icomp: ',icomp,' m33a: ',m33a,
C     + ' pperp: ',pperp,' m13a=',m13a
              if (nopp.eq.0) then
                m33b=p2*gcur/f1n1p2
                m13b=m33b*pperp
                m33b=m33b*m33b
              endif
C*** now matrix elements are defined
C*** now put together with factor of gcur/f1n1p^2 and J_alpha1 * J_alpha2
C*** recall:
C  Now CJ(1) =  J_-1+b(rho), CJ(2) = J_b(rho),    CJ(3)=J_1+b(rho)
C  Now CWRK(1) = J_1-b(rho), CWRK(2)=J_-b(rho), CWRK(3)=J_-1-b(rho)
              C1=gcur/f1n1psq
              if (nopp.eq.0) then
                S1=gcur/f1n1p2sq
              else
                S1=0.0D0
              endif
  
              fg11a=m11*CWRK(3)*CJ(3)
              fg11b=S1*fg11a
              fg11a=C1*fg11a
  
              fg12a=-m11*CWRK(1)*CJ(3)
              fg12b=S1*fg12a
              fg12a=C1*fg12a
  
              fg13a=CWRK(2)*CJ(3)
              if (msg.ge.2 .and. msg.ne.3) then
                print *,'Inner point ',bind,' of ',PIN
                print *,'CWRK(2)*CJ(3) ',fg13a,' m13a=',m13a
              endif
              fg13b=-icomp*m13b*S1*fg13a
              fg13a=-icomp*m13a*C1*fg13a
  
              fg22a=m11*CWRK(1)*CJ(1)
              fg22b=S1*fg22a
              fg22a=C1*fg22a
  
              fg23a=CWRK(1)*CJ(2)
              fg23b=-icomp*m13b*S1*fg23a
              fg23a=-icomp*m13a*C1*fg23a

              fg33a=CWRK(2)*CJ(2)
              fg33b=m33b*S1*fg33a
              fg33a=m33a*C1*fg33a
              if (msg.ge.2 .and. msg.ne.3) 
     +  print *,'fg11a=',fg11a,' fg13a=',fg13a
C******** Now save these for later use in the extra integral and the
C******** infinite series
              nsvc=nsvc+1
              nsave(nsvc)=gcur
              fn11a(nsvc)=fg11a
              fn11b(nsvc)=fg11b
              fn12a(nsvc)=fg12a
              fn12b(nsvc)=fg12b
              fn13a(nsvc)=fg13a
              fn13b(nsvc)=fg13b
              fn22a(nsvc)=fg22a
              fn22b(nsvc)=fg22b
              fn23a(nsvc)=fg23a
              fn23b(nsvc)=fg23b
              fn33a(nsvc)=fg33a
              fn33b(nsvc)=fg33b
              
            end if   ! of calculating new f0(g) (if gcur.ne.gtest)
CC   xieval=xi(bother)/Bp+Ap;
CCC Now I need to evaluate the 6 bessel functions needed for this pbar and omegastar
CCC (i.e. this p,p2 and xieval). if xieval is integral, need to evaluate everything twice
CCC and take average, to avoid 0/0.
C
C  J bessel functions needed at rho=n2b*fb*pperp, and order
C -1-b, -b, 1-b   where b is xieval (always positive)
C -1+b, b, 1+b

C  alpha1=[-1,1,0;1,1,1;0,1,0];  %% -omegastar
C  alpha2=[1,1,1;1,-1,0;1,0,0];  %% +omegastar
C            print *,'-------------------f1n1p*f1n1p: ',f1n1p*f1n1p
 23         m11=(1-psq)*xieval*xieval/f1n1psq-1
            pperp=dsqrt(m11)
            m11=-m11  ! used as -pperp^2 in the matrix elements
CC%%% pperp is the same whether p is p or pprime
            rho=n2b*fb*pperp
            if (rho.ne.rho) then
               print *,'Error: rho.ne.rho'
               print *,'n2b=',n2b,' fb=',fb,' pperp=',pperp
               print *,'Setting n2=1'
cSm               n2=COMPLEX(1.0D0,0.0D0)
               n2=dcmplx(1.0D0,0.0D0)
               dispfun=0.0D0
               return
            endif
C  [bj1,err1]=besselj(alpha1(i,j)-b,rho);
C  [bj2,err2]=besselj(alpha2(i,j)+b,rho);
            if (xieval .eq. dnint(xieval)) then
              print *,'IIIIII integral xieval: ',xieval
              xieval=xieval-db
              if (xieval .lt. 0) then
                print *,'ERROR: xieval (omegastar) < 0. xieval=',xieval
C** Note: this should not happen with Gauss-Legendre formula, because xieval will be 0
C* only when xieval=l=0 (i.e. at integration limit).
                stop
              end if
              recount=2
              go to 23
CC%%% if integral, take average of values just around the integral value
CC%%% to avoid evaluating 0/0 without doing messy derivative of bessel
CC%%% function with respect to both order and argument
            else
CCC              PI=X01AAF()  ! NAG pi storage
              IFAIL=0
              rmg=0
C compute besselj functions  --- see bessjcomplex.f,
C calculate both J_nu and Y_nu functions and get J_-nu from the combination
CCC***** SPECIAL CASE IF xieval is less than 1 !!!!
              if (rho .eq. 0.0D0) then
C*** if xieval=1, then need J_0(0), J_1(0), J_2(0)   =  1, 0, 0
C****                       J_0(0), J_-1(0), J_-2(0)  = 1, 0, 0
C*** else any J_nu(0) =0 if nu is not 0
                if (xieval .eq. 1.0D0) then
                  CWRK(1)=1.0D0
                  CJ(1)=1.0D0
                else
                  CWRK(1)=0.0D0
                  CJ(1)=0.0D0
                end if
                CWRK(2)=0.0D0
                CWRK(3)=0.0D0
                CJ(2)=0.0D0
                CJ(3)=0.0D0
              else
                if (-1+xieval .lt. 0 .and. msg .ge. 1) then
                  print *,' 0 < xieval < 1, calculate special'
CC***** NAG S17DCF is for Bessel Y function, S17DEF is for Bessel J function
COLD NAG                  CALL S17DCF(xieval,rho,2,'U',CY,NZ,CWRK,IFAIL)
cSm                  CALL zbesy(realpart(rho),imagpart(rho),xieval,1,2,
                  CALL zbesy(dreal(rho),dimag(rho),xieval,1,2,
     +             cjr,cji,NZ,cwrkr,cwrki,IFAIL)
            if (IFAIL .ne. 0 .and. IFAIL .ne. 3) then
              print *,'WARNING: BesselY failed for nu(xieval)=',
     + xieval,' rho=',rho,' nperp=',n2b,
     + ' %%% During inner integral, -1+xieval<0 %%%',
     + ' Stopping root finder for this root, and setting nperp to 1.'
              n2=dcmplx(1.0D0,0.0D0)
              dispfun=dcmplx(0.0D0,0.0D0)
              return
            endif
                  DO I=1,2
                     CY(I)=DCMPLX(cjr(I),cji(I))
                  ENDDO
C  Now CY(1) = Y_b(rho), CY(2)=Y_1+b(rho)
                  CY(3)=CY(2)
                  CY(2)=CY(1)
                  IFAIL=0
C                  CALL S17DCF(1-xieval,rho,1,'U',CY,NZ,CWRK,IFAIL)
cSm                  CALL zbesy(realpart(rho),imagpart(rho),1-xieval,1,1,
                  CALL zbesy(dreal(rho),dimag(rho),1-xieval,1,1,
     +             cjr,cji,NZ,cwrkr,cwrki,IFAIL)
            if (IFAIL .ne. 0 .and. IFAIL .ne. 3) then
              print *,'WARNING: BesselY failed for nu(1-xieval)=',
     + 1-xieval,' rho=',rho,' nperp=',n2b,
     + ' ^^^ During inner integral, -1+xieval<0 ^^^',
     + ' Stopping root finder for this root, and setting nperp to 1.'
              n2=dcmplx(1.0D0,0.0D0)
              dispfun=dcmplx(0.0D0,0.0D0)
              return
            endif
                  CY(1)=DCMPLX(cjr(1),cji(1))
C  Now CY(1) = Y_1-b(rho), CY(2) = Y_b(rho), CY(3)=Y_1+b(rho)
                  IFAIL=0
CC***** NAG S17DEF is for Bessel J function
COLD NAG                  CALL S17DEF(xieval,rho,2,'U',CJ,NZ,IFAIL)
cSm                  CALL zbesj(realpart(rho),imagpart(rho),xieval,1,2,
                  CALL zbesj(dreal(rho),dimag(rho),xieval,1,2,
     +             cjr,cji,NZ,IFAIL)
                  if (IFAIL .ne. 0 .and. IFAIL .ne. 3) then
                print *,'WARNING: BesselJ failed for nu=',xieval,
     + ' rho=',rho,' nperp=',n2b,
     + ' *** During inner integral, -1+xieval<0 ***',
     + ' Stopping root finder for this root, and setting nperp to 1.'
                    n2=dcmplx(1.0D0,0.0D0)
                    dispfun=dcmplx(0.0D0,0.0D0)
                    return
                  endif
                  DO I=1,2
                     CJ(I)=DCMPLX(cjr(I),cji(I))
                  ENDDO
C  Now CJ(1) = J_b(rho), CJ(2)=J_1+b(rho)
                  CJ(3)=CJ(2)
                  CJ(2)=CJ(1)
                  IFAIL=0
COLD NAG                  CALL S17DEF(1-xieval,rho,1,'U',CJ,NZ,IFAIL)
cSm                  CALL zbesj(realpart(rho),imagpart(rho),1-xieval,1,1,
                  CALL zbesj(dreal(rho),dimag(rho),1-xieval,1,1,
     +             cjr,cji,NZ,IFAIL)
                  if (IFAIL .ne. 0 .and. IFAIL .ne. 3) then
                print *,'WARNING: BesselJ failed for nu=',1-xieval,
     + ' rho=',rho,' nperp=',n2b,
     + ' @@@ During inner integral, -1+xieval<0 @@@',
     + ' Stopping root finder for this root, and setting nperp to 1.'
                    n2=dcmplx(1.0D0,0.0D0)
                    dispfun=dcmplx(0.0D0,0.0D0)
                    return
                  endif
                  CJ(1)=DCMPLX(cjr(1),cji(1))
C  Now CJ(1) = J_1-b(rho), CJ(2) = J_b(rho), CJ(3)=J_1+b(rho)
CCC Now calculate J_-nu from J_nu and Y_nu, put into variables CWRK
                  C1=DCOS(xieval*PI)
                  S1=DSIN(xieval*PI)
CC      WRITE (NOUT,*) 'PI: ',PI,' FNU: ',FNU,' Cos: ',C1,' Sin: ',S1
                  DO I=2,3
                    CWRK(I)=CJ(I)*C1-CY(I)*S1
                    C1=-C1
                    S1=-S1
                  ENDDO
C after DO loop, C1 and S1 are back to Cos(xieval*PI) and Sin(xieval*PI)
C                C1=DCOS((1-xieval)*PI)= -cos(xieval*PI) so just use -C1 in relation below
C                S1=DSIN((1-xieval)*PI)= sin(xieval*PI) so just use S1 in relation below
CC      WRITE (NOUT,*) 'PI: ',PI,' FNU: ',FNU,' Cos: ',C1,' Sin: ',S1
                  Z1=-CJ(1)*C1-CY(1)*S1  !CWRK(1)=J_-1+b
CCC Now swap the first one in CWRK and CJ to be as expected below
                  CWRK(1)=CJ(1)
                  CJ(1)=Z1
C  Now CJ(1) =  J_-1+b(rho), CJ(2) = J_b(rho),    CJ(3)=J_1+b(rho)
C  Now CWRK(1) = J_1-b(rho), CWRK(2)=J_-b(rho), CWRK(3)=J_-1-b(rho)


                else  ! xieval >= 1 (the usual case)

C***** S17DCF is Bessel Y function
                  IFAIL=1
COLD NAG                  CALL S17DCF(-1+xieval,rho,3,'U',CY,NZ,CWRK,IFAIL)
C   zbesy is SLATEC Bessel Y function -- first 1 is for unscaled output,
C   second number, 3, is for number of Bessel functions to return in series
cSm                  CALL zbesy(realpart(rho),imagpart(rho),-1+xieval,1,3,
                   CALL zbesy(dreal(rho),dimag(rho),-1+xieval,1,3,
     +             cjr,cji,NZ,cwrkr,cwrki,IFAIL)
C**** NOTE: if fails, handle just below using bess_rmg_complex function
                  DO I=1,3
                     CY(I)=DCMPLX(cjr(I),cji(I))
                  ENDDO

                  if (IFAIL .ne. 0) then
C***** NOTE: failure on NAG and SLATEC routines here during inner integral
C****  is almost always due to large nu but small rho -- perfect
C**** for series expansion. So use my own bessel function of the reduced series (with
C*** the (z/2)^nu dependence taken out -- it will cancel in the products anyway.
C** CHECK THAT |z|/rho is small enough to allow good calculation of
C the series expansion for both +/- xieval. The +nu part is more restrictive
C so check only positive nu, with the minimum nu. nu=-1+xieval, z=rho
C just use boundary for nu>3 (NAG probably wouldn't fail with nu<3)
C    pslim=0.15911+4.1*sqrt(nu-3);  /* fit minus 2 for extra accuracy */
CC now if zmag < pslim then OK to do positive series limit
C                    if (CDABS(rho) .ge. 
C     + 0.15911D0+4.1D0*DSQRT(-4.0D0+xieval)) then
C                      print *,'ERROR: SLATEC routine for Bessel ',
C     + 'function failed, but rho/nu is too large for series limit.'
C                      print *,'nu=',-1+xieval,' rho=',rho,' nperp=',n2b
CC                      stop
C                    endif
C***** With bess_iw_complex.f, calculate all 6 dielectric tensor element
C**** bessel products.
                    if (msg.ge.1) then
COLD NAG                      print *,'bessel Y (S17DCF) failed IFAIL=',IFAIL
                      print *,'bessel Y (zbesy) failed IFAIL=',IFAIL
                      print *,'nu=',xieval,' rho=',rho
                      print *,'Calculating with bess_iw'
                    endif
C**** Note: xieval will not be an integer, so using the hypergeometric series expansion should not encounter problems
C      call bess_iw_complex(nu,z,A11,A12,A13,A22,A23,A33)
                    call bess_iw_complex(xieval,rho,CJ(1),CJ(2),
     +  CJ(3),CWRK(1),CWRK(2),CWRK(3))

C*** this returns the products already needed for the dielectric tensor elements.
                    rmg=1
                    if (msg.ge.1) then
                       print *,'A11,A12,A13 ',CJ
                       print *,'A22,A23,A33 ',CWRK
                    endif
                  else  ! regular SLATEC routine zbesy worked
C  Now CY(1) = Y_-1+b(rho), CY(2) = Y_b(rho), CY(3)=Y_1+b(rho)
                    IFAIL=1
COLD NAG                    CALL S17DEF(-1+xieval,rho,3,'U',CJ,NZ,IFAIL)
cSm                    CALL zbesj(realpart(rho),imagpart(rho),-1+xieval,
                    CALL zbesj(dreal(rho),dimag(rho),-1+xieval,
     +                   1,3,CJR, CJI, nz,IFAIL)
                    if (IFAIL .ne. 0 .and. IFAIL .ne. 3) then
                       print *,'bessel J (zbesj) failed IFAIL=',IFAIL
                       print *,'(Yet zbesy was OK).'
                       print *,'nu_start=',-1+xieval,' rho=',rho
                       stop
                    end if
                    DO I=1,3
                      CJ(I)=DCMPLX(CJR(I),CJI(I))
                    END DO
C  Now CJ(1) = J_-1+b(rho), CJ(2) = J_b(rho), CJ(3)=J_1+b(rho)
CCC Now calculate J_-nu from J_nu and Y_nu, put into variables CWRK
                    C1=DCOS((-1+xieval)*PI)
                    S1=DSIN((-1+xieval)*PI)
CC      WRITE (NOUT,*) 'PI: ',PI,' FNU: ',FNU,' Cos: ',C1,' Sin: ',S1
                    DO I=1,3
                      CWRK(I)=CJ(I)*C1-CY(I)*S1
                      C1=-C1
                      S1=-S1
                    ENDDO
                  end if
                end if
              end if
C  Now CWRK(1) = J_1-b(rho), CWRK(2)=J_-b(rho), CWRK(3)=J_-1-b(rho)
CCC /export/soft/nag/examples/source   
CCCC f77 bessjcomplex.f -lnag -o bessjcomplex.exe
C              print *,'nu1: ',-1+xieval,' J(3series) ',CJ
C              print *,'-nu1: ',1-xieval,' J(3series) ',CWRK
            end if
C**** Now bessel functions are calculated. Calculate integrand for G(r)
C**** work from the lowest level upward
CC            m11=-pperp*pperp  ! m12 is just -m11 (NOTE: m11 is assigned above)
            m33a=p*xieval/f1n1p    ! needs to be squared after use in next statement
            m13a=m33a*pperp   !m23 is the same as m13
            m33a=m33a*m33a
            if (nopp.eq.0) then
              m33b=p2*xieval/f1n1p2
              m13b=m33b*pperp
              m33b=m33b*m33b
            endif
C*** now matrix elements are defined
C*** now put together with factor of xieval/f1n1p^2 and J_alpha1 * J_alpha2
C*** recall:
C  Now CJ(1) =  J_-1+b(rho), CJ(2) = J_b(rho),    CJ(3)=J_1+b(rho)
C  Now CWRK(1) = J_1-b(rho), CWRK(2)=J_-b(rho), CWRK(3)=J_-1-b(rho)
C  unless bess_iw_complex has been called, in which case the products
C  are already done, see notebook #3 for more details.
            C1=xieval/f1n1psq
            if (nopp.eq.0) then
              S1=xieval/f1n1p2sq
            else
              S1=0.0D0
            endif
            T1=PI*xieval
C pperp is available as a DOUBLE variable
C Z1 is still 4/rho^2
            if (msg.ge.2 .and. msg.ne.3) then
               print *,'m11: ',m11,' CWRK(3): ',CWRK(3)
               print *,'CJ(3): ',CJ(3)
            endif
C**** the A11, A12, etc. correspond to CJ(1), CJ(2), CJ(3), etc.
C      call bess_iw_complex(nu,z,A11,A12,A13,A22,A23,A33)
C                    call bess_iw_complex(xieval,rho,CJ(1),CJ(2),CJ(3),
C     +  CWRK(1),CWRK(2),CWRK(3))
            if (rmg.eq.1) then
               f11a=m11*CJ(1)
            else
               f11a=m11*CWRK(3)*CJ(3) ! type 1 in notebook 2
            endif
            if (msg.ge.2) print *,'f11a(1): ',f11a
            f11b=S1*f11a
            f11a=C1*f11a
            if (msg.ge.2) print *,'f11a(2): ',f11a

            if (rmg.eq.1) then
              f12a=-m11*CJ(2)
            else
              f12a=-m11*CWRK(1)*CJ(3)  ! type 2 in notebook 2
            end if
            f12b=S1*f12a
            f12a=C1*f12a

            if (rmg.eq.1) then
              f13a=CJ(3)
            else
              f13a=CWRK(2)*CJ(3)  ! type 3 in notebook 2
            end if
            f13b=-icomp*m13b*S1*f13a
            f13a=-icomp*m13a*C1*f13a

            if (rmg.eq.1) then
              f22a=m11*CWRK(1)
            else
              f22a=m11*CWRK(1)*CJ(1)  ! type 4 in notebook 2
            end if
            f22b=S1*f22a
            f22a=C1*f22a

            if (rmg.eq.1) then
              f23a=CWRK(2)
            else
              f23a=CWRK(1)*CJ(2)   ! type 5 in notebook 2
            end if
            f23b=-icomp*m13b*S1*f23a
            f23a=-icomp*m13a*C1*f23a

            if (rmg.eq.1) then
              f33a=CWRK(3)
            else
              f33a=CWRK(2)*CJ(2)   ! type 6 in notebook 2
            end if
            f33b=m33b*S1*f33a
            f33a=m33a*C1*f33a
C also keep track of g and fgs
C****** now calculate integrands
            C1=dcos(T1)    ! T1 is PI*xieval
            C1=C1*C1
            S1=wlag(bind)/sin(T1)   ! include weight for Gauss-Laguerre integration
            T1=exp(Bp*(xieval-gcur))*C1
C*** WATCH OUT: here for extremely low temperatures this could become infinity !

C**** multiply by weights (maxwellian, expfac, wlag, etc.) outside of loop
C*** Actually, 14 Mar 2003, changed to multiplying here inside loop
            S1=S1*expfac*outweight
            if (recount .gt. 0) then
              S1=S1/2.0D0
              recount=recount-1
            end if
            inte11=inte11+(f11a-fg11a*T1)*S1
            inte12=inte12+(f12a-fg12a*T1)*S1
            inte13=inte13+(f13a-fg13a*T1)*S1
            inte22=inte22+(f22a-fg22a*T1)*S1
            inte23=inte23+(f23a-fg23a*T1)*S1
            inte33=inte33+(f33a-fg33a*T1)*S1
            if (msg.ge.2 .and. msg.ne.3) then
              print *,'Inner point ',bind,' of ',PIN
              print *,'f11a=',f11a,' fg11a=',fg11a
              print *,'T1=(exp(Bp*(xieval-gcur)))*cos^2 ',T1,
     + ' C1=',C1
              print *,'inte11=',inte11
            endif
CC%%% Then with pprime
            if (p2fac .ne. 0) then
              T1=exp(Bp2*(xieval-gcur))*C1
C**** S1 is already the exponential factor and the outer Gauss-Laguerre weight
              S1=S1*p2fac*exp(xieval*(Bp-Bp2))
C**** NOTE: could be a problem when T1*S1 should be finite
C**** but alone the xieval-gcur part overflows or something.
              inte11=inte11+(f11b-fg11b*T1)*S1
              inte12=inte12+(f12b-fg12b*T1)*S1
              inte13=inte13+(f13b-fg13b*T1)*S1
              inte22=inte22+(f22b-fg22b*T1)*S1
              inte23=inte23+(f23b-fg23b*T1)*S1
              inte33=inte33+(f33b-fg33b*T1)*S1
              if (msg.ge.2 .and. msg.ne.3) then
                print *,'Bp2 ',Bp2,' xieval ',xieval,' gcur ',gcur
                print *,'C1 ',C1,' p2fac ',p2fac,' Bp ',Bp
                print *,'arg of exp in S1: xieval*(Bp-Bp2) ',
     +  xieval*(Bp-Bp2)
                print *,'T1=(exp(Bp2*(xieval-gcur)))*cos^2 ',T1,
     + ' S1new=',S1
                print *,'inte11=',inte11
              endif
            end if              
C**** HERE: problem with going too far in inner integral with low temperatures,
C*** 64 points are too many, do something about becoming NaN
C*** THERE CAN BE SEVERAL PROBLEMS WITH THIS SECTION WHEN ALPHA IS TOO LARGE
C*** (i.e. too low temperatures) arising from the large arguments of the exponentials
C** Easy solution: (may take a few more points to evaluate) -- don't use the
C*** p2fac terms, just go from -1 to 1 in the outer integrand, not from -n1 to 1
C*** with the transformation of the outer integral.
C****** There could also be similar problems with this section if n1 is too large
C***** (too close to +/- 1) -- so also just do the -1 to 1 without the outer integral
C**** transformation for large n1.
C***
C*** Some research into the possible range of xieval*(Bp-Bp2) shows that
C*** the factor involving pbar and pbarprime and n1 can be limited to about
C*** no larger than 5 (see iwdispresearch.m for more details)
C%%%%%%% ---------- End of the function to return integrand

C            print *,'p: ',p,' xieval: ',xieval,' pperp: ',pperp
C            print *,'rho: ',rho
            if (recount .gt. 0) then
              xieval=xieval+2*db
              go to 23
            end if
 20       continue   ! inner loop, bind
CCC**** Multiply by expfac and outweight before adding to sum of outer loop
CC**** from testing, it seems that perhaps it would be better to multiply
CC*** things inside the inner loop by C1, to avoid overflows
CCC -- although this increases the number of flops, it prevents overflows
CCC and could actually increase the accuracy (less round-off perhaps)
CC so move the multiplication by expfac*outweight inside the loop.
          i11=i11+inte11
          i12=i12+inte12
          i13=i13+inte13
          i22=i22+inte22
          i23=i23+inte23
          i33=i33+inte33
CC*** 14 Mar 2003 -- put factor of expfac*outweight inside integral.
          if (msg.ge.2) then
            print *,'Outer loop point p: ',p
            print *,'i11=',i11,' i12=',i12,' i13=',i13
            print *,'i22=',i22,' i23=',i23,' i33=',i33
C            print *,'C1 (expfac*outweight): ',C1
            print *,'inte11=',inte11,' inte12=',inte12,' inte13=',inte13
            print *,'inte22=',inte22,' inte23=',inte23,' inte33=',inte33
            if (i11 .ne. i11) stop
          endif
CCC Now add other parts -- the initial finite integral and the infinite series that comes
CCC from the poles around integral omegastar (xieval)
C%%%%%------------------------------------------
C%%% Now do the integral from a to b
C%% Same integral for all components, because f(nmin) is a constant
C%% see above
          if (dmod(nmin+1,2.0D0) .eq. 0) then
            S1=1.0D0
          else
            S1=-1.0D0
          end if
C**** now S1 is -1^(nmin+1), so -S1 is -1^nmin
          if (l .lt. nmin-0.5) then
C %%% limits for leftover principal value integral
C %%% at beginning of integration range
            C1=-(dlog(dtan(PI*l/2)*S1)+dcos(PI*l))/PI
          else
            C1=-(dlog(-dtan(PI*(nmin-l/2))*S1)+
     +  dcos(PI*(2*nmin-l)))/PI
C%%% exact integral of cos^2(pi*x)/sin(pi*x) from lima to limb
          end if
Cinteab=inte*maxwell;  %% fijreg1_weiss has maxwell factor in it
C***** now first saved in nsave should be nmin
C**** unless Gauss-Laguerre points are relatively spread out, in which
C**** case I'll have to recalculate the ones that are missing.
          if (nsvc .lt. 1 .or. nsave(1) .gt. nmin) then
C*** this will rarely happen, unless PIN is not so many points,
C**** in which case it might happen -- the first saved integer
C**** evaluation of the integrand is beyond nmin. In this case,
C**** the series below will have to evaluate more points.
             icntstart=1
          else
C*** This should be the more common case, you can start looking for
C*** saved integers after nmin
             icntstart=2
          end if
C*** this next term takes care of the f0(nmin)*exp(-alphaprime*nmin)*(integral from a to b - i (-1)^nmin
C*** the f0(nmin) is saved in fn11a(1) etc.
          Z1=outweight*exp(-Bp*nmin)*(C1+icomp*S1)
CCC**** Multiply by outweight before adding to sum from integral which already has these factors
          i11=i11+Z1*fn11a(1)
          i12=i12+Z1*fn12a(1)
          i13=i13+Z1*fn13a(1)
          i22=i22+Z1*fn22a(1)
          i23=i23+Z1*fn23a(1)
          i33=i33+Z1*fn33a(1)
C**** Now for the pprime part
          if (p2fac.ne.0) then
            Z1=outweight*exp(-Bp2*nmin)*(C1+icomp*S1)*p2fac
CCC**** Multiply by outweight before adding to sum from integral which already has these factors
            i11=i11+Z1*fn11b(1)
            i12=i12+Z1*fn12b(1)
            i13=i13+Z1*fn13b(1)
            i22=i22+Z1*fn22b(1)
            i23=i23+Z1*fn23b(1)
            i33=i33+Z1*fn33b(1)
          endif
C**** Now add rest of infinite series, using those saved in fn11a etc. when possible
C          print *,'nmin=',nmin
C          do icnt=1,nsvc
C            print *,'nsave(',icnt,')=',nsave(icnt)
C            print *,'fn11a(icnt)=',fn11a(icnt),
C     + ' fn13a(icnt)=',fn13a(icnt)
C          enddo
C          read *
C%%% Now do convergent series. -i Sum n=nmin+1 to infinity, (-1)^n f0(n) exp(-alpha/f1n1p*n)
C%%% if nperp is real, then this is the only imaginary part for f0 real, or the only real
C%%% part for f0 imaginary  (the n=nmin term is covered already above)
          gtest=S1   ! store -1^(nmin+1) in gtest for use with each series calculation
          DO icnt=1,MAXSERIES
C*** reset nsav2 to all zeros (none saved yet for this rho)
             nsav2(icnt)=0
          ENDDO

******
C---------- Series for inte11 -------------
C
C *** NOTE: for the series parts below, probably no need to go beyond
C  about 5 or 6 digits of accuracy, which is probably all the accuracy
C  attainable with the Gaussian integration anyway.
C
C*** NOTE: with only a few terms in the series, this method of
C  saving and looking up can actually slow down a bit, but for
C  long series, it can help.

          doneseries=0
          Z1=0.0D0
          inte11=0.0D0
          f11a=inte11
          C1=nmin+1.0D0
CC          print *,'inte11C1=',C1
          icnt=icntstart   ! start off at nmin+1 location in nsave
 31       if (nsave(icnt) .lt. int(C1)) then
            icnt=icnt+1
            if (icnt .le. nsvc) then
               go to 31
            end if
          end if
C*** now nsave(icnt) >= int(C1)
          if (lowtemp(ccnt).eq.1) then
            S2=-Bp*C1+alpha(ccnt)
            S3=-Bp2*C1+alpha(ccnt)
          else
            S2=-Bp*C1
            S3=-Bp2*C1
          endif
          if (nopp.eq.1) S3=0
C*** If saved from integral (rounding omegastar), use it
          if (icnt .le. nsvc .and. nsave(icnt) .eq. int(C1)) then
C**** use previously saved value of f0 at the integer value nsave(icnt)
            Z1=S1*(fn11a(icnt)*exp(S2)+
     +  p2fac*fn11b(icnt)*exp(S3))
            if (msg.ge.2) print *,'!!!!Using saved Z1, inte11 C1=',C1
          else  ! must construct f0 at this integer value
C***** NOTE: these sections of special construction for integer values
C*** that were not saved in the course of the calculation of the inner
C*** is now saved for the other 5 series, contrary to previous versions.
C*** switching to low PIN leads to fewer points that were saved before,
C*** and series can go up to hundreds of points, so calculate and save
C*** bessel functions as you go.
C  for f11 -- need J_n+1
C      f12 --      J_n-1 and J_n+1
C      f13 --      J_n and J_n+1
C      f22 --      J_n-1
C      f23 --      J_n-1 and J_n
C      f33 --      J_n
C**** the Bessel functions needed might range from nmin up to thousands
C
C            print *,'Construct f0 for n=',C1
CCC need to calculate fo(C1) -- same as below but for integer xieval
            gcur=C1
C*** gcur is an integer omegastar
            m11=(1.0D0-psq)*gcur*gcur/f1n1psq-1.0D0
            pperp=dsqrt(m11)
C            print *,'m11(pperpsq)=',m11,' pperp=',pperp
            m11=-m11  ! used as -pperp^2 in the matrix elements
CC%%% pperp is the same whether p is p or pprime
            rho=n2b*fb*pperp
            scaled=0
            IFAIL=-1
C compute besselj functions  --- see bessjcomplex.f,
C calculate only J_nu and get J_-nu from the relation
C J_-nu = -1^nu J_nu  for integral nu, since gcur is always integral
            if (msg.ge.2) then
              print *,'(inte11) zbesj with order ',1+gcur,' z=',rho
            endif
COLD NAG            CALL S17DEF(1+gcur,rho,1,'U',CJ,NZ,IFAIL)
            nsav2ind=NINT(gcur-nmin)+2  ! for 1+gcur
            if (nsav2(nsav2ind).eq.1) then
C***** use previously calculated and now saved J_integer(rho)
               CJ(1)=bessjsav(nsav2ind)
               if (bjsavscaled(nsav2ind)) scaled=1
            else
cSm              CALL zbesj(realpart(rho),imagpart(rho),1+gcur,1,1,
              CALL zbesj(dreal(rho),dimag(rho),1+gcur,1,1,
     +                 CJR, CJI, nz,IFAIL)
C*** J_nu where nu is 1+gcur
              CJ(1)=DCMPLX(CJR(1),CJI(1))
              if (IFAIL .ne. 0 .and. IFAIL .ne. 3) then
C**** Before giving up, try calculating scaled version of J_n
C*** because sometimes in these series, it is the large imaginary part
C*** that causes it to fail.
cSm                 CALL zbesj(realpart(rho),imagpart(rho),1+gcur,2,1,
                 CALL zbesj(dreal(rho),dimag(rho),1+gcur,2,1,
     +    CJR, CJI, nz, IFAIL)
                 CJ(1)=DCMPLX(CJR(1),CJI(1))
                 if (IFAIL .ne. 0) then
               print *,'WARNING: infinite series for damping terms ',
     + 'does not seem to be converging. Probably too large imaginary ',
     + 'nperp. nperp=',n2b,' Stopping root finder for this root, ',
     + ' and setting nperp to 1. (~~~ zbesj failed ~~~)'
         print *,'rho=',rho,' nu=1+gcur=',1+gcur,' zbesj result:'
         print *,CJ(1),' IFAIL=',IFAIL,' icnt=',icnt,' Z1=',Z1,
     + ' inte11=',inte11
               n2=dcmplx(1.0D0,0.0D0)
               dispfun=dcmplx(0.0D0,0.0D0)
               return
                 else
                  scaled=1
                 endif
              endif
C**** After testing some scenarios where the root finder started
C*** wandering into heavily damped areas (i.e. Im nperp > Re nperp)
C*** and finding that the infinite series from the damping terms
C*** doesn't converge -- I decided to just pretend it found the root
C*** and abort
              bessjsav(nsav2ind)=CJ(1)
              nsav2(nsav2ind)=1
              if (scaled.eq.1) then
                 bjsavscaled(nsav2ind)=.TRUE.
              else
                 bjsavscaled(nsav2ind)=.FALSE.
              endif
            endif

C  Now CJ(1)=J_1+g(rho)
CCC Now calculate J_-nu from J_-nu = -1 ^ nu J_nu
            if (dmod(1+gcur,2.0D0) .eq. 0) then
              T1=1.0D0
            else
              T1=-1.0D0
            end if
C            print *,'rho(gs): ',rho
C            print *,'nu1(gs): ',-1+gcur,' J(3series) ',CJ
C            print *,'-nu1(gs): ',1-gcur,' J(3series) ',CWRK
C**** Now bessel functions are calculated.
C*** now put together with factor of gcur/f1n1p^2 and J_alpha1 * J_alpha2
            fg11a=m11*CJ(1)*T1
C**** Because of overflow sometimes, I had to move the second
C*** factor of CJ(1) down below
C            if (C1 .gt. 100) then
C              print *,'fg11a: ',fg11a
CCC              if (C1 .gt. 140) stop
C            endif
            if (nopp.eq.0) then
              fg11b=gcur/f1n1p2sq*fg11a
            else
              fg11b=0.0D0
            endif
            fg11a=gcur/f1n1psq*fg11a
            if (scaled .eq. 1) then
C*** need to multiply by a factor of exp(|Im rho|) for each CJ(1)
C*** S2 and S3 should be negative, to help avoid overflow.
              print *,'S2+2*DABS(DIMAG(RHO)) ',S2+2*DABS(DIMAG(RHO))
              print *,'S1=',S1,' CJ(1)=',CJ(1),' S3=',S3
              Z1=S1*fg11a*exp(S2+2*DABS(DIMAG(RHO)))*CJ(1)
              if (p2fac.ne.0) Z1=Z1+S1*p2fac*fg11b*exp(S3+
     +   2*DABS(DIMAG(RHO)))*CJ(1)
              print *,'Z1(from scaled) ',Z1
              print *,'f=',f,' n1=',n1,' te=',MEC2/alpha(1)
              print *,'INTE11'
              if (Z1 .ne. Z1) return
C***** END here
            else
              Z1=S1*(fg11a*exp(S2)*CJ(1)+
     + p2fac*fg11b*exp(S3)*CJ(1))
C              if (C1 .gt. 100) then
C                print *,'Z1(unscaled) ',Z1,' fg11a ',fg11a,
C     + ' Bp*C1 ',Bp*C1,' exp(-Bp*C1) ',exp(-Bp*C1),' p2fac ',
C     + p2fac,' fg11b ',fg11b,' Bp2*C1 ',Bp2*C1,' exp(-Bp2*C1) ',
C     + exp(-Bp2*C1)
C              endif
            endif
              if (msg.ge.2) then
              print *,'(inte11series) S1 ',S1,' fg11a ',fg11a,
     +  ' S2 ',S2,' CJ(1) ',CJ(1),' p2fac ',p2fac,' fg11b ',fg11b,
     +  ' S3 ',S3
              endif
          end if  !end of block to construct f0 at integer n
C          if (msgon.eq.1) then
C            write (*,32) C1,Z1,inte11
C 32         format ('***C1=',f6.1,' Z1=',g12.5,',',g12.5,
C     +  ' inte11=',g12.5,',',g12.5)
C          endif
          inte11=inte11+Z1
C**OLD more restrictive
C  f11a is inte11 before the addition of Z1
          if (inte11 .eq. f11a .or. Z1 .eq. DCMPLX(0.D0,0.D0)
     + .or. ABS(Z1/inte11) .lt. SERIESTOL) then
            doneseries=1
          else
            if (msg.ge.2) then
              print *,'(inte11series) inte11=',inte11,' Z1=',Z1,
     +  ' f11a=',f11a,' C1=',C1
            endif
            if (msg.eq.-1 .and. f.lt.5.07D0 .and. 
cSm     +  realpart(n2).gt.40.0D0) then
     +  dreal(n2).gt.40.0D0) then
               print *,'C1=',C1,' Z1=',Z1,' inte11=',inte11,
     + ' CJ(1)=',CJ(1)
            endif
            f11a=inte11
            C1=C1+1.0D0
            S1=-S1
          end if
C          print *,'inte11end C1=',C1
C          if (C1.gt.100) then
C            print *,'C1=',C1,' Z1=',Z1,' inte11=',inte11
C            print *,' ---doneseries: ',doneseries
CC            read (10) dummy,dummy
C          endif
          if (doneseries .eq. 0) then
            go to 31
          end if
CC**** now inte11 has the series contribution, add to i11 with Gauss-Laguerre weight
CC*** and Maxwellian factor (in variable outweight)
          i11=i11-icomp*inte11*outweight
C          print *,'--Done with inte11 series --'
C          if (C1.gt.100) stop
C---------- Series for inte12 -------------
          S1=gtest  ! restore -1^(nmin+1) to S1
          doneseries=0          
          Z1=0.0D0
          inte12=0.0D0
          f12a=inte12
          C1=nmin+1.0D0
          icnt=icntstart   ! start off at nmin+1 location in nsave
 33       if (icnt .le. nsvc .and. nsave(icnt) .lt. int(C1)) then
            icnt=icnt+1
            if (icnt .le. nsvc) then
               go to 33
            end if
          end if
          if (lowtemp(ccnt).eq.1) then
            S2=-Bp*C1+alpha(ccnt)
            S3=-Bp2*C1+alpha(ccnt)
          else
            S2=-Bp*C1
            S3=-Bp2*C1
          endif
          if (nopp.eq.1) S3=0
C*** now nsave(icnt) >= int(C1)
          scaled=0
          scaled2=0
          if (nsave(icnt) .eq. int(C1)) then
C**** use previously saved value of f0 at the integer value nsave(icnt)
            Z1=S1*(fn12a(icnt)*exp(S2)+
     + p2fac*fn12b(icnt)*exp(S3))
          else  ! must construct f0 at this integer value
C            print *,'Construct f0 for n=',C1
CCC need to calculate fo(C1) -- same as below but for integer xieval
            gcur=C1
            m11=(1.0D0-psq)*gcur*gcur/f1n1psq-1.0D0
            pperp=dsqrt(m11)
C            print *,'m11(pperpsq)=',m11,' pperp=',pperp
            m11=-m11  ! used as -pperp^2 in the matrix elements
CC%%% pperp is the same whether p is p or pprime
            rho=n2b*fb*pperp

            nsav2ind=NINT(gcur-nmin)  ! for -1+gcur
            if (nsav2(nsav2ind).eq.1 .and. nsav2(nsav2ind+2)
     + .eq.1) then
C***** use previously calculated and now saved J_integer(rho)
               CJ(1)=bessjsav(nsav2ind)
               if (bjsavscaled(nsav2ind)) scaled=1
               CJ(3)=bessjsav(nsav2ind+2)
               if (bjsavscaled(nsav2ind+2)) scaled2=1
            else
              IFAIL=0
C compute besselj functions  --- see bessjcomplex.f,
C calculate only J_nu and get J_-nu from the relation
C J_-nu = -1^nu J_nu  for integral nu, since gcur is always integral
CCC***** SPECIAL CASE IF gcur is less than 1, but this will never occur
COLD NAG            CALL S17DEF(-1+gcur,rho,3,'U',CJ,NZ,IFAIL)
C*** for zbesj, after rho and nu, 1 is for unscaled, 3 is for # of bessel func.
cSm              CALL zbesj(realpart(rho),imagpart(rho),-1+gcur,1,3,
              CALL zbesj(dreal(rho),dimag(rho),-1+gcur,1,3,
     +                 CJR, CJI, nz,IFAIL)
              if (IFAIL .ne. 0 .and. IFAIL .ne. 3) then
C**** Before giving up, try calculating scaled version of J_n
C*** because sometimes in these series, it is the large imaginary part
C*** that causes it to fail.
cSm                 CALL zbesj(realpart(rho),imagpart(rho),-1+gcur,2,3,
                 CALL zbesj(dreal(rho),dimag(rho),-1+gcur,2,3,
     +    CJR, CJI, nz, IFAIL)
                 if (IFAIL .ne. 0) then
                print *,'WARNING: BesselJ failed for nu(-1+gcur)=',
     + -1+gcur,' rho=',rho,' nperp=',n2b,' inte12=',inte12,
     + ' Stopping root finder for this root, and setting nperp to 1.'
                n2=dcmplx(1.0D0,0.0D0)
                dispfun=dcmplx(0.0D0,0.0D0)
                return
                 else
                  scaled=1
                  scaled2=1
                 endif
              endif
              DO I=1,3
                 CJ(I)=DCMPLX(CJR(I),CJI(I))
                 bessjsav(nsav2ind+I-1)=CJ(I)
                 nsav2(nsav2ind+I-1)=1
                 if (scaled.eq.1) then
                    bjsavscaled(nsav2ind+I-1)=.TRUE.
                 else
                    bjsavscaled(nsav2ind+I-1)=.FALSE.
                 endif
              END DO
C  Now CJ(1) = J_-1+g(rho), CJ(2) = J_g(rho), CJ(3)=J_1+g(rho)
CCC Now calculate J_-nu from J_-nu = -1 ^ nu J_nu
            endif
            if (dmod(-1+gcur,2.0D0) .eq. 0) then
              T1=1.0D0
            else
              T1=-1.0D0
            end if
            CWRK(1)=CJ(1)*T1
C  Now CWRK(1) = J_1-g(rho) * (-1)^(n-1)
C**** Now bessel functions are calculated.
C*** now put together with factor of gcur/f1n1p^2 and J_alpha1 * J_alpha2
            fg12a=-m11*CWRK(1)
C*** to avoid possible overflow, put CJ(3) down below, after exp(S2)
            if (nopp.eq.0) then
              fg12b=gcur/f1n1p2sq*fg12a
            else
              fg12b=0.0D0
            endif
            fg12a=gcur/f1n1psq*fg12a
            I=scaled+scaled2    ! need either 1 or 2 factors of exp(|Im rho|)
            if (scaled .eq. 1) then
C*** need to multiply by a factor of exp(|Im rho|) for each CJ(1)
C*** S2 and S3 should be negative, to help avoid overflow.
              print *,'S2+I*DABS(DIMAG(RHO)) ',S2+I*DABS(DIMAG(RHO))
              print *,'S1=',S1,' CJ(3)=',CJ(3),' S3=',S3
              Z1=S1*fg12a*exp(S2+I*DABS(DIMAG(RHO)))*CJ(3)
              if (p2fac.ne.0) Z1=Z1+S1*p2fac*fg12b*exp(S3+
     +   I*DABS(DIMAG(RHO)))*CJ(3)
              print *,'Z1(from scaled) ',Z1
              print *,'f=',f,' n1=',n1,' te=',MEC2/alpha(1)
              print *,'INTE12'
              if (Z1 .ne. Z1) return
            else
              Z1=S1*(fg12a*exp(S2)*CJ(3)+
     + p2fac*fg12b*exp(S3)*CJ(3))
            endif
          end if  ! of block to construct Z1 from integer n
          inte12=inte12+Z1
          if (inte12 .eq. f12a .or. Z1 .eq. DCMPLX(0.0D0,0.0D0)
     + .or. ABS(Z1/inte12) .lt. SERIESTOL) then
            doneseries=1
          else
            f12a=inte12
            C1=C1+1.0D0
            S1=-S1
          end if
C          print *,'C1=',C1,' Z1=',Z1,' inte12=',inte12
          if (doneseries .eq. 0) then
            go to 33
          end if
CC**** now inte11 has the series contribution, add to i11 with Gauss-Laguerre weight
CC*** and Maxwellian factor (in variable outweight)
          i12=i12-icomp*inte12*outweight
C---------- Series for inte13 -------------
          S1=gtest  ! restore -1^(nmin+1) to S1
          doneseries=0          
          Z1=0.0D0
          scaled=0
          scaled2=0
          inte13=0.0D0
          f13a=inte13
          C1=nmin+1.0D0
          icnt=icntstart   ! start off at nmin+1 location in nsave
 35       if (icnt .le. nsvc .and. nsave(icnt) .lt. int(C1)) then
            icnt=icnt+1
            if (icnt .le. nsvc) then
               go to 35
            end if
          end if
          if (lowtemp(ccnt).eq.1) then
            S2=-Bp*C1+alpha(ccnt)
            S3=-Bp2*C1+alpha(ccnt)
          else
            S2=-Bp*C1
            S3=-Bp2*C1
          endif
          if (nopp.eq.1) S3=0
C*** now nsave(icnt) >= int(C1)
          if (nsave(icnt) .eq. int(C1)) then
C**** use previously saved value of f0 at the integer value nsave(icnt)
            Z1=S1*(fn13a(icnt)*exp(S2)+
     + p2fac*fn13b(icnt)*exp(S3))
          else  ! must construct f0 at this integer value
C            print *,'Construct f0 for n=',C1
CCC need to calculate fo(C1) -- same as below but for integer xieval
            gcur=C1
            m11=(1.0D0-psq)*gcur*gcur/f1n1psq-1.0D0
            pperp=dsqrt(m11)
C            print *,'m11(pperpsq)=',m11,' pperp=',pperp
            m11=-m11  ! used as -pperp^2 in the matrix elements
CC%%% pperp is the same whether p is p or pprime
            rho=n2b*fb*pperp
            IFAIL=0
C compute besselj functions  --- see bessjcomplex.f,
C calculate only J_nu and get J_-nu from the relation
C J_-nu = -1^nu J_nu  for integral nu, since gcur is always integral
COLD NAG            CALL S17DEF(-1+gcur,rho,3,'U',CJ,NZ,IFAIL)
            nsav2ind=NINT(gcur-nmin)  ! for -1+gcur
            if (nsav2(nsav2ind+1).eq.1 .and. nsav2(nsav2ind+2)
     + .eq.1) then
C***** use previously calculated and now saved J_integer(rho)
               CJ(2)=bessjsav(nsav2ind+1)
               if (bjsavscaled(nsav2ind+1)) scaled=1
               CJ(3)=bessjsav(nsav2ind+2)
               if (bjsavscaled(nsav2ind+2)) scaled2=1
            else
cSm              CALL zbesj(realpart(rho),imagpart(rho),-1+gcur,1,3,
              CALL zbesj(dreal(rho),dimag(rho),-1+gcur,1,3,
     +                 CJR, CJI, nz,IFAIL)
              if (IFAIL .ne. 0 .and. IFAIL .ne. 3) then

C**** Before giving up, try calculating scaled version of J_n
C*** because sometimes in these series, it is the large imaginary part
C*** that causes it to fail.
cSm                 CALL zbesj(realpart(rho),imagpart(rho),-1+gcur,2,3,
                 CALL zbesj(dreal(rho),dimag(rho),-1+gcur,2,3,
     +    CJR, CJI, nz, IFAIL)
                 if (IFAIL .ne. 0) then
                print *,'WARNING: BesselJ failed for nu(-1+gcur)=',
     + -1+gcur,' rho=',rho,' nperp=',n2b,' inte13=',inte13,
     + ' Stopping root finder for this root, and setting nperp to 1.'
                n2=dcmplx(1.0D0,0.0D0)
                dispfun=dcmplx(0.0D0,0.0D0)
                return
                 else
                  scaled=1
                  scaled2=1
                 endif
              endif
              DO I=1,3
                 CJ(I)=DCMPLX(CJR(I),CJI(I))
                 bessjsav(nsav2ind+I-1)=CJ(I)
                 nsav2(nsav2ind+I-1)=1
                 if (scaled.eq.1) then
                    bjsavscaled(nsav2ind+I-1)=.TRUE.
                 else
                    bjsavscaled(nsav2ind+I-1)=.FALSE.
                 endif
              END DO
C  Now CJ(1) = J_-1+g(rho), CJ(2) = J_g(rho), CJ(3)=J_1+g(rho)
CCC Now calculate J_-nu from J_-nu = -1 ^ nu J_nu
            endif
            if (dmod(-1+gcur,2.0D0) .eq. 0) then
              T1=1.0D0
            else
              T1=-1.0D0
            end if
            CWRK(2)=-T1*CJ(2)
C  Now CWRK(2)=(-1)^(gcur)*J_g(rho)=J_-g(rho)
CCC /export/soft/nag/examples/source   
CCCC f77 bessjcomplex.f -lnag -o bessjcomplex.exe
C            print *,'rho(gs): ',rho
C            print *,'nu1(gs): ',-1+gcur,' J(3series) ',CJ
C            print *,'-nu1(gs): ',1-gcur,' J(3series) ',CWRK
C**** Now bessel functions are calculated. Calculate integrand for G(r)
C**** work from the lowest level upward
CC            m11=-pperp*pperp  ! m12 is just -m11 (NOTE: m11 is assigned above)
            m33a=p*gcur/f1n1p    ! needs to be squared after use in next statement
            if (nopp.eq.0) then
              m33b=p2*gcur/f1n1p2
            else
              m33b=0.0D0
            endif
            m13a=m33a*pperp   !m23 is the same as m13
            m13b=m33b*pperp
C*** now matrix elements are defined
C*** now put together with factor of gcur/f1n1p^2 and J_alpha1 * J_alpha2
C*** recall:
C  Now CJ(1) =  J_-1+b(rho), CJ(2) = J_b(rho),    CJ(3)=J_1+b(rho)
C  Now CWRK(1) = J_1-b(rho), CWRK(2)=J_-b(rho), CWRK(3)=J_-1-b(rho)
            fg13a=CWRK(2)
C**** put CJ(3) factor below after exp(-...) to avoid overflow
            if (nopp.eq.0) then
              fg13b=-icomp*m13b*gcur/f1n1p2sq*fg13a
            else
              fg13b=0.0D0
            endif
            fg13a=-icomp*m13a*gcur/f1n1psq*fg13a
            I=scaled+scaled2    ! need either 1 or 2 factors of exp(|Im rho|)
            if (scaled .eq. 1) then
C*** need to multiply by a factor of exp(|Im rho|) for each CJ(1)
C*** S2 and S3 should be negative, to help avoid overflow.
              print *,'S2+I*DABS(DIMAG(RHO)) ',S2+I*DABS(DIMAG(RHO))
              print *,'S1=',S1,' CJ(3)=',CJ(3),' S3=',S3
              Z1=S1*fg13a*exp(S2+I*DABS(DIMAG(RHO)))*CJ(3)
              if (p2fac.ne.0) Z1=Z1+S1*p2fac*fg13b*exp(S3+
     +   I*DABS(DIMAG(RHO)))*CJ(3)
              print *,'Z1(from scaled) ',Z1
              print *,'f=',f,' n1=',n1,' te=',MEC2/alpha(1)
              print *,'INTE13'
              if (Z1 .ne. Z1) return
            else
              Z1=S1*(fg13a*exp(S2)*CJ(3)+
     + p2fac*fg13b*exp(S3)*CJ(3))
            endif
          end if
          inte13=inte13+Z1
          if (inte13 .eq. f13a .or. Z1 .eq. DCMPLX(0.0D0,0.0D0)
     + .or. ABS(Z1/inte13) .lt. SERIESTOL) then
            doneseries=1
          else
            f13a=inte13
            C1=C1+1.0D0
            S1=-S1
          end if
C          print *,'C1=',C1,' Z1=',Z1,' inte13=',inte13
          if (doneseries .eq. 0) then
            go to 35
          end if
CC**** now inte13 has the series contribution, add to i13 with Gauss-Laguerre weight
CC*** and Maxwellian factor (in variable outweight)
          i13=i13-icomp*inte13*outweight
C---------- Series for inte22 -------------
          S1=gtest  ! restore -1^(nmin+1) to S1
          doneseries=0          
          Z1=0.0D0
          scaled=0
          scaled2=0
          inte22=0.0D0
          f22a=inte22
          C1=nmin+1.0D0
          icnt=icntstart   ! start off at nmin+1 location in nsave
 37       if (icnt .le. nsvc .and. nsave(icnt) .lt. int(C1)) then
            icnt=icnt+1
            if (icnt .le. nsvc) then
               go to 37
            end if
          end if
          if (lowtemp(ccnt).eq.1) then
            S2=-Bp*C1+alpha(ccnt)
            S3=-Bp2*C1+alpha(ccnt)
          else
            S2=-Bp*C1
            S3=-Bp2*C1
          endif
          if (nopp.eq.1) S3=0
C*** now nsave(icnt) >= int(C1)
          if (nsave(icnt) .eq. int(C1)) then
C**** use previously saved value of f0 at the integer value nsave(icnt)
            Z1=S1*(fn22a(icnt)*exp(S2)+
     + p2fac*fn22b(icnt)*exp(S3))
          else  ! must construct f0 at this integer value
C            print *,'Construct f0 for n=',C1
CCC need to calculate fo(C1) -- same as below but for integer xieval
            gcur=C1
            m11=(1.0D0-psq)*gcur*gcur/f1n1psq-1.0D0
            pperp=dsqrt(m11)
C            print *,'m11(pperpsq)=',m11,' pperp=',pperp
            m11=-m11  ! used as -pperp^2 in the matrix elements
CC%%% pperp is the same whether p is p or pprime
            rho=n2b*fb*pperp

            nsav2ind=NINT(gcur-nmin)  ! for -1+gcur
            if (nsav2(nsav2ind).eq.1) then
C***** use previously calculated and now saved J_integer(rho)
               CJ(1)=bessjsav(nsav2ind)
               if (bjsavscaled(nsav2ind)) scaled=1
            else
              IFAIL=0
C compute besselj functions  --- see bessjcomplex.f,
C calculate only J_nu and get J_-nu from the relation
C J_-nu = -1^nu J_nu  for integral nu, since gcur is always integral
COLD NAG            CALL S17DEF(-1+gcur,rho,1,'U',CJ,NZ,IFAIL)
cSm              CALL zbesj(realpart(rho),imagpart(rho),-1+gcur,1,1,
              CALL zbesj(dreal(rho),dimag(rho),-1+gcur,1,1,
     +                 CJR, CJI, nz,IFAIL)
              if (IFAIL .ne. 0 .and. IFAIL .ne. 3) then
C**** Before giving up, try calculating scaled version of J_n
C*** because sometimes in these series, it is the large imaginary part
C*** that causes it to fail.
cSm                 CALL zbesj(realpart(rho),imagpart(rho),-1+gcur,2,1,
                 CALL zbesj(dreal(rho),dimag(rho),-1+gcur,2,1,
     +    CJR, CJI, nz, IFAIL)
                 if (IFAIL .ne. 0) then
                print *,'WARNING: BesselJ failed for nu(-1+gcur)=',
     + -1+gcur,' rho=',rho,' nperp=',n2b,' inte22=',inte22,
     + ' Stopping root finder for this root, and setting nperp to 1.'
                n2=dcmplx(1.0D0,0.0D0)
                dispfun=dcmplx(0.0D0,0.0D0)
                return
                 else
                  scaled=1
                 endif
              endif
              CJ(1)=DCMPLX(CJR(1),CJI(1))
C  Now CJ(1) = J_-1+g(rho)
CCC Now calculate J_-nu from J_-nu = -1 ^ nu J_nu
              nsav2(nsav2ind)=1
              bessjsav(nsav2ind)=CJ(1)
              if (scaled.eq.1) then
                 bjsavscaled(nsav2ind)=.TRUE.
              else
                 bjsavscaled(nsav2ind)=.FALSE.
              endif
            endif
            if (dmod(-1+gcur,2.0D0) .eq. 0) then
              T1=1.0D0
            else
              T1=-1.0D0
            end if
            CWRK(1)=CJ(1)*T1
C  Now CWRK(1) = J_1-g(rho)
C**** Now bessel functions are calculated.
C*** now put together with factor of gcur/f1n1p^2 and J_alpha1 * J_alpha2
            fg22a=m11*CWRK(1)
C*** put factor of CJ(1) below to avoid overflow
            if (nopp.eq.0) then
              fg22b=gcur/f1n1p2sq*fg22a
            else
              fg22b=0.0D0
            endif
            fg22a=gcur/f1n1psq*fg22a

            if (scaled .eq. 1) then
C*** need to multiply by a factor of exp(|Im rho|) for each CJ(1)
C*** S2 and S3 should be negative, to help avoid overflow.
              print *,'S2+2*DABS(DIMAG(RHO)) ',S2+2*DABS(DIMAG(RHO))
              print *,'S1=',S1,' CJ(1)=',CJ(1),' S3=',S3
              Z1=S1*fg22a*exp(S2+2*DABS(DIMAG(RHO)))*CJ(1)
              if (p2fac.ne.0) Z1=Z1+S1*p2fac*fg22b*exp(S3+
     +   2*DABS(DIMAG(RHO)))*CJ(1)
              print *,'Z1(from scaled) ',Z1
              print *,'f=',f,' n1=',n1,' te=',MEC2/alpha(1)
              print *,'INTE22'
              if (Z1 .ne. Z1) return
            else
              Z1=S1*(fg22a*exp(S2)*CJ(1)+
     + p2fac*fg22b*exp(S3)*CJ(1))
            endif
          end if
          inte22=inte22+Z1
          if (inte22 .eq. f22a .or. Z1 .eq. DCMPLX(0.0D0,0.0D0)
     + .or. ABS(Z1/inte22) .lt. SERIESTOL) then
            doneseries=1
          else
            f22a=inte22
            C1=C1+1.0D0
            S1=-S1
          end if
C          print *,'C1=',C1,' Z1=',Z1,' inte22=',inte22
          if (doneseries .eq. 0) then
            go to 37
          end if
CC**** now inte22 has the series contribution, add to i22 with Gauss-Laguerre weight
CC*** and Maxwellian factor (in variable outweight)
          i22=i22-icomp*inte22*outweight
C---------- Series for inte23 -------------
          S1=gtest  ! restore -1^(nmin+1) to S1
          doneseries=0          
          Z1=0.0D0
          scaled=0
          scaled2=0
          inte23=0.0D0
          f23a=inte23
          C1=nmin+1.0D0
          icnt=icntstart   ! start off at nmin+1 location in nsave
 39       if (icnt .le. nsvc .and. nsave(icnt) .lt. int(C1)) then
            icnt=icnt+1
            if (icnt .le. nsvc) then
               go to 39
            end if
          end if
          if (lowtemp(ccnt).eq.1) then
            S2=-Bp*C1+alpha(ccnt)
            S3=-Bp2*C1+alpha(ccnt)
          else
            S2=-Bp*C1
            S3=-Bp2*C1
          endif
          if (nopp.eq.1) S3=0
C*** now nsave(icnt) >= int(C1)
          if (nsave(icnt) .eq. int(C1)) then
C**** use previously saved value of f0 at the integer value nsave(icnt)
            Z1=S1*(fn23a(icnt)*exp(S2)+
     + p2fac*fn23b(icnt)*exp(S3))
          else  ! must construct f0 at this integer value
C            print *,'Construct f0 for n=',C1
CCC need to calculate fo(C1) -- same as below but for integer xieval
            gcur=C1
            m11=(1.0D0-psq)*gcur*gcur/f1n1psq-1.0D0
            pperp=dsqrt(m11)
C            print *,'m11(pperpsq)=',m11,' pperp=',pperp
            m11=-m11  ! used as -pperp^2 in the matrix elements
CC%%% pperp is the same whether p is p or pprime
            rho=n2b*fb*pperp

            nsav2ind=NINT(gcur-nmin)  ! for -1+gcur
            if (nsav2(nsav2ind).eq.1 .and. nsav2(nsav2ind+1)
     + .eq.1) then
C***** use previously calculated and now saved J_integer(rho)
               CJ(1)=bessjsav(nsav2ind)
               if (bjsavscaled(nsav2ind)) scaled=1
               CJ(2)=bessjsav(nsav2ind+1)
               if (bjsavscaled(nsav2ind+1)) scaled2=1
            else
              IFAIL=0
C compute besselj functions  --- see bessjcomplex.f,
C calculate only J_nu and get J_-nu from the relation
C J_-nu = -1^nu J_nu  for integral nu, since gcur is always integral
COLD NAG            CALL S17DEF(-1+gcur,rho,2,'U',CJ,NZ,IFAIL)
C*** zbesj, 1 is for unscaled, 2 is for getting back 2 bessel funcs.
cSm              CALL zbesj(realpart(rho),imagpart(rho),-1+gcur,1,2,
              CALL zbesj(dreal(rho),dimag(rho),-1+gcur,1,2,
     +                 CJR, CJI, nz,IFAIL)
              if (IFAIL .ne. 0 .and. IFAIL .ne. 3) then


C**** Before giving up, try calculating scaled version of J_n
C*** because sometimes in these series, it is the large imaginary part
C*** that causes it to fail.
cSm                 CALL zbesj(realpart(rho),imagpart(rho),-1+gcur,2,2,
                 CALL zbesj(dreal(rho),dimag(rho),-1+gcur,2,2,
     +    CJR, CJI, nz, IFAIL)
                 if (IFAIL .ne. 0) then
                print *,'WARNING: BesselJ failed for nu(-1+gcur)=',
     + -1+gcur,' rho=',rho,' nperp=',n2b,' inte23=',inte23,
     + ' Stopping root finder for this root, and setting nperp to 1.'
                n2=dcmplx(1.0D0,0.0D0)
                dispfun=dcmplx(0.0D0,0.0D0)
                return
                 else
                  scaled=1
                  scaled2=1
                 endif
              endif
              DO I=1,2
                 CJ(I)=DCMPLX(CJR(I),CJI(I))
                 bessjsav(nsav2ind+I-1)=CJ(I)
                 nsav2(nsav2ind+I-1)=1
                 if (scaled.eq.1) then
                    bjsavscaled(nsav2ind+I-1)=.TRUE.
                 else
                    bjsavscaled(nsav2ind+I-1)=.FALSE.
                 endif
              END DO
            endif
C  Now CJ(1) = J_-1+g(rho), CJ(2) = J_g(rho)
CCC Now calculate J_-nu from J_-nu = -1 ^ nu J_nu
            if (dmod(-1+gcur,2.0D0) .eq. 0) then
              T1=1.0D0
            else
              T1=-1.0D0
            end if
            CWRK(1)=CJ(1)*T1
C  Now CWRK(1) = J_1-g(rho)
CCC /export/soft/nag/examples/source   
CCCC f77 bessjcomplex.f -lnag -o bessjcomplex.exe
C            print *,'rho(gs): ',rho
C            print *,'nu1(gs): ',-1+gcur,' J(3series) ',CJ
C            print *,'-nu1(gs): ',1-gcur,' J(3series) ',CWRK
C**** Now bessel functions are calculated. Calculate integrand for G(r)
C**** work from the lowest level upward
CC            m11=-pperp*pperp  ! m12 is just -m11 (NOTE: m11 is assigned above)
            m33a=p*gcur/f1n1p    ! needs to be squared after use in next statement
            if (nopp.eq.0) then
              m33b=p2*gcur/f1n1p2
            else
              m33b=0.0D0
            endif
            m13a=m33a*pperp   !m23 is the same as m13
C            print *,'icomp: ',icomp,' m33a: ',m33a,
C     + ' pperp: ',pperp,' m13a=',m13a
            m13b=m33b*pperp
C*** now matrix elements are defined
C*** now put together with factor of gcur/f1n1p^2 and J_alpha1 * J_alpha2
            fg23a=CWRK(1)
C*** put factor of CJ(2) below to avoid overflow
            if (nopp.eq.0) then
              fg23b=-icomp*m13b*gcur/f1n1p2sq*fg23a
            else
              fg23b=0.0D0
            endif
            fg23a=-icomp*m13a*gcur/f1n1psq*fg23a

            I=scaled+scaled2    ! need either 1 or 2 factors of exp(|Im rho|)
            if (scaled .eq. 1) then
C*** need to multiply by a factor of exp(|Im rho|) for each CJ(1)
C*** S2 and S3 should be negative, to help avoid overflow.
              print *,'S2+I*DABS(DIMAG(RHO)) ',S2+I*DABS(DIMAG(RHO))
              print *,'S1=',S1,' CJ(2)=',CJ(2),' S3=',S3
              Z1=S1*fg23a*exp(S2+I*DABS(DIMAG(RHO)))*CJ(2)
              if (p2fac.ne.0) Z1=Z1+S1*p2fac*fg23b*exp(S3+
     +   I*DABS(DIMAG(RHO)))*CJ(2)
              print *,'Z1(from scaled) ',Z1
              print *,'f=',f,' n1=',n1,' te=',MEC2/alpha(1)
              print *,'INTE23'
              if (Z1 .ne. Z1) return
           else
              Z1=S1*(fg23a*exp(S2)*CJ(2)+
     + p2fac*fg23b*exp(S3)*CJ(2))
            endif
          end if
          inte23=inte23+Z1
          if (inte23 .eq. f23a .or. Z1 .eq. DCMPLX(0.0D0,0.0D0)
     + .or. ABS(Z1/inte23) .lt. SERIESTOL) then
            doneseries=1
          else
            f23a=inte23
            C1=C1+1.0D0
            S1=-S1
          end if
C          print *,'C1=',C1,' Z1=',Z1,' inte23=',inte23
          if (doneseries .eq. 0) then
            go to 39
          end if
CC**** now inte23 has the series contribution, add to i23 with Gauss-Laguerre weight
CC*** and Maxwellian factor (in variable outweight)
          i23=i23-icomp*inte23*outweight
C---------- Series for inte33 -------------
          S1=gtest  ! restore -1^(nmin+1) to S1
          doneseries=0          
          Z1=0.0D0
          scaled=0
          scaled2=0
          inte33=0.0D0
          f33a=inte33
          C1=nmin+1.0D0
          icnt=icntstart   ! start off at nmin+1 location in nsave
 41       if (icnt .le. nsvc .and. nsave(icnt) .lt. int(C1)) then
            icnt=icnt+1
            if (icnt .le. nsvc) then
               go to 41
            end if
          end if
          if (lowtemp(ccnt).eq.1) then
            S2=-Bp*C1+alpha(ccnt)
            S3=-Bp2*C1+alpha(ccnt)
          else
            S2=-Bp*C1
            S3=-Bp2*C1
          endif
          if (nopp.eq.1) S3=0
C*** now nsave(icnt) >= int(C1)
          if (nsave(icnt) .eq. int(C1)) then
C**** use previously saved value of f0 at the integer value nsave(icnt)
            Z1=S1*(fn33a(icnt)*exp(S2)+
     + p2fac*fn33b(icnt)*exp(S3))
          else  ! must construct f0 at this integer value
C            print *,'Construct f0 for n=',C1
CCC need to calculate fo(C1) -- same as below but for integer xieval
            gcur=C1
            m11=(1.0D0-psq)*gcur*gcur/f1n1psq-1.0D0
            pperp=dsqrt(m11)
C            print *,'m11(pperpsq)=',m11,' pperp=',pperp
            m11=-m11  ! used as -pperp^2 in the matrix elements
CC%%% pperp is the same whether p is p or pprime
            rho=n2b*fb*pperp

            nsav2ind=NINT(gcur-nmin)  ! for -1+gcur
            if (nsav2(nsav2ind+1).eq.1) then
C***** use previously calculated and now saved J_integer(rho)
               CJ(1)=bessjsav(nsav2ind+1)
               if (bjsavscaled(nsav2ind+1)) scaled=1
            else
              IFAIL=0
C compute besselj functions  --- see bessjcomplex.f,
C calculate only J_nu and get J_-nu from the relation
C J_-nu = -1^nu J_nu  for integral nu, since gcur is always integral
COLD NAG            CALL S17DEF(gcur,rho,1,'U',CJ,NZ,IFAIL)
C*** zbesj, 1 is for unscaled, 1 is for getting back 1 bessel func
cSm              CALL zbesj(realpart(rho),imagpart(rho),gcur,1,1,
              CALL zbesj(dreal(rho),dimag(rho),gcur,1,1,
     +                 CJR, CJI, nz,IFAIL)
              if (IFAIL .ne. 0 .and. IFAIL .ne. 3) then
C**** Before giving up, try calculating scaled version of J_n
C*** because sometimes in these series, it is the large imaginary part
C*** that causes it to fail.
cSm                 CALL zbesj(realpart(rho),imagpart(rho),gcur,2,1,
                 CALL zbesj(dreal(rho),dimag(rho),gcur,2,1,
     +    CJR, CJI, nz, IFAIL)
                 if (IFAIL .ne. 0) then
                print *,'WARNING: BesselJ failed for nu(gcur)=',
     + gcur,' rho=',rho,' nperp=',n2b,' inte33=',inte33,
     + ' Stopping root finder for this root, and setting nperp to 1.'
                n2=dcmplx(1.0D0,0.0D0)
                dispfun=dcmplx(0.0D0,0.0D0)
                return
                 else
                  scaled=1
                 endif
              endif
              CJ(1)=DCMPLX(CJR(1),CJI(1))
C  Now CJ(1) = J_-1+g(rho)
CCC Now calculate J_-nu from J_-nu = -1 ^ nu J_nu
              nsav2(nsav2ind+1)=1
              bessjsav(nsav2ind+1)=CJ(1)
              if (scaled.eq.1) then
                 bjsavscaled(nsav2ind+1)=.TRUE.
              else
                 bjsavscaled(nsav2ind+1)=.FALSE.
              endif
C  Now CJ(1) = J_g(rho)
CCC Now calculate J_-nu from J_-nu = -1 ^ nu J_nu
            endif
            if (dmod(gcur,2.0D0) .eq. 0) then
              T1=1.0D0
            else
              T1=-1.0D0
            end if
            CWRK(1)=CJ(1)*T1
C  Now CWRK(1) = J_-g(rho)
C**** Now bessel functions are calculated.
CC            m11=-pperp*pperp  ! m12 is just -m11 (NOTE: m11 is assigned above)
            m33a=p*gcur/f1n1p    ! needs to be squared after use in next statement
            m33a=m33a*m33a
            if (nopp.eq.0) then
              m33b=p2*gcur/f1n1p2
              m33b=m33b*m33b
            endif
C*** now matrix elements are defined
C*** now put together with factor of gcur/f1n1p^2 and J_alpha1 * J_alpha2
            fg33a=CWRK(1)
C**** put factor of CJ(1) below to avoid overflow
            if (nopp.eq.0) then
              fg33b=m33b*gcur/f1n1p2sq*fg33a
            else
              fg33b=0.0D0
            endif
            fg33a=m33a*gcur/f1n1psq*fg33a
            if (scaled .eq. 1) then
C*** need to multiply by a factor of exp(|Im rho|) for each CJ(1)
C*** S2 and S3 should be negative, to help avoid overflow.
              print *,'S2+2*DABS(DIMAG(RHO)) ',S2+2*DABS(DIMAG(RHO))
              print *,'S1=',S1,' CJ(1)=',CJ(1),' S3=',S3
              Z1=S1*fg33a*exp(S2+2*DABS(DIMAG(RHO)))*CJ(1)
              if (p2fac.ne.0) Z1=Z1+S1*p2fac*fg33b*exp(S3+
     +   2*DABS(DIMAG(RHO)))*CJ(1)
              print *,'Z1(from scaled) ',Z1
              print *,'f=',f,' n1=',n1,' te=',MEC2/alpha(1)
              print *,'INTE33'
              if (Z1 .ne. Z1) return
            else
              Z1=S1*(fg33a*exp(S2)*CJ(1)+
     + p2fac*fg33b*exp(S3)*CJ(1))
            endif
          end if
          inte33=inte33+Z1
          if (inte33 .eq. f33a .or. Z1 .eq. DCMPLX(0.0D0,0.0D0)
     + .or. ABS(Z1/inte33) .lt. SERIESTOL) then
            doneseries=1
          else
            f33a=inte33
            C1=C1+1.0D0
            S1=-S1
          end if
C          print *,'C1=',C1,' Z1=',Z1,' inte33=',inte33
          if (doneseries .eq. 0) then
            go to 41
          end if
CC**** now inte33 has the series contribution, add to i33 with Gauss-Laguerre weight
CC*** and Maxwellian factor (in variable outweight)
          i33=i33-icomp*inte33*outweight
C------------------------- Outer integral term is done ----------
C**** if expfac is 0, then the other parts will be practically zero also, because of
C**** the exp(-Bp*nmin or n which is almost Ap)
        end if  ! of if expfac .ne. 0
CCC%%% for our case: A=l, B=alpha/fb*(1+n1b*p)
C      read *
 5    continue   ! outer integral over pbar and pbarprime

CCC***** now the i11 etc. contain the (sum of ) sigma elements
CC sig11 = sigxx, sig12 = sigxy, sig13 = sigxz, sig22=sigyy, sig23=sigyz, sig33=sigzz
      C1=fb*ompesqom(ccnt)
CC      COMPLEX*16 fg11a,fg11b,fg12a,fg12b,fg13a,fg13b,fg22a,fg22b
CC      COMPLEX*16 fg23a,fg23b,fg33a,fg33b
CCC re-use variables fg11a, fg11b, etc. for use in calculating
CCC normal and shifted sigma components.
      fg11a=(i11+2*i12+i22)/4.0D0  ! in sig11
      fg12a=icomp*(i11-i22)/4.0D0  ! in sig12
      fg13a=(i13+i23)/2.0D0  ! in sig13
      fg22a=(i11-2*i12+i22)/4.0D0  !in sig22
      fg23a=-icomp*(i13-i23)/2.0D0  !in sig23
C  i33 doesn't need a transformation from + and - coordinates
      if (p0(ccnt).eq.0.0D0) then
        sig11=sig11+C1*fg11a
        sig12=sig12+C1*fg12a
        sig13=sig13+C1*fg13a
        sig22=sig22+C1*fg22a
        sig23=sig23+C1*fg23a
        sig33=sig33+C1*i33
      else
C**** Use Lorentz-transformed elements from Bornatici et al
C  1985 PPCF paper, plus transformation from Weiss' coordinates (ky=kperp)
C to more normal where kx=kperp
C         gam0=DSQRT(1.0D0+p0(ccnt)*p0(ccnt))
C         v0=p0(ccnt)/gam0
C         onen1v0=1.0D0-n1*v0
C         n1b=(n1-v0)/onen1v0
        T1=gam0*onen1v0
        S1=T1*T1
        fg33b=n2*p0(ccnt)  !reuse fg33b as nperp*p0
        sig11=sig11+S1*C1*fg11a   !normal: sigyy
        sig12=sig12+S1*C1*fg12a   !normal: sigxy
        sig13=sig13+C1*T1*(fg13a+fg33b*fg12a)    !normal: -sigyz, Weiss: sigxz
        sig22=sig22+S1*C1*fg22a   !normal: sigxx, Weiss: sigyy
        sig23=sig23+C1*T1*(fg23a+fg33b*fg22a)  !normal: sigxz, Weiss: sigyz
        sig33=sig33+C1*(i33+2*fg33b*fg23a+fg33b*fg33b*fg22a)
      endif ! of p0 non zero
      endif  ! of whether there is a non-zero plasma density for this species
 3    continue  ! components (up to 3) of main distribution function

C*** Turn f11a, etc. from sigma (conductivity) to epsilon (dielectric)
C*** by adding identity matrix - 4*pi*ompesqom*f* sigmas
      C1=4*PI   ! the multiplication by ompesqom and fb is done inside the component loop
C**** NOTE: THIS IS IN WEISS' COORDINATES WHERE nvector=(0,nperp,nparallel) i.e. (x,y,z)
CC***** reuse the f11a etc. variables to mean eps_ij  (dielectric tensor elements)
      f11a=1.0D0-C1*sig11  ! eps_xx
      f22a=1.0D0-C1*sig22  ! eps_yy
      f33a=1.0D0-C1*sig33  ! eps_zz
      f12a=-C1*sig12       ! eps_xy
      f13a=-C1*sig13       ! eps_xz
      f23a=-C1*sig23       ! eps_yz
C******* the preceeding are in Weiss' coordinates, where k = (0,k_perp,k_para)
C**** now, the immediately following are for standard output with genray and
C**** in combination with abhay's trubnikov code, where k= (k_perp,0,k_para)
C**** so there are a few rearrangements
C*** Put these elements in eps(1:3,1:3)
      eps(1,1)=f22a   !eps_yy(stix)=eps_xx(weiss)
      eps(1,2)=f12a   !eps_xy(stix)=eps_xy(weiss)
      eps(2,1)=-f12a  !eps_yx(stix)=-eps_xy(stix)=-eps_xy(weiss)
      eps(1,3)=f23a   !eps_xz(stix)=eps_yz(weiss)
      eps(3,1)=f23a   !eps_zx(stix)=eps_xz(stix)=eps_yz(weiss)
      eps(2,2)=f11a   !eps_yy(stix)=eps_xx(weiss)
      eps(2,3)=-f13a  !eps_yz(stix)=-eps_xz(weiss)
      eps(3,2)=f13a   !eps_zy(stix)=-eps_yz(stix)=eps_xz(weiss)
      eps(3,3)=f33a   !eps_zz(stix)=eps_zz(weiss)

      if (msg.ge.2) then
        print *,'### eps_xx=',f11a
        print *,'### eps_xy=',f12a
        print *,'### eps_xz=',f13a
        print *,'### eps_yy=',f22a
        print *,'### eps_yz=',f23a
        print *,'### eps_zz=',f33a
      endif
C*** now place combinations that appear in the determinant of the dispersion tensor in unused variables
      i22=n2*n2           ! nperp^2
      inte12=n2*n1+f23a   ! nperp*npara + eps_yz
      inte22=f22a-n1sq    ! eps_yy - npara^2
      inte11=f11a-i22-n1sq  ! eps_xx-n^2
      inte33=f33a-i22     ! eps_zz - nperp^2

C**** Now calculate determinant of dispersion tensor
      fg11a=inte22*inte33-inte12*inte12
      i11=inte11*fg11a              ! 1st term
      fg12a=f12a*inte33-f13a*inte12
      i11=i11+f12a*fg12a            ! add 2nd term
      fg13a=f12a*inte12-f13a*inte22
      i11=i11-f13a*fg13a            ! subtract 3rd term

C      sigxyz=ompesqom*f*[sigxx,sigxy,sigxz;-sigxy,sigyy,sigyz;-sigxz,sigyz,sigzz]
C      K=[1,0,0;0,1,0;0,0,1]-4*pi*sigxyz
C      Z1=n2*n2+n1sq
C%%% Dispersion tensor D(nperp=n2,npara=n1)
C      Dtens=K+[-Z1,0,0;0,-n1sq,n1*n2;0,n1*n2,-n2*n2]

C%%% Dispersion function is determinant of dispersion tensor
C      dispfun=det(Dtens);

C      FVEC(1)=DBLE(i11)
C      FVEC(2)=DIMAG(i11)
CC Now returning FVEC - real and imaginary parts of determinant
      dispfun=i11
      if (msg.ge.1) print *,'Det D =',i11
      if (dispfun .ne. dispfun) then
C*** This tests for NaN, if so, then return 0 and quit root solving
        print *,'WARNING: Det=NaN. (dispfun=',dispfun,') ',
     + ' Possibly Bessel Function routines ',
     + 'or exponentials failed. Possibly too large imaginary ',
     + 'nperp. nperp=',n2,' Stopping root finder for this root, ',
     + ' and setting nperp to 1.'
        n2=dcmplx(1.0D0,0.0D0)
        dispfun=dcmplx(0.0D0,0.0D0)
      endif

      return
      END  !FUNCTION DISPFUN (begins l. 403)
C%%----- now return determinant of dispersion tensor in dispfun ----
C--------------------------------------------------
CCC************ BESS_IW_COMPLEX *******
CCC E. Nelson-Melby (from C program for bessel functions)

C  This is a series expansion, works for  rho/nu small (less than 1
C  of course, but best for even smaller), and also converges best for
C  fairly large nu.
C
C  Examples: fails when rho is about the same size or larger than nu.
C      nu=2386.2368236, z=(2386.236236,10.2368236)  fails (QNAN, IND)
C      nu=2346.2368, z=(2000.235632,12.230623)  OK
C      nu=153.36262, z=(300.2365,2.326823) FAILS 
C
C  This evaluates all 6 bessel function products needed for the evaluation
C  of the Isaac Weiss fully relativistic electron-frequency range
C  dielectric tensor elements.
C
      subroutine bess_iw_complex(nu,z,A11,A12,A13,A22,A23,A33)
C
C
      COMPLEX*16 z,A11,A12,A13,A22,A23,A33,zsqf,b,oldA33,a
      INTEGER k,kmax
      DOUBLE PRECISION knu,nu,pinu,PI,sinpinu,nup,num,nusq
      INTEGER msg

      DATA PI /3.1415926535897932384626433832795D0/

      common /debug/ msg

C  /*** in the expansion, zr and zi always appears as (zr^2-zi^2)/4 or zr*zi/2  ***/
C      if (debug.eq.1) print *,'IN BESS_IW_COMPLEX nu=',nu,' z=',z
      zsqf=-z*z/4.0D0
C      if (debug.eq.1) print *,'-zsq/4: ',zsqf
      if (nu.EQ.DNINT(nu)) then
         print *,'WARNING: called bess_iw_complex with integer nu '
         print *,'z= ',z,' nu=',nu
         A11=(0.0D0,0.0D0)
         A12=(0.0D0,0.0D0)
         A13=(0.0D0,0.0D0)
         A22=(0.0D0,0.0D0)
         A23=(0.0D0,0.0D0)
         A33=(0.0D0,0.0D0)
         return
      endif

      b=(1.0D0,0.0D0)     !this b_0 in the series
      pinu=PI*nu
      sinpinu=SIN(pinu)   !if this is zero, then this whole series doesn't work
C*** there is another way of dealing with the series with integer nu
      a=b*sinpinu/pinu     !this is b*1/f_0 in the series
C Now a is a_0
      nusq=nu*nu
      nup=1.0D0+nu
      num=1.0D0-nu
      A33=a
      A22=a*nu/num
      A12=a/(num*nup)
      A11=-a*nu/nup
      A23=a/num
      A13=a/nup

      kmax=4000  !was 2000, but for Pentium 4, 4000 is fine
C   /* maximum number of terms to do, to avoid infinite loops */

C*** Do tolerance checking on sum on A33

      DO k=1,kmax
         a=a*zsqf*2.0D0*(2.0D0*k-1.0D0)/(k*(k*k*1.0D0-nusq))
C** k*k*1.0D0 to make sure it is converted to double precision
C*** Now a is the next a_k in the series, using the previous one
         oldA33=A33
         A33=A33+a
         A22=A22+a*(nu+k)/(num+k)
         A12=A12+a*(2.0D0*k+2.0D0)*(2.0D0*k+1.0D0)/
     +  ((k+2.0D0)*(k+1.0D0)*(k+num)*(k+nup))
         A11=A11+a*(k-nu)/(k+nup)
         A23=A23+a*(2.0D0*k+1.0D0)/((k+1.0D0)*(k+num))
         A13=A13+a*(2.0D0*k+1)/((k+1.0D0)*(k+nup))

C        print *,'k=',k,' tr=',tr,' ti=',ti
C  /** Took out tolerance checking -- just let it get too small to add **/

C  /* if adding to sum does not change sum, terms are too small, quit */
         if (oldA33.eq.A33) then
           go to 30
         end if
      ENDDO
 30   if ((msg.ge.2 .and. msg.ne.3) .and. k.ge.kmax) then
        print *,'** WARNING: quit bess_rmseries at k>=kmax: k=',
     + k,' oldA33=',oldA33,' A33=',A33
      end if
C**** now multiply a few elements by extra powers of z/2
      A12=-A12*zsqf   ! z^2/4 is -zsqf
      A23=A23*z/2.0D0
      A13=A13*z/2.0D0
      return
      END
