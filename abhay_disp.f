C**** NOTE: to change resolution of integration grid, change
C**** paramters PIN and/or POUT both here and in iwfdispfuneps, below

C****** Extra parameters that can be assigned before calling:
C****
C**** xi_lo,xi_hi,nxi  (float)

      SUBROUTINE Disp_Ram(T_e,N_parallel,X_e,Y_e,
     +  N_perp,eps,D)

      !implicit double precision (a-h, o-z)
CENM
C** Note: by default, only use Trubnikov method
      IMPLICIT NONE

      DOUBLE PRECISION T_e,N_parallel,X_e,Y_e
      DOUBLE COMPLEX   N_perp,eps(1:3,1:3),D,ci,fr_func
      DOUBLE PRECISION twopi,clight,xi_lo,dxi,arg,xi_hi
      DOUBLE PRECISION cnpa,alpha,omega,vt,besk2
      DOUBLE PRECISION dbsk0e, dbsk1e, besk_a
      INTEGER nxi, kpr, irule0
C
      include 'param.i'
      include 'one.i'
C*** in one.i is also included navg, diff_err, errabs0, errrel0
      EXTERNAL dbsk0e, dbsk1e, besk_a, fr_func
C
C#####################################################################
C    evaluating the relativistic dielectric tensor elements
C#####################################################################
C
C---------------------------------------------------------------------
C
C VERSION 3:  July 19, 2002
C             (getting sophisticated by evaluating the
C              non-relativistic tensor separately)
C             July 26, 2002
C             (fixed the non-relativistic routines which now
C              seem to be working)
C             November, 2002
C             (getting more advanced)
C             December 26, 2002
C             (fixing the integrands to include asymptotic forms)
C             June 5, 2003
C             (incorporated the weiss' way of doing things)
C             June 9, 2003
C             (beginning to incorporate solving of dispersion
C              relation at different positions in a plasma)
C             July 27, 2003
C             (extended to arbitrary number of roots -- use the
C              n_rts to change the number of roots in various
C              routines)
C
C       July 27, 2004 -- Eric Nelson-Melby. Modifications to IMSL calls
C       to use SLATEC routines (freely available, non-commercial). Also
C       put in comments to myself to explain what is going on in code.
C
C       July 28, 2004 -- E.N-M. Added D0 to double constants (G77 otherwise
C         would interpret it as single precision)
C       3 Nov 2004 -- add output of some critical input parameters
C       compile with: g77 %1.f -mpentium -lslatec -o %1.exe
C
C       Nov 18, 2004 - E.N-M. changed all Bessel function K_2(x) to correct
C        expression (sign error before: K_0 - 2/x K_1, should be K_0 + 2/x K_1)
C       USES file inr1 for input to guide what roots to search for
C       in what parameter space.
C *** NOTE: G77 needs the have namelist name right after opening $ (no space)
C OUTPUT FILES:
C      open(5,file='nr1.m',status='replace')
C  non-rel roots.
C      open(6,file='d1.m',status='replace')
C  non-rel, and fully rel roots, warm roots, dielectric elements
C  and keeping track of parameter that varies. Small file.
C      open(7,file='p1.m',status='replace')
C        write (7,101) r*rp,wpe,wce,dsqrt(wpe**2+wce**2),
C     %                1.d0/wce,vte,bt*bt0,den*den0,te*te0
C physical parameters
C      open(8,file='r1.m',status='replace')
C  fully-rel roots.  (MAIN RESULTS OF CODE IN THIS FILE)
CENM use unit 0 for error messages (stderr)
C  this takes care of all umach(-3,9) calls.
C      open(0,file='errmsg1',status='replace')
C      open(10,file='out1',status='replace')
C  integration parameters (sum, limits) (for Weiss method only right now,
C  can result in a LARGE file)
C
C     March, 2006 -- E.N-M. added fr_func_noeps for use in Muller which
C      expects to call a function f(x), without a second parameter.
C---------------------------------------------------------------------
C
C

      common /trub1/ xi_lo,dxi,nxi
CENM 2Sep05 -- diff_err,navg,errabs0,errel0 are all in the common block
C in one.i now (because they can be set by the user)
C      common /trub2/ diff_err,navg
C      common /trerr/ errabs0,errrel0,irule0
      common /trerr/ irule0
C
      common /wave0/ cnpa
      common /para9/ alpha,omega,vt,besk2
C
      common /prnt/  kpr
C
C.... absolute and relative errors and the rule for integration ......
C.. navg = number of intervals for averaging the integration result ..
C.... diff_err = error tolerance between three of these averages .....
C
C      namelist /trint1/ ktrub,xi_lo,xi_hi,nxi
C      namelist /trint2/ errabs0,errrel0,irule0
C      namelist /trint3/ navg,diff_err


      twopi=6.283185307179586476925287d0
      clight=2.99792458d10

C********** DEFAULT VALUES OF OPTIONAL PARAMETERS ***************
C integration limits and subdivisions
      xi_lo=0.0d0
      xi_hi=1.d8
      nxi=1000000
      irule0=6
C absolute and relative errors for integration

cENM050218 I think I have a good balance now between speed and accuracy.
cENM050218 I'll keep looking at this somewhat, making sure that it is not
cENM050218 losing too much accuracy for the EBWs and other O and X mode runs.
cENM050218 Perhaps the code itself could choose what would be most 
cENM050218 appropriate.
cENM050218
cENM050218      errabs0=1.0d-7
cENM050218      errrel0=1.0d-7
cENM050218      irule0=6
cENM050218      navg=25
cENM050218     diff_err=1.0d-6
c
c
      if (relres.eq.1) then
        errabs0=1.0d-4
        errrel0=1.0d-4
        navg=3
        diff_err=0.1d0
C*** if relres is 0, then it was never set, and use the medium option by default
      elseif (relres.eq.2 .or. relres.eq.0) then
        errabs0=1.0d-5
        errrel0=1.0d-5
        navg=12
        diff_err=1.d-3
      elseif (relres.eq.3) then
        errabs0=1.0d-6
        errrel0=1.0d-6
        navg=25
        diff_err=1.d-6
      elseif (relres.eq.4) then
        if (errabs0.eq.0 .and. errrel0.eq.0 .and. navg.eq.0 .and. 
     +    diff_err.eq.0) then
           WRITE(*,*) 'ERROR -- id=14, and relres=4, but errabs0,',
     +       'errrel0,navg, and diff_err are all 0. '
           WRITE(*,*) 'Set them manually in genray.dat'
           STOP
        endif
      endif         

!---------------------------------------------------------
! For use with Abhay Ram's dispersion relation (id=14) parameters to
! control the integration routine for the Trubnikov integral.
!
! relres offers 4 choices for the resolution to use for the relativistic
! dispersion relation using the Trubnikov integral:
! relres=1 low resolution, errabs0=1.d-4, errrel0=1.d-4, navg=3, diff_err=0.1
!       =2 medium res., errabs0=1.d-5, errrel0=1.d-5, navg=12, diff_err=1.d-3
!       =3 high res., errabs0=1.d-6, errrel0=1.d-6, navg=25, diff_err=1.d-6
!       =4 user-defined res., set the following parameters manually:
! The Trubnikov one-dimensional (complex) integral is performed by splitting
! up a region from 0 to 1.d8 into 10^6 pieces, and each piece is integrated
! using the SLATEC adaptive quadrature routine dqag. errabs0 and errrel0 are
! the absolute and relative error tolerances passed directly to dqaq.
! Then the adjacent pieces are compared (it is an oscillatory integrand)
! and using navg number of pieces, when the average difference between them
! are less then diff_err, the integration is presumed finished (Thus it may
! finish long before the upper limit of 1.d8).
!
! errabs0 - absolute error for dqag integration routine
! errrel0 - relative error for dqag integration routine
! navg - number of adjacent integration intervals to use in comparison
! diff_err - error tolerance using navg pieces, when the average difference
!         is less than diff_err, then the integration is done.
! To decide when one should use the low, medium, or high resolution 
! integration, here are some suggestions based on the behavior of the
! Trubnikov integrand: The integrand converges more slowly, and hence
! the resolutions should be set higher, for low electron temperature,
! low (i.e. near zero) magnitude of n_parallel, and for low (near or
! below the fundamental cyclotron frequency) frequency. 
! Examples: n_parallel = -0.05, Te=400 eV, omega/omega_ce=0.4 to 1.2, 
!    it was necessary to use errabs0=1.d-5,errrel0=1.d-5,navg=20,diff_err=1.d-5
!    to be completely converged.  By changing Te to 4000 eV, it was sufficient
!    to use 1.d-4,1.d-4,15 and 1.d-4. 
!   An easy case: n_parallel=0.3, Te=7 keV, omega/omega_ce=2.4 to 2.7, 
!    complete convergence already at errabs0=1.d-4,errrel0=1.d-4,navg=2,
!    diff_err=0.5
!   An intermediate case: n_parallel=0.1, omega/omega_ce=1.0, Te=300 eV
!    errabs0=1.d-5,errrel0=1.d-5,navg=12,diff_err=1.d-3 was OK.
!--------------------------------------------------------
!C****************************************************************

C...... the integration grid for the trubnikov form (if needed) ......
C
      dxi=(xi_hi-xi_lo)/dfloat(nxi)
      ci=dcmplx(0.d0,1.d0)
C
C*** Need to set up te, ompesqom, p0, f, n1
C in inr1: alpha,omega,vt,besk2
C alpha = omega_pe^2/omega_ce^2 = omega_pe^2/omega^2 * omega^2/omega_ce^2
C omega = omega/omega_ce, where omega is 2*pi*frequency (Hz)
C vt = thermal velocity / c0
      alpha=X_e/Y_e**2  ! X_e = omega_pe^2/omega^2
C T_e is in keV
      vt=4.19396657d7/2.99792458d10*dsqrt(T_e*1.d3)

      cnpa=N_parallel
      omega=1.0D0/Y_e    ! Y_e is omega_ce/omega

      arg=1./vt**2
CENM  if arg.lt.50 is when Te is more than about 10.22 keV
      if (arg.lt.50.d0) then
CENM        call umach(-3,9)
CENM IMSL->SLATEC dbsk0e and dbsk1e into dbsk0e and dbsk1e
        besk2=dbsk0e(arg)+2.*dbsk1e(arg)/arg
CENM just let SLATEC handle errors        call erset(0,1,1)
      else
C if lower temp than 10 keV, use asymptotic version of K2
        besk2=besk_a(2,arg)
      end if

      D=fr_func(N_perp,eps)

      END


C---------------------------------------------------------------------
C
C............ the fully-relativistic dispersion function .............
C
      double complex function fr_func(roots,eps)
C
CENM      !implicit double precision (a-h, o-z)
      IMPLICIT NONE
C
      double complex roots,eps(1:3,1:3)
      double complex d11,d12,d13,d22,d23,d33
      double complex x11,x12,x13,x22,x23,x33
C
      double complex ep(6),cnpe
      double precision cnpa
C
      common /wave0/ cnpa
      common /wave1/ cnpe
C
      cnpe=roots
C      print *,'^^^^ in fr_func: cnpa=',cnpa,' cnpe=',cnpe
C
C*ENM- by default, only use Trubnikov method
C      if (ktrub.eq.1) then
      call trub(ep)
C      else
C        call weiss(ep)
C      end if
C
      x11=ep(1)
      x12=ep(2)
      x13=ep(3)
      x22=ep(4)
      x23=ep(5)
      x33=ep(6)

C*** Put these elements in eps(1:3,1:3)
      eps(1,1)=1.0D0+x11
      eps(1,2)=x12
      eps(2,1)=-x12
      eps(1,3)=x13
      eps(3,1)=x13
      eps(2,2)=1.0D0+x22
      eps(2,3)=x23
      eps(3,2)=-x23
      eps(3,3)=1.0D0+x33

C
C........ the dielectric tensor elements (satisfying D.E = 0) ........
C
      d11=1.0D0-cnpa**2+x11
      d12=x12
      d13=cnpe*cnpa+x13
      d22=1.0D0-cnpe**2-cnpa**2+x22
      d23=x23
      d33=1.0D0-cnpe**2+x33
C
      fr_func=d11*(d22*d33+d23*d23)+
     %        d12*(d23*d13+d12*d33)+
     %        d13*(d12*d23-d22*d13)
C
C
C      print *,'*fr_func=',fr_func,' cnpe=',cnpe
      return
      end
C
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
C............ the fully-relativistic dispersion function .............
C
CENM -- Same as fr_func, but without the eps in the call,
C       in order to work with iabsorp=12 which uses complex
C       root finding (Muller algorithm).
      double complex function fr_func_noeps(roots)
C
CENM      !implicit double precision (a-h, o-z)
      IMPLICIT NONE
C
      double complex roots
      double complex d11,d12,d13,d22,d23,d33
      double complex x11,x12,x13,x22,x23,x33
C
      double complex ep(6),cnpe
      double precision cnpa
C
      common /wave0/ cnpa
      common /wave1/ cnpe
C
      cnpe=roots
C      print *,'^^^^ in fr_func_noeps: cnpa=',cnpa,' cnpe=',cnpe
C
C*ENM- by default, only use Trubnikov method
C      if (ktrub.eq.1) then
      call trub(ep)
C      else
C        call weiss(ep)
C      end if
C
      x11=ep(1)
      x12=ep(2)
      x13=ep(3)
      x22=ep(4)
      x23=ep(5)
      x33=ep(6)

C
C........ the dielectric tensor elements (satisfying D.E = 0) ........
C
      d11=1.0D0-cnpa**2+x11
      d12=x12
      d13=cnpe*cnpa+x13
      d22=1.0D0-cnpe**2-cnpa**2+x22
      d23=x23
      d33=1.0D0-cnpe**2+x33
C
      fr_func_noeps=d11*(d22*d33+d23*d23)+
     %        d12*(d23*d13+d12*d33)+
     %        d13*(d12*d23-d22*d13)
C
C
C      print *,'*fr_func=',fr_func,' cnpe=',cnpe
      return
      end
C
C---------------------------------------------------------------------

C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
C..... evaluating the trubnikov form of the relativistic tensor ......
C
      subroutine trub(epsilon)
C
      implicit double precision (a-h, o-z)
C
      double precision stot(3),sumr(0:1000)
C
      double complex dcmplx,csum
C
      double complex epsilon(6),co(0:1)
      double precision work
      integer ier,iwork,last,lenw,limit,neval
CENM this is based on limit=500, lenw>4*limit
      dimension iwork(500),work(2500)
C
      include 'param.i'
      include 'one.i'

      common /trub1/ xi_lo,dxi,nxi
CENM 2Sep05 -- diff_err,navg,errabs0,errel0 are all in the common block
C in one.i now (because they can be set by the user)
C      common /trub2/ diff_err,navg
C      common /trerr/ errabs0,errrel0,irule0
      common /trerr/ irule0
      common /trubc/ kel,kreal
      common /prnt/  kpr
C
CENM      external dqdag,ftrub,umach,erset
      external dqag,ftrub
C
      co(0)=dcmplx(0.d0,1.d0)
      co(1)=dcmplx(1.d0,0.d0)
C
      do 1 j=1,6
        kel=j
        csum=dcmplx(0.d0,0.d0)
C
        do 2 k=0,1
          kreal=k
C
C.. ntr keeps track of the number of averages (3) that are compared ..
C
          ntr=0
C
          xi_0=xi_lo
          xi_1=xi_0+dxi
C
          sum=0.d0
C
          do 3 i=1,nxi
C
            errabs=errabs0
            errrel=errrel0
            irule=irule0
C
CENM turn dqdaq (IMSL) into similar dqag (SLATEC)
CENM            call umach(-3,9)
C            call dqdag(ftrub,xi_0,xi_1,errabs,errrel,irule,
C     %                 t_real,errest)
            limit = 500
C limit is maximum number of subintervals in (xi_0,xi_1)
CENM  IMSL dqdag uses 500 for max number of subintervals
            lenw = limit*4
            call dqag(ftrub,xi_0,xi_1,errabs,errrel,irule,
     *         t_real,errest,neval,ier,limit,lenw,last,iwork,work)

CENM let dqag do error handling            call erset(0,1,1)
C
            sum=sum+t_real
C
            jj=mod(i-1,navg)
            sumr(jj)=sum
            if (jj.eq.navg-1) then
              ntr=ntr+1
              stot(ntr)=0.d0
              do 4 ij=0,navg-1
                stot(ntr)=stot(ntr)+sumr(ij)
    4         continue
              stot(ntr)=stot(ntr)/dfloat(navg)
              if (ntr.eq.3) then
                diff1=dabs(dabs(stot(1))-dabs(stot(2)))
                diff2=dabs(dabs(stot(2))-dabs(stot(3)))
                if (diff1.lt.diff_err.and.diff2.lt.diff_err) then
                  sum_real=stot(3)
                  go to 5
                else
                  stot(1)=stot(2)
                  stot(2)=stot(3)
                  ntr=2
                end if
              end if
            end if
C
            xi_0=xi_1
            xi_1=xi_0+dxi
C
    3     continue
C
    5     continue
C
          if (kpr.eq.1) then
            if (kreal.eq.0) then
              write (11,201) i,sum
  201         format("imag integral used ",i9,
     %               "steps to get ",g18.10)
            else
              write (11,202) i,sum
  202         format("real integral used ",i9,
     %               "steps to get ",g18.10)
            end if
          end if
C
          csum=csum+co(k)*sum
    2   continue
C
        epsilon(kel)=csum
C
    1 continue
C
C
      return
      end
C
C---------------------------------------------------------------------
C
      double precision function ftrub(xi)
C
      implicit double precision (a-h, o-z)
C
      double complex dcmplx,ci,factor
      double complex tx_11,tx_12,tx_13
      double complex tx_22,tx_23,tx_33
C
      common /para9/ alpha,omega,vt,besk2
      common /trubc/ kel,kreal
C
      external tx_11,tx_12,tx_13
      external tx_22,tx_23,tx_33
C
      ci=dcmplx(0.d0,1.d0)
C
C...... factor is the multiplier in the dielectric tensor form .......
C.............. factor = i * {(wpe/wce)**2} * (wce/w0) ...............
C
      factor=ci*alpha/omega
C
      if (kreal.eq.0) then
        if (kel.eq.1) ftrub=dimag(factor*tx_11(xi))
        if (kel.eq.2) ftrub=dimag(factor*tx_12(xi))
        if (kel.eq.3) ftrub=dimag(factor*tx_13(xi))
        if (kel.eq.4) ftrub=dimag(factor*tx_22(xi))
        if (kel.eq.5) ftrub=dimag(factor*tx_23(xi))
        if (kel.eq.6) ftrub=dimag(factor*tx_33(xi))
      else
        if (kel.eq.1) ftrub=dreal(factor*tx_11(xi))
        if (kel.eq.2) ftrub=dreal(factor*tx_12(xi))
        if (kel.eq.3) ftrub=dreal(factor*tx_13(xi))
        if (kel.eq.4) ftrub=dreal(factor*tx_22(xi))
        if (kel.eq.5) ftrub=dreal(factor*tx_23(xi))
        if (kel.eq.6) ftrub=dreal(factor*tx_33(xi))
      end if
C
C
      return
      end
C
C---------------------------------------------------------------------
C
      double complex function tx_11(xi)
C
      implicit double precision (a-h, o-z)
C
      double complex dcmplx,cdexp,cdsqrt
      double complex cnpe
      double complex ci,r,sr,t2
      double complex ck2,ck3,factor
      double complex ce,cfac,cbesk_a
C
      double complex ckbs(2)
      double precision cyr(2),cyi(2)
      integer nz,ierr
C
      common /para9/ alpha,omega,vt,besk2
      common /wave0/ cnpa
      common /wave1/ cnpe
C
CENM      external dcbks,umach,erset
      external zbesk
C
      ci=dcmplx(0.d0,1.d0)
C
      twopi=6.283185307179586476925287d0
      xim=dmod(xi,twopi)
C
C... r is the actual r multiplied by (vte/c)**4 and K_2({c/vt}**2) ...
C.... this takes care of the normalization factors in the tensor .....
C... sr is the square-root of r (argument of the Bessel functions) ...
C............ dom = omega*(vte/c)**2 = (w/wce)*(vte/c)**2 ............
C
      dom=omega*vt**2
      r=(1.d0-ci*xi*dom)**2+2.d0*(cnpe*dom)**2*(1.d0-dcos(xim))+
     %  (cnpa*dom*xi)**2
CC      r=(1./vt**2-ci*xi*omega)**2+
CC     %  2.*(cnpe*omega)**2*(1.-dcos(xim))+
CC     %  (cnpa*omega*xi)**2
      sr=cdsqrt(r)/vt**2
      r=r*besk2
C
      if (dreal(sr).lt.50.d0) then
C
        numb=2
        xnu=2.d0
CENM        call umach(-3,9)
CENM change IMSL routine dcbks into SLATEC routine zbesk
C        call dcbks(xnu,sr,numb,ckbs)
        call zbesk(dreal(sr),dimag(sr),xnu,1,numb,
     +    cyr,cyi,nz,ierr)
        do i=1,2
           ckbs(i)=dcmplx(cyr(i),cyi(i))
        enddo

CENM        call erset(0,1,1)
C
        ck2=ckbs(1)*dexp(1.d0/vt**2)
        ck3=ckbs(2)*dexp(1.d0/vt**2)
C
      else
C
        ce=1.d0/vt**2-sr
        cfac=cdexp(ce)
        ck2=cbesk_a(2,sr)*cfac
        ck3=cbesk_a(3,sr)*cfac
C
      end if
C
      t1=dcos(xim)
C
      co=omega**2
C
      t2=co*(cnpe*dsin(xim))**2
C
CC      tx_11=ck2/(r*vt**4*besk2)*t1-ck3/(r*sr*vt**4*besk2)*t2
      tx_11=ck2/r*t1-ck3/(r*sr)*t2
C
C
      return
      end
C
C---------------------------------------------------------------------
C
      double complex function tx_12(xi)
C
      implicit double precision (a-h, o-z)
C
      double complex dcmplx,cdexp,cdsqrt
      double complex cnpe
      double complex ci,r,sr,t2
      double complex ck2,ck3
      double complex ce,cfac,cbesk_a
C
      double complex ckbs(2)
      double precision cyr(2),cyi(2)
      integer nz,ierr
C
      common /para9/ alpha,omega,vt,besk2
      common /wave0/ cnpa
      common /wave1/ cnpe
C
CENM      external dcbks,umach,erset
      external zbesk
C
      ci=dcmplx(0.d0,1.d0)
C
      twopi=6.283185307179586476925287d0
      xim=dmod(xi,twopi)
C
C... r is the actual r multiplied by (vte/c)**4 and K_2({c/vt}**2) ...
C.... this takes care of the normalization factors in the tensor .....
C... sr is the square-root of r (argument of the Bessel functions) ...
C............ dom = omega*(vte/c)**2 = (w/wce)*(vte/c)**2 ............
C
      dom=omega*vt**2
      r=(1.d0-ci*xi*dom)**2+2.d0*(cnpe*dom)**2*(1.d0-dcos(xim))+
     %  (cnpa*dom*xi)**2
CC      r=(1./vt**2-ci*xi*omega)**2+
CC     %  2.*(cnpe*omega)**2*(1.-dcos(xim))+
CC     %  (cnpa*omega*xi)**2
      sr=cdsqrt(r)/vt**2
      r=r*besk2
C
      if (dreal(sr).lt.50.d0) then
C
        numb=2
        xnu=2.d0
CENM        call umach(-3,9)
CENM change IMSL routine dcbks into SLATEC routine zbesk
C        call dcbks(xnu,sr,numb,ckbs)
        call zbesk(dreal(sr),dimag(sr),xnu,1,numb,
     +    cyr,cyi,nz,ierr)
        do i=1,2
           ckbs(i)=dcmplx(cyr(i),cyi(i))
        enddo
CENM        call erset(0,1,1)
C
        ck2=ckbs(1)*dexp(1.d0/vt**2)
        ck3=ckbs(2)*dexp(1.d0/vt**2)
C
      else
C
        ce=1.d0/vt**2-sr
        cfac=cdexp(ce)
        ck2=cbesk_a(2,sr)*cfac
        ck3=cbesk_a(3,sr)*cfac
C
      end if
C
      t1=-dsin(xim)
C
      co=omega**2
C
      t2=-co*cnpe**2*dsin(xim)*(1.d0-dcos(xim))
C
CC      tx_12=ck2/(r*vt**4*besk2)*t1-ck3/(r*sr*vt**4*besk2)*t2
      tx_12=ck2/r*t1-ck3/(r*sr)*t2
C
C
      return
      end
C
C---------------------------------------------------------------------
C
      double complex function tx_13(xi)
C
      implicit double precision (a-h, o-z)
C
      double complex dcmplx,cdexp,cdsqrt
      double complex cnpe
      double complex ci,r,sr,t2
      double complex ck3
      double complex ce,cfac,cbesk_a
C
      double complex ckbs(2)
      double precision cyr(2),cyi(2)
      integer nz,ierr
C
      common /para9/ alpha,omega,vt,besk2
      common /wave0/ cnpa
      common /wave1/ cnpe
C
CENM      external dcbks,umach,erset
      external zbesk
C
      ci=dcmplx(0.d0,1.d0)
C
      twopi=6.283185307179586476925287d0
      xim=dmod(xi,twopi)
C
C... r is the actual r multiplied by (vte/c)**4 and K_2({c/vt}**2) ...
C.... this takes care of the normalization factors in the tensor .....
C... sr is the square-root of r (argument of the Bessel functions) ...
C............ dom = omega*(vte/c)**2 = (w/wce)*(vte/c)**2 ............
C
      dom=omega*vt**2
      r=(1.d0-ci*xi*dom)**2+2.d0*(cnpe*dom)**2*(1.d0-dcos(xim))+
     %  (cnpa*dom*xi)**2
CC      r=(1./vt**2-ci*xi*omega)**2+
CC     %  2.*(cnpe*omega)**2*(1.-dcos(xim))+
CC     %  (cnpa*omega*xi)**2
      sr=cdsqrt(r)/vt**2
      r=r*besk2
C
      if (dreal(sr).lt.50.d0) then
C
        numb=2
        xnu=2.d0
CENM        call umach(-3,9)
CENM change IMSL routine dcbks into SLATEC routine zbesk
C        call dcbks(xnu,sr,numb,ckbs)
        call zbesk(dreal(sr),dimag(sr),xnu,1,numb,
     +    cyr,cyi,nz,ierr)
        do i=1,2
           ckbs(i)=dcmplx(cyr(i),cyi(i))
        enddo
CENM        call erset(0,1,1)
C
        ck2=ckbs(1)*dexp(1.d0/vt**2)
        ck3=ckbs(2)*dexp(1.d0/vt**2)
C
      else
C
        ce=1.d0/vt**2-sr
        cfac=cdexp(ce)
        ck2=cbesk_a(2,sr)*cfac
        ck3=cbesk_a(3,sr)*cfac
C
      end if
C
      co=omega**2
C
      t2=co*cnpe*cnpa*xi*dsin(xim)
C
CC      tx_13=-ck3/(r*sr*vt**4*besk2)*t2
      tx_13=-ck3/(r*sr)*t2
C
C
      return
      end
C
C---------------------------------------------------------------------
C
      double complex function tx_22(xi)
C
      implicit double precision (a-h, o-z)
C
      double complex dcmplx,cdexp,cdsqrt
      double complex cnpe
      double complex ci,r,sr,t2
      double complex ck2,ck3
      double complex ce,cfac,cbesk_a
C
      double complex ckbs(2)
      double precision cyr(2),cyi(2)
      integer nz,ierr
C
      common /para9/ alpha,omega,vt,besk2
      common /wave0/ cnpa
      common /wave1/ cnpe
C
CENM      external dcbks,umach,erset
      external zbesk
C
      ci=dcmplx(0.d0,1.d0)
C
      twopi=6.283185307179586476925287d0
      xim=dmod(xi,twopi)
C
C... r is the actual r multiplied by (vte/c)**4 and K_2({c/vt}**2) ...
C.... this takes care of the normalization factors in the tensor .....
C... sr is the square-root of r (argument of the Bessel functions) ...
C............ dom = omega*(vte/c)**2 = (w/wce)*(vte/c)**2 ............
C
      dom=omega*vt**2
      r=(1.d0-ci*xi*dom)**2+2.d0*(cnpe*dom)**2*(1.d0-dcos(xim))+
     %  (cnpa*dom*xi)**2
CC      r=(1./vt**2-ci*xi*omega)**2+
CC     %  2.*(cnpe*omega)**2*(1.-dcos(xim))+
CC     %  (cnpa*omega*xi)**2
      sr=cdsqrt(r)/vt**2
      r=r*besk2
C
      if (dreal(sr).lt.50.d0) then
C
        numb=2
        xnu=2.d0
CENM        call umach(-3,9)
CENM change IMSL routine dcbks into SLATEC routine zbesk
C        call dcbks(xnu,sr,numb,ckbs)
        call zbesk(dreal(sr),dimag(sr),xnu,1,numb,
     +    cyr,cyi,nz,ierr)
        do i=1,2
           ckbs(i)=dcmplx(cyr(i),cyi(i))
        enddo
CENM        call erset(0,1,1)
C
        ck2=ckbs(1)*dexp(1.d0/vt**2)
        ck3=ckbs(2)*dexp(1.d0/vt**2)
C
      else
C
        ce=1.d0/vt**2-sr
        cfac=cdexp(ce)
        ck2=cbesk_a(2,sr)*cfac
        ck3=cbesk_a(3,sr)*cfac
C
      end if
C
      t1=dcos(xim)
C
      co=omega**2
C
      t2=-co*cnpe**2*(1.d0-dcos(xim))**2
C
CC      tx_22=ck2/(r*vt**4*besk2)*t1-ck3/(r*sr*vt**4*besk2)*t2
      tx_22=ck2/r*t1-ck3/(r*sr)*t2
C
C
      return
      end
C
C---------------------------------------------------------------------
C
      double complex function tx_23(xi)
C
      implicit double precision (a-h, o-z)
C
      double complex dcmplx,cdexp,cdsqrt
      double complex cnpe
      double complex ci,r,sr,t2
      double complex ck3
      double complex ce,cfac,cbesk_a
C
      double complex ckbs(2)
      double precision cyr(2),cyi(2)
      integer nz,ierr
C
      common /para9/ alpha,omega,vt,besk2
      common /wave0/ cnpa
      common /wave1/ cnpe
C
CENM      external dcbks,umach,erset
      external zbesk
C
      ci=dcmplx(0.d0,1.d0)
C
      twopi=6.283185307179586476925287d0
      xim=dmod(xi,twopi)
C
C... r is the actual r multiplied by (vte/c)**4 and K_2({c/vt}**2) ...
C.... this takes care of the normalization factors in the tensor .....
C... sr is the square-root of r (argument of the Bessel functions) ...
C............ dom = omega*(vte/c)**2 = (w/wce)*(vte/c)**2 ............
C
      dom=omega*vt**2
      r=(1.d0-ci*xi*dom)**2+2.d0*(cnpe*dom)**2*(1.d0-dcos(xim))+
     %  (cnpa*dom*xi)**2
CC      r=(1./vt**2-ci*xi*omega)**2+
CC     %  2.*(cnpe*omega)**2*(1.-dcos(xim))+
CC     %  (cnpa*omega*xi)**2
      sr=cdsqrt(r)/vt**2
      r=r*besk2
C
      if (dreal(sr).lt.50.d0) then
C
        numb=2
        xnu=2.d0
CENM        call umach(-3,9)
CENM change IMSL routine dcbks into SLATEC routine zbesk
C        call dcbks(xnu,sr,numb,ckbs)
        call zbesk(dreal(sr),dimag(sr),xnu,1,numb,
     +    cyr,cyi,nz,ierr)
        do i=1,2
           ckbs(i)=dcmplx(cyr(i),cyi(i))
        enddo
CENM        call erset(0,1,1)
C
        ck2=ckbs(1)*dexp(1.d0/vt**2)
        ck3=ckbs(2)*dexp(1.d0/vt**2)
C
      else
C
        ce=1.d0/vt**2-sr
        cfac=cdexp(ce)
        ck2=cbesk_a(2,sr)*cfac
        ck3=cbesk_a(3,sr)*cfac
C
      end if
C
      co=omega**2
C
      t2=co*cnpe*cnpa*xi*(1.d0-dcos(xim))
C
CC      tx_23=-ck3/(r*sr*vt**4*besk2)*t2
      tx_23=-ck3/(r*sr)*t2
C
C
      return
      end
C
C---------------------------------------------------------------------
C
      double complex function tx_33(xi)
C
      implicit double precision (a-h, o-z)
C
      double complex dcmplx,cdexp,cdsqrt
      double complex cnpe
      double complex ci,r,sr,t2
      double complex ck2,ck3
      double complex ce,cfac,cbesk_a
C
      double complex ckbs(2)
      double precision cyr(2),cyi(2)
      integer nz,ierr
C
      common /para9/ alpha,omega,vt,besk2
      common /wave0/ cnpa
      common /wave1/ cnpe
C
CENM      external dcbks,umach,erset
      external zbesk
C
      ci=dcmplx(0.d0,1.d0)
C
      twopi=6.283185307179586476925287d0
      xim=dmod(xi,twopi)
C
C... r is the actual r multiplied by (vte/c)**4 and K_2({c/vt}**2) ...
C.... this takes care of the normalization factors in the tensor .....
C... sr is the square-root of r (argument of the Bessel functions) ...
C............ dom = omega*(vte/c)**2 = (w/wce)*(vte/c)**2 ............
C
      dom=omega*vt**2
      r=(1.d0-ci*xi*dom)**2+2.d0*(cnpe*dom)**2*(1.d0-dcos(xim))+
     %  (cnpa*dom*xi)**2
CC      r=(1./vt**2-ci*xi*omega)**2+
CC     %  2.*(cnpe*omega)**2*(1.-dcos(xim))+
CC     %  (cnpa*omega*xi)**2
      sr=cdsqrt(r)/vt**2
      r=r*besk2
C
      if (dreal(sr).lt.50.d0) then
C
        numb=2
        xnu=2.d0
CENM        call umach(-3,9)
CENM change IMSL routine dcbks into SLATEC routine zbesk
C        call dcbks(xnu,sr,numb,ckbs)
        call zbesk(dreal(sr),dimag(sr),xnu,1,numb,
     +    cyr,cyi,nz,ierr)
        do i=1,2
           ckbs(i)=dcmplx(cyr(i),cyi(i))
        enddo
CENM        call erset(0,1,1)
C
        ck2=ckbs(1)*dexp(1.d0/vt**2)
        ck3=ckbs(2)*dexp(1.d0/vt**2)
C
      else
C
        ce=1.d0/vt**2-sr
        cfac=cdexp(ce)
        ck2=cbesk_a(2,sr)*cfac
        ck3=cbesk_a(3,sr)*cfac
C
      end if
C
      co=omega**2
C
      t2=co*(cnpa*xi)**2
C
CC      tx_33=ck2/(r*vt**4*besk2)-ck3/(r*sr*vt**4*besk2)*t2
      tx_33=ck2/r-ck3/(r*sr)*t2
C
C
      return
      end
C
C---------------------------------------------------------------------
C
C... asymptotic form for the modified Bessel function (second kind) ...
C.......... this is for real argument and is exp(x)*K(nu,x) ..........
C
      double precision function besk_a(nu,z)
C
      implicit none !double precision (a-h,o-z)
C
      integer nu
      real*8 mu,fac
      real*8 z
C
C.......................... fac=sqrt(pi/2) ...........................
C
      fac=1.253314137315500251207883D0
C
CENM  See A&S 9.7.2
      mu=dfloat(4*nu**2)
      besk_a=fac/dsqrt(z)*(1.+(mu-1.d0)/(8.d0*z)+
     %             (mu-1.d0)*(mu-9.d0)/(128.d0*z**2)+
     %             (mu-1.d0)*(mu-9.d0)*(mu-25.d0)/(3072.d0*z**3))
C
C
      return
      end
C
C---------------------------------------------------------------------
C
C... asymptotic form for the modified Bessel function (second kind) ...
C........ this is for complex argument and is exp(x)*K(nu,x) .........
C
      double complex function cbesk_a(nu,z)
C
      implicit none !double precision (a-h,o-z)
      integer nu
      real*8 mu,fac
      double complex z !,cdsqrt
C
C.......................... fac=sqrt(pi/2) ...........................
C
      fac=1.253314137315500251207883D0
C
      mu=dfloat(4*nu**2)
      cbesk_a=fac/cdsqrt(z)*(1.+(mu-1.d0)/(8.d0*z)+
     %            (mu-1.d0)*(mu-9.d0)/(128.d0*z**2)+
     %            (mu-1.d0)*(mu-9.d0)*(mu-25.d0)/(3072.d0*z**3))
C
C
      return
      end
C
C---------------------------------------------------------------------

      subroutine 
     &relativist_absorp_projection_method_det
     &(T_e,N_parallel,X_e,Y_e,N_perp_re,eps,D,N_perp_im)
c----------------------------------------------------------------
c     calculates N_perp_im=Im(N_perp) by the projection method
c     for the relativistic dispersion function D=determinat.
c     N_perp_im=Im(D)/(dD/dRe(N_operp))
c     It uses Abhay-Ram (if id=14) or Nelson-Melby (if id=11) complex relativistic tensor
c-------------------------------------------------------------
      implicit none

c-----input
      real*8 
     &T_e,        !electron temperature in KeV
     &N_parallel, !parallel refractive index
     &X_e,         !(omega_pe/omega)**2
     &Y_e,         !(omega_ce/omega),
     &N_perp_re   !Real part of the perpendicular refractive index

c-----output
      real*8 N_perp_im  !Imaginary part of the perpendicular refractive index
      complex*16 eps    !relativistic dielectric tensor

c-----local
      complex*16 D,D_p,D_m      !dispersion function
      real*8 step,N_perp_re_p,N_perp_re_m,dD_dNperp_re
      
      step=1.d-7
   
      if(N_perp_re.gt.0.d0) then
         N_perp_re_p=N_perp_re*(1.d0+step)
         N_perp_re_m=N_perp_re*(1.d0-step)
      else
         N_perp_re_m=N_perp_re
         N_perp_re_p=N_perp_re*(1.d0+step)
      endif

      call Disp_combined(T_e,N_parallel,X_e,Y_e,N_perp_re_p,eps,D_p)
      call Disp_combined(T_e,N_parallel,X_e,Y_e,N_perp_re_m,eps,D_m)
      dD_dNperp_re= DREAL((D_p-D_m)/(N_perp_re_p-N_perp_re_m))

      call Disp_combined(T_e,N_parallel,X_e,Y_e,N_perp_re,eps,D)
      N_perp_im=dabs(dimag(D)/dD_dNperp_re)
!      write(*,*)'in relativist_absorp_projection_method_det'
!      write(*,*)'T_e,N_parallel,X_e,Y_e,N_perp_re',
!     &T_e,N_parallel,X_e,Y_e,N_perp_re
!      write(*,*)'D,dimag(D),dD_dNperp_re,N_perp_im',
!     & D,dimag(D),dD_dNperp_re,N_perp_im

      return
      end
