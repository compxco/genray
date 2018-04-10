C-ENM 9 June 2004: Put together Bessel function routines needed
C for Isaac Weiss implementation of relativistic plasma wave solver
C  ZBESJ, ZBESY, DBESK0, DBESK1, four bessel function routines
C
*DECK ZBESJ
      SUBROUTINE ZBESJ (ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
C***BEGIN PROLOGUE  ZBESJ
C***PURPOSE  Compute a sequence of the Bessel functions J(a,z) for
C            complex argument z and real nonnegative orders a=b,b+1,
C            b+2,... where b>0.  A scaling option is available to
C            help avoid overflow.
C***LIBRARY   SLATEC
C***CATEGORY  C10A4
C***TYPE      COMPLEX (CBESJ-C, ZBESJ-C)
C***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
C             BESSEL FUNCTIONS OF THE FIRST KIND, J BESSEL FUNCTIONS
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         On KODE=1, ZBESJ computes an N member sequence of complex
C         Bessel functions CY(L)=J(FNU+L-1,Z) for real nonnegative
C         orders FNU+L-1, L=1,...,N and complex Z in the cut plane
C         -pi<arg(Z)<=pi where Z=ZR+i*ZI.  On KODE=2, CBESJ returns
C         the scaled functions
C
C            CY(L) = exp(-abs(Y))*J(FNU+L-1,Z),  L=1,...,N and Y=Im(Z)
C
C         which remove the exponential growth in both the upper and
C         lower half planes as Z goes to infinity.  Definitions and
C         notation are found in the NBS Handbook of Mathematical
C         Functions (Ref. 1).
C
C         Input
C           ZR     - DOUBLE PRECISION real part of argument Z
C           ZI     - DOUBLE PRECISION imag part of argument Z
C           FNU    - DOUBLE PRECISION initial order, FNU>=0
C           KODE   - A parameter to indicate the scaling option
C                    KODE=1  returns
C                            CY(L)=J(FNU+L-1,Z), L=1,...,N
C                        =2  returns
C                            CY(L)=J(FNU+L-1,Z)*exp(-abs(Y)), L=1,...,N
C                            where Y=Im(Z)
C           N      - Number of terms in the sequence, N>=1
C
C         Output
C           CYR    - DOUBLE PRECISION real part of result vector
C           CYI    - DOUBLE PRECISION imag part of result vector
C           NZ     - Number of underflows set to zero
C                    NZ=0    Normal return
C                    NZ>0    CY(L)=0, L=N-NZ+1,...,N
C           IERR   - Error flag
C                    IERR=0  Normal return     - COMPUTATION COMPLETED
C                    IERR=1  Input error       - NO COMPUTATION
C                    IERR=2  Overflow          - NO COMPUTATION
C                            (Im(Z) too large on KODE=1)
C                    IERR=3  Precision warning - COMPUTATION COMPLETED
C                            (Result has half precision or less
C                            because abs(Z) or FNU+N-1 is large)
C                    IERR=4  Precision error   - NO COMPUTATION
C                            (Result has no precision because
C                            abs(Z) or FNU+N-1 is too large)
C                    IERR=5  Algorithmic error - NO COMPUTATION
C                            (Termination condition not met)
C
C *Long Description:
C
C         The computation is carried out by the formulae
C
C            J(a,z) = exp( a*pi*i/2)*I(a,-i*z),  Im(z)>=0
C
C            J(a,z) = exp(-a*pi*i/2)*I(a, i*z),  Im(z)<0
C
C         where the I Bessel function is computed as described in the
C         prologue to CBESI.
C
C         For negative orders, the formula
C
C            J(-a,z) = J(a,z)*cos(a*pi) - Y(a,z)*sin(a*pi)
C
C         can be used.  However, for large orders close to integers, the
C         the function changes radically.  When a is a large positive
C         integer, the magnitude of J(-a,z)=J(a,z)*cos(a*pi) is a
C         large negative power of ten.  But when a is not an integer,
C         Y(a,z) dominates in magnitude with a large positive power of
C         ten and the most that the second term can be reduced is by
C         unit roundoff from the coefficient.  Thus, wide changes can
C         occur within unit roundoff of a large integer for a.  Here,
C         large means a>abs(z).
C
C         In most complex variable computation, one must evaluate ele-
C         mentary functions.  When the magnitude of Z or FNU+N-1 is
C         large, losses of significance by argument reduction occur.
C         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
C         losses exceeding half precision are likely and an error flag
C         IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is double
C         precision unit roundoff limited to 18 digits precision.  Also,
C         if either is larger than U2=0.5/UR, then all significance is
C         lost and IERR=4.  In order to use the INT function, arguments
C         must be further restricted not to exceed the largest machine
C         integer, U3=I1MACH(9).  Thus, the magnitude of Z and FNU+N-1
C         is restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, and
C         U3 approximate 2.0E+3, 4.2E+6, 2.1E+9 in single precision
C         and 4.7E+7, 2.3E+15 and 2.1E+9 in double precision.  This
C         makes U2 limiting in single precision and U3 limiting in
C         double precision.  This means that one can expect to retain,
C         in the worst cases on IEEE machines, no digits in single pre-
C         cision and only 6 digits in double precision.  Similar con-
C         siderations hold for other machines.
C
C         The approximate relative error in the magnitude of a complex
C         Bessel function can be expressed as P*10**S where P=MAX(UNIT
C         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre-
C         sents the increase in error due to argument reduction in the
C         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))),
C         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF
C         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may
C         have only absolute accuracy.  This is most likely to occur
C         when one component (in magnitude) is larger than the other by
C         several orders of magnitude.  If one component is 10**K larger
C         than the other, then one can expect only MAX(ABS(LOG10(P))-K,
C         0) significant digits; or, stated another way, when K exceeds
C         the exponent of P, no significant digits remain in the smaller
C         component.  However, the phase angle retains absolute accuracy
C         because, in complex arithmetic with precision P, the smaller
C         component will not (as a rule) decrease below P times the
C         magnitude of the larger component.  In these extreme cases,
C         the principal phase angle is on the order of +P, -P, PI/2-P,
C         or -PI/2+P.
C
C***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
C                 matical Functions, National Bureau of Standards
C                 Applied Mathematics Series 55, U. S. Department
C                 of Commerce, Tenth Printing (1972) or later.
C               2. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument, Report SAND83-0086, Sandia National
C                 Laboratories, Albuquerque, NM, May 1983.
C               3. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument and Large Order, Report SAND83-0643,
C                 Sandia National Laboratories, Albuquerque, NM, May
C                 1983.
C               4. D. E. Amos, A Subroutine Package for Bessel Functions
C                 of a Complex Argument and Nonnegative Order, Report
C                 SAND85-1018, Sandia National Laboratory, Albuquerque,
C                 NM, May 1985.
C               5. D. E. Amos, A portable package for Bessel functions
C                 of a complex argument and nonnegative order, ACM
C                 Transactions on Mathematical Software, 12 (September
C                 1986), pp. 265-273.
C
C***ROUTINES CALLED  D1MACH, I1MACH, ZABS, ZBINU
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   890801  REVISION DATE from Version 3.2
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   920128  Category corrected.  (WRB)
C   920811  Prologue revised.  (DWL)
C***END PROLOGUE  ZBESJ
C
C     COMPLEX CI,CSGN,CY,Z,ZN
      DOUBLE PRECISION AA, ALIM, ARG, CII, CSGNI, CSGNR, CYI, CYR, DIG,
     * ELIM, FNU, FNUL, HPI, RL, R1M5, STR, TOL, ZI, ZNI, ZNR, ZR,
     * D1MACH, BB, FN, AZ, ZABS, ASCLE, RTOL, ATOL, STI
      INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, N, NL, NZ, I1MACH
      DIMENSION CYR(N), CYI(N)
      EXTERNAL ZABS
      DATA HPI /1.57079632679489662D0/
C
C***FIRST EXECUTABLE STATEMENT  ZBESJ
      IERR = 0
      NZ=0
      IF (FNU.LT.0.0D0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
      TOL = MAX(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN(ABS(K1),ABS(K2))
      ELIM = 2.303D0*(K*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*K1
      DIG = MIN(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + MAX(-AA,-41.45D0)
      RL = 1.2D0*DIG + 3.0D0
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
C-----------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
      AZ = ZABS(ZR,ZI)
      FN = FNU+(N-1)
      AA = 0.5D0/TOL
      BB = I1MACH(9)*0.5D0
      AA = MIN(AA,BB)
      IF (AZ.GT.AA) GO TO 260
      IF (FN.GT.AA) GO TO 260
      AA = SQRT(AA)
      IF (AZ.GT.AA) IERR=3
      IF (FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     CALCULATE CSGN=EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      CII = 1.0D0
      INU = FNU
      INUH = INU/2
      IR = INU - 2*INUH
      ARG = (FNU-(INU-IR))*HPI
      CSGNR = COS(ARG)
      CSGNI = SIN(ARG)
      IF (MOD(INUH,2).EQ.0) GO TO 40
      CSGNR = -CSGNR
      CSGNI = -CSGNI
   40 CONTINUE
C-----------------------------------------------------------------------
C     ZN IS IN THE RIGHT HALF PLANE
C-----------------------------------------------------------------------
      ZNR = ZI
      ZNI = -ZR
      IF (ZI.GE.0.0D0) GO TO 50
      ZNR = -ZNR
      ZNI = -ZNI
      CSGNI = -CSGNI
      CII = -CII
   50 CONTINUE

cSAP091104
c      write(*,*)'in ZBESJ before ZBINU(ZNR, ZNI, FNU',ZNR, ZNI, FNU

      CALL ZBINU(ZNR, ZNI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL,
     * ELIM, ALIM)
      IF (NZ.LT.0) GO TO 130
      NL = N - NZ
      IF (NL.EQ.0) RETURN
      RTOL = 1.0D0/TOL
      ASCLE = D1MACH(1)*RTOL*1.0D+3
      DO 60 I=1,NL
C       STR = CYR(I)*CSGNR - CYI(I)*CSGNI
C       CYI(I) = CYR(I)*CSGNI + CYI(I)*CSGNR
C       CYR(I) = STR
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0D0
        IF (MAX(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 55
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   55   CONTINUE
        STR = AA*CSGNR - BB*CSGNI
        STI = AA*CSGNI + BB*CSGNR
        CYR(I) = STR*ATOL
        CYI(I) = STI*ATOL
        STR = -CSGNI*CII
        CSGNI = CSGNR*CII
        CSGNR = STR
   60 CONTINUE
      RETURN
  130 CONTINUE
      IF(NZ.EQ.(-2)) GO TO 140
      NZ = 0
      IERR = 2
      RETURN
  140 CONTINUE
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
C-ENM   added ZBESY here
*DECK ZBESY
      SUBROUTINE ZBESY (ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, CWRKR,
     +   CWRKI, IERR)
C***BEGIN PROLOGUE  ZBESY
C***PURPOSE  Compute a sequence of the Bessel functions Y(a,z) for
C            complex argument z and real nonnegative orders a=b,b+1,
C            b+2,... where b>0.  A scaling option is available to
C            help avoid overflow.
C***LIBRARY   SLATEC
C***CATEGORY  C10A4
C***TYPE      COMPLEX (CBESY-C, ZBESY-C)
C***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
C             BESSEL FUNCTIONS OF SECOND KIND, WEBER'S FUNCTION,
C             Y BESSEL FUNCTIONS
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         On KODE=1, ZBESY computes an N member sequence of complex
C         Bessel functions CY(L)=Y(FNU+L-1,Z) for real nonnegative
C         orders FNU+L-1, L=1,...,N and complex Z in the cut plane
C         -pi<arg(Z)<=pi where Z=ZR+i*ZI.  On KODE=2, CBESY returns
C         the scaled functions
C
C            CY(L) = exp(-abs(Y))*Y(FNU+L-1,Z),  L=1,...,N, Y=Im(Z)
C
C         which remove the exponential growth in both the upper and
C         lower half planes as Z goes to infinity.  Definitions and
C         notation are found in the NBS Handbook of Mathematical
C         Functions (Ref. 1).
C
C         Input
C           ZR     - DOUBLE PRECISION real part of nonzero argument Z
C           ZI     - DOUBLE PRECISION imag part of nonzero argument Z
C           FNU    - DOUBLE PRECISION initial order, FNU>=0
C           KODE   - A parameter to indicate the scaling option
C                    KODE=1  returns
C                            CY(L)=Y(FNU+L-1,Z), L=1,...,N
C                        =2  returns
C                            CY(L)=Y(FNU+L-1,Z)*exp(-abs(Y)), L=1,...,N
C                            where Y=Im(Z)
C           N      - Number of terms in the sequence, N>=1
C           CWRKR  - DOUBLE PRECISION work vector of dimension N
C           CWRKI  - DOUBLE PRECISION work vector of dimension N
C
C         Output
C           CYR    - DOUBLE PRECISION real part of result vector
C           CYI    - DOUBLE PRECISION imag part of result vector
C           NZ     - Number of underflows set to zero
C                    NZ=0    Normal return
C                    NZ>0    CY(L)=0 for NZ values of L, usually on
C                            KODE=2 (the underflows may not be in an
C                            uninterrupted sequence)
C           IERR   - Error flag
C                    IERR=0  Normal return     - COMPUTATION COMPLETED
C                    IERR=1  Input error       - NO COMPUTATION
C                    IERR=2  Overflow          - NO COMPUTATION
C                            (abs(Z) too small and/or FNU+N-1
C                            too large)
C                    IERR=3  Precision warning - COMPUTATION COMPLETED
C                            (Result has half precision or less
C                            because abs(Z) or FNU+N-1 is large)
C                    IERR=4  Precision error   - NO COMPUTATION
C                            (Result has no precision because
C                            abs(Z) or FNU+N-1 is too large)
C                    IERR=5  Algorithmic error - NO COMPUTATION
C                            (Termination condition not met)
C
C *Long Description:
C
C         The computation is carried out by the formula
C
C            Y(a,z) = (H(1,a,z) - H(2,a,z))/(2*i)
C
C         where the Hankel functions are computed as described in CBESH.
C
C         For negative orders, the formula
C
C            Y(-a,z) = Y(a,z)*cos(a*pi) + J(a,z)*sin(a*pi)
C
C         can be used.  However, for large orders close to half odd
C         integers the function changes radically.  When a is a large
C         positive half odd integer, the magnitude of Y(-a,z)=J(a,z)*
C         sin(a*pi) is a large negative power of ten.  But when a is
C         not a half odd integer, Y(a,z) dominates in magnitude with a
C         large positive power of ten and the most that the second term
C         can be reduced is by unit roundoff from the coefficient.
C         Thus,  wide changes can occur within unit roundoff of a large
C         half odd integer.  Here, large means a>abs(z).
C
C         In most complex variable computation, one must evaluate ele-
C         mentary functions.  When the magnitude of Z or FNU+N-1 is
C         large, losses of significance by argument reduction occur.
C         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
C         losses exceeding half precision are likely and an error flag
C         IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is double
C         precision unit roundoff limited to 18 digits precision.  Also,
C         if either is larger than U2=0.5/UR, then all significance is
C         lost and IERR=4.  In order to use the INT function, arguments
C         must be further restricted not to exceed the largest machine
C         integer, U3=I1MACH(9).  Thus, the magnitude of Z and FNU+N-1
C         is restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, and
C         U3 approximate 2.0E+3, 4.2E+6, 2.1E+9 in single precision
C         and 4.7E+7, 2.3E+15 and 2.1E+9 in double precision.  This
C         makes U2 limiting in single precision and U3 limiting in
C         double precision.  This means that one can expect to retain,
C         in the worst cases on IEEE machines, no digits in single pre-
C         cision and only 6 digits in double precision.  Similar con-
C         siderations hold for other machines.
C
C         The approximate relative error in the magnitude of a complex
C         Bessel function can be expressed as P*10**S where P=MAX(UNIT
C         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre-
C         sents the increase in error due to argument reduction in the
C         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))),
C         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF
C         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may
C         have only absolute accuracy.  This is most likely to occur
C         when one component (in magnitude) is larger than the other by
C         several orders of magnitude.  If one component is 10**K larger
C         than the other, then one can expect only MAX(ABS(LOG10(P))-K,
C         0) significant digits; or, stated another way, when K exceeds
C         the exponent of P, no significant digits remain in the smaller
C         component.  However, the phase angle retains absolute accuracy
C         because, in complex arithmetic with precision P, the smaller
C         component will not (as a rule) decrease below P times the
C         magnitude of the larger component.  In these extreme cases,
C         the principal phase angle is on the order of +P, -P, PI/2-P,
C         or -PI/2+P.
C
C***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
C                 matical Functions, National Bureau of Standards
C                 Applied Mathematics Series 55, U. S. Department
C                 of Commerce, Tenth Printing (1972) or later.
C               2. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument, Report SAND83-0086, Sandia National
C                 Laboratories, Albuquerque, NM, May 1983.
C               3. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument and Large Order, Report SAND83-0643,
C                 Sandia National Laboratories, Albuquerque, NM, May
C                 1983.
C               4. D. E. Amos, A Subroutine Package for Bessel Functions
C                 of a Complex Argument and Nonnegative Order, Report
C                 SAND85-1018, Sandia National Laboratory, Albuquerque,
C                 NM, May 1985.
C               5. D. E. Amos, A portable package for Bessel functions
C                 of a complex argument and nonnegative order, ACM
C                 Transactions on Mathematical Software, 12 (September
C                 1986), pp. 265-273.
C
C***ROUTINES CALLED  D1MACH, I1MACH, ZBESH
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   890801  REVISION DATE from Version 3.2
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   920128  Category corrected.  (WRB)
C   920811  Prologue revised.  (DWL)
C***END PROLOGUE  ZBESY
C
C     COMPLEX CWRK,CY,C1,C2,EX,HCI,Z,ZU,ZV
      DOUBLE PRECISION CWRKI, CWRKR, CYI, CYR, C1I, C1R, C2I, C2R,
     * ELIM, EXI, EXR, EY, FNU, HCII, STI, STR, TAY, ZI, ZR,
     * D1MACH, ASCLE, RTOL, ATOL, AA, BB, TOL, R1M5
      INTEGER I, IERR, K, KODE, K1, K2, N, NZ, NZ1, NZ2, I1MACH
      DIMENSION CYR(N), CYI(N), CWRKR(N), CWRKI(N)
C***FIRST EXECUTABLE STATEMENT  ZBESY
      IERR = 0
      NZ=0
      IF (ZR.EQ.0.0D0 .AND. ZI.EQ.0.0D0) IERR=1
      IF (FNU.LT.0.0D0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      HCII = 0.5D0
      CALL ZBESH(ZR, ZI, FNU, KODE, 1, N, CYR, CYI, NZ1, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 170
      CALL ZBESH(ZR, ZI, FNU, KODE, 2, N, CWRKR, CWRKI, NZ2, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 170
      NZ = MIN(NZ1,NZ2)
      IF (KODE.EQ.2) GO TO 60
      DO 50 I=1,N
        STR = CWRKR(I) - CYR(I)
        STI = CWRKI(I) - CYI(I)
        CYR(I) = -STI*HCII
        CYI(I) = STR*HCII
   50 CONTINUE
      RETURN
   60 CONTINUE
      TOL = MAX(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      K = MIN(ABS(K1),ABS(K2))
      R1M5 = D1MACH(5)
C-----------------------------------------------------------------------
C     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT
C-----------------------------------------------------------------------
      ELIM = 2.303D0*(K*R1M5-3.0D0)
      EXR = COS(ZR)
      EXI = SIN(ZR)
      EY = 0.0D0
      TAY = ABS(ZI+ZI)
      IF (TAY.LT.ELIM) EY = EXP(-TAY)
      IF (ZI.LT.0.0D0) GO TO 90
      C1R = EXR*EY
      C1I = EXI*EY
      C2R = EXR
      C2I = -EXI
   70 CONTINUE
      NZ = 0
      RTOL = 1.0D0/TOL
      ASCLE = D1MACH(1)*RTOL*1.0D+3
      DO 80 I=1,N
C       STR = C1R*CYR(I) - C1I*CYI(I)
C       STI = C1R*CYI(I) + C1I*CYR(I)
C       STR = -STR + C2R*CWRKR(I) - C2I*CWRKI(I)
C       STI = -STI + C2R*CWRKI(I) + C2I*CWRKR(I)
C       CYR(I) = -STI*HCII
C       CYI(I) = STR*HCII
        AA = CWRKR(I)
        BB = CWRKI(I)
        ATOL = 1.0D0
        IF (MAX(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 75
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   75   CONTINUE
        STR = (AA*C2R - BB*C2I)*ATOL
        STI = (AA*C2I + BB*C2R)*ATOL
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0D0
        IF (MAX(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 85
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   85   CONTINUE
        STR = STR - (AA*C1R - BB*C1I)*ATOL
        STI = STI - (AA*C1I + BB*C1R)*ATOL
        CYR(I) = -STI*HCII
        CYI(I) =  STR*HCII
        IF (STR.EQ.0.0D0 .AND. STI.EQ.0.0D0 .AND. EY.EQ.0.0D0) NZ = NZ
     *   + 1
   80 CONTINUE
      RETURN
   90 CONTINUE
      C1R = EXR
      C1I = EXI
      C2R = EXR*EY
      C2I = -EXI*EY
      GO TO 70
  170 CONTINUE
      NZ = 0
      RETURN
      END

       DOUBLE PRECISION FUNCTION ZABS (ZR, ZI)
C***BEGIN PROLOGUE  ZABS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZABS-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE
C     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI)
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZABS
      DOUBLE PRECISION ZR, ZI, U, V, Q, S
C***FIRST EXECUTABLE STATEMENT  ZABS
      U = ABS(ZR)
      V = ABS(ZI)
      S = U + V
C-----------------------------------------------------------------------
C     S*1.0D0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A
C     TRUE FLOATING ZERO
C-----------------------------------------------------------------------
      S = S*1.0D+0
      IF (S.EQ.0.0D+0) GO TO 20
      IF (U.GT.V) GO TO 10
      Q = U/V
      ZABS = V*SQRT(1.D+0+Q*Q)
      RETURN
   10 Q = V/U
      ZABS = U*SQRT(1.D+0+Q*Q)
      RETURN
   20 ZABS = 0.0D+0
      RETURN
      END
   
*DECK ZBINU
      SUBROUTINE ZBINU (ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL,
     +   TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZBINU
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK and ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (CBINU-A, ZBINU-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBIRY
C***ROUTINES CALLED  ZABS, ZASYI, ZBUNI, ZMLRI, ZSERI, ZUOIK, ZWRSK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZBINU
      DOUBLE PRECISION ALIM, AZ, CWI, CWR, CYI, CYR, DFNU, ELIM, FNU,
     * FNUL, RL, TOL, ZEROI, ZEROR, ZI, ZR, ZABS
      INTEGER I, INW, KODE, N, NLAST, NN, NUI, NW, NZ
      DIMENSION CYR(N), CYI(N), CWR(2), CWI(2)
      EXTERNAL ZABS
      DATA ZEROR,ZEROI / 0.0D0, 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZBINU
      NZ = 0
      AZ = ZABS(ZR,ZI)
      NN = N
      DFNU = FNU + (N-1)
      IF (AZ.LE.2.0D0) GO TO 10
      IF (AZ*AZ*0.25D0.GT.DFNU+1.0D0) GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES
C-----------------------------------------------------------------------
cSAP091104
c      write(*,*)'in UBROUTINE ZBINU ZR, ZI, FNU',ZR, ZI, FNU

      CALL ZSERI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
      INW = ABS(NW)
      NZ = NZ + INW
      NN = NN - INW
      IF (NN.EQ.0) RETURN
      IF (NW.GE.0) GO TO 120
      DFNU = FNU + (NN-1)
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 40
      IF (DFNU.LE.1.0D0) GO TO 30
      IF (AZ+AZ.LT.DFNU*DFNU) GO TO 50
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z
C-----------------------------------------------------------------------
   30 CONTINUE
      CALL ZASYI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, RL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
   40 CONTINUE
      IF (DFNU.LE.1.0D0) GO TO 70
   50 CONTINUE
C-----------------------------------------------------------------------
C     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, NN, CYR, CYI, NW, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      NN = NN - NW
      IF (NN.EQ.0) RETURN
      DFNU = FNU+(NN-1)
      IF (DFNU.GT.FNUL) GO TO 110
      IF (AZ.GT.FNUL) GO TO 110
   60 CONTINUE
      IF (AZ.GT.RL) GO TO 80
   70 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES
C-----------------------------------------------------------------------
      CALL ZMLRI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL)
      IF(NW.LT.0) GO TO 130
      GO TO 120
   80 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
C-----------------------------------------------------------------------
      CALL ZUOIK(ZR, ZI, FNU, KODE, 2, 2, CWR, CWI, NW, TOL, ELIM,
     * ALIM)
      IF (NW.GE.0) GO TO 100
      NZ = NN
      DO 90 I=1,NN
        CYR(I) = ZEROR
        CYI(I) = ZEROI
   90 CONTINUE
      RETURN
  100 CONTINUE
      IF (NW.GT.0) GO TO 130
      CALL ZWRSK(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, CWR, CWI, TOL,
     * ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
  110 CONTINUE
C-----------------------------------------------------------------------
C     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
C-----------------------------------------------------------------------
      NUI = FNUL-DFNU + 1
      NUI = MAX(NUI,0)
      CALL ZBUNI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, NUI, NLAST, FNUL,
     * TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      IF (NLAST.EQ.0) GO TO 120
      NN = NLAST
      GO TO 60
  120 CONTINUE
      RETURN
  130 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
C-ENM added ZBESH here
*DECK ZBESH
      SUBROUTINE ZBESH (ZR, ZI, FNU, KODE, M, N, CYR, CYI, NZ, IERR)
C***BEGIN PROLOGUE  ZBESH
C***PURPOSE  Compute a sequence of the Hankel functions H(m,a,z)
C            for superscript m=1 or 2, real nonnegative orders a=b,
C            b+1,... where b>0, and nonzero complex argument z.  A
C            scaling option is available to help avoid overflow.
C***LIBRARY   SLATEC
C***CATEGORY  C10A4
C***TYPE      COMPLEX (CBESH-C, ZBESH-C)
C***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
C             BESSEL FUNCTIONS OF THE THIRD KIND, H BESSEL FUNCTIONS,
C             HANKEL FUNCTIONS
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         On KODE=1, ZBESH computes an N member sequence of complex
C         Hankel (Bessel) functions CY(L)=H(M,FNU+L-1,Z) for super-
C         script M=1 or 2, real nonnegative orders FNU+L-1, L=1,...,
C         N, and complex nonzero Z in the cut plane -pi<arg(Z)<=pi
C         where Z=ZR+i*ZI.  On KODE=2, CBESH returns the scaled
C         functions
C
C            CY(L) = H(M,FNU+L-1,Z)*exp(-(3-2*M)*Z*i),  i**2=-1
C
C         which removes the exponential behavior in both the upper
C         and lower half planes.  Definitions and notation are found
C         in the NBS Handbook of Mathematical Functions (Ref. 1).
C
C         Input
C           ZR     - DOUBLE PRECISION real part of nonzero argument Z
C           ZI     - DOUBLE PRECISION imag part of nonzero argument Z
C           FNU    - DOUBLE PRECISION initial order, FNU>=0
C           KODE   - A parameter to indicate the scaling option
C                    KODE=1  returns
C                            CY(L)=H(M,FNU+L-1,Z), L=1,...,N
C                        =2  returns
C                            CY(L)=H(M,FNU+L-1,Z)*exp(-(3-2M)*Z*i),
C                            L=1,...,N
C           M      - Superscript of Hankel function, M=1 or 2
C           N      - Number of terms in the sequence, N>=1
C
C         Output
C           CYR    - DOUBLE PRECISION real part of result vector
C           CYI    - DOUBLE PRECISION imag part of result vector
C           NZ     - Number of underflows set to zero
C                    NZ=0    Normal return
C                    NZ>0    CY(L)=0 for NZ values of L (if M=1 and
C                            Im(Z)>0 or if M=2 and Im(Z)<0, then
C                            CY(L)=0 for L=1,...,NZ; in the com-
C                            plementary half planes, the underflows
C                            may not be in an uninterrupted sequence)
C           IERR   - Error flag
C                    IERR=0  Normal return     - COMPUTATION COMPLETED
C                    IERR=1  Input error       - NO COMPUTATION
C                    IERR=2  Overflow          - NO COMPUTATION
C                            (abs(Z) too small and/or FNU+N-1
C                            too large)
C                    IERR=3  Precision warning - COMPUTATION COMPLETED
C                            (Result has half precision or less
C                            because abs(Z) or FNU+N-1 is large)
C                    IERR=4  Precision error   - NO COMPUTATION
C                            (Result has no precision because
C                            abs(Z) or FNU+N-1 is too large)
C                    IERR=5  Algorithmic error - NO COMPUTATION
C                            (Termination condition not met)
C
C *Long Description:
C
C         The computation is carried out by the formula
C
C            H(m,a,z) = (1/t)*exp(-a*t)*K(a,z*exp(-t))
C                   t = (3-2*m)*i*pi/2
C
C         where the K Bessel function is computed as described in the
C         prologue to CBESK.
C
C         Exponential decay of H(m,a,z) occurs in the upper half z
C         plane for m=1 and the lower half z plane for m=2.  Exponential
C         growth occurs in the complementary half planes.  Scaling
C         by exp(-(3-2*m)*z*i) removes the exponential behavior in the
C         whole z plane as z goes to infinity.
C
C         For negative orders, the formula
C
C            H(m,-a,z) = H(m,a,z)*exp((3-2*m)*a*pi*i)
C
C         can be used.
C
C         In most complex variable computation, one must evaluate ele-
C         mentary functions.  When the magnitude of Z or FNU+N-1 is
C         large, losses of significance by argument reduction occur.
C         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
C         losses exceeding half precision are likely and an error flag
C         IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is double
C         precision unit roundoff limited to 18 digits precision.  Also,
C         if either is larger than U2=0.5/UR, then all significance is
C         lost and IERR=4.  In order to use the INT function, arguments
C         must be further restricted not to exceed the largest machine
C         integer, U3=I1MACH(9).  Thus, the magnitude of Z and FNU+N-1
C         is restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, and
C         U3 approximate 2.0E+3, 4.2E+6, 2.1E+9 in single precision
C         and 4.7E+7, 2.3E+15 and 2.1E+9 in double precision.  This
C         makes U2 limiting in single precision and U3 limiting in
C         double precision.  This means that one can expect to retain,
C         in the worst cases on IEEE machines, no digits in single pre-
C         cision and only 6 digits in double precision.  Similar con-
C         siderations hold for other machines.
C
C         The approximate relative error in the magnitude of a complex
C         Bessel function can be expressed as P*10**S where P=MAX(UNIT
C         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre-
C         sents the increase in error due to argument reduction in the
C         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))),
C         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF
C         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may
C         have only absolute accuracy.  This is most likely to occur
C         when one component (in magnitude) is larger than the other by
C         several orders of magnitude.  If one component is 10**K larger
C         than the other, then one can expect only MAX(ABS(LOG10(P))-K,
C         0) significant digits; or, stated another way, when K exceeds
C         the exponent of P, no significant digits remain in the smaller
C         component.  However, the phase angle retains absolute accuracy
C         because, in complex arithmetic with precision P, the smaller
C         component will not (as a rule) decrease below P times the
C         magnitude of the larger component.  In these extreme cases,
C         the principal phase angle is on the order of +P, -P, PI/2-P,
C         or -PI/2+P.
C
C***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
C                 matical Functions, National Bureau of Standards
C                 Applied Mathematics Series 55, U. S. Department
C                 of Commerce, Tenth Printing (1972) or later.
C               2. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument, Report SAND83-0086, Sandia National
C                 Laboratories, Albuquerque, NM, May 1983.
C               3. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument and Large Order, Report SAND83-0643,
C                 Sandia National Laboratories, Albuquerque, NM, May
C                 1983.
C               4. D. E. Amos, A Subroutine Package for Bessel Functions
C                 of a Complex Argument and Nonnegative Order, Report
C                 SAND85-1018, Sandia National Laboratory, Albuquerque,
C                 NM, May 1985.
C               5. D. E. Amos, A portable package for Bessel functions
C                 of a complex argument and nonnegative order, ACM
C                 Transactions on Mathematical Software, 12 (September
C                 1986), pp. 265-273.
C
C***ROUTINES CALLED  D1MACH, I1MACH, ZABS, ZACON, ZBKNU, ZBUNK, ZUOIK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   890801  REVISION DATE from Version 3.2
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   920128  Category corrected.  (WRB)
C   920811  Prologue revised.  (DWL)
C***END PROLOGUE  ZBESH
C
C     COMPLEX CY,Z,ZN,ZT,CSGN
      DOUBLE PRECISION AA, ALIM, ALN, ARG, AZ, CYI, CYR, DIG, ELIM,
     * FMM, FN, FNU, FNUL, HPI, RHPI, RL, R1M5, SGN, STR, TOL, UFL, ZI,
     * ZNI, ZNR, ZR, ZTI, D1MACH, ZABS, BB, ASCLE, RTOL, ATOL, STI,
     * CSGNR, CSGNI
      INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, M,
     * MM, MR, N, NN, NUF, NW, NZ, I1MACH
      DIMENSION CYR(N), CYI(N)
      EXTERNAL ZABS
C
      DATA HPI /1.57079632679489662D0/
C
C***FIRST EXECUTABLE STATEMENT  ZBESH
      IERR = 0
      NZ=0
      IF (ZR.EQ.0.0D0 .AND. ZI.EQ.0.0D0) IERR=1
      IF (FNU.LT.0.0D0) IERR=1
      IF (M.LT.1 .OR. M.GT.2) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      NN = N
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      TOL = MAX(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN(ABS(K1),ABS(K2))
      ELIM = 2.303D0*(K*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*K1
      DIG = MIN(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + MAX(-AA,-41.45D0)
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
      RL = 1.2D0*DIG + 3.0D0
      FN = FNU + (NN-1)
      MM = 3 - M - M
      FMM = MM
      ZNR = FMM*ZI
      ZNI = -FMM*ZR
C-----------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
      AZ = ZABS(ZR,ZI)
      AA = 0.5D0/TOL
      BB = I1MACH(9)*0.5D0
      AA = MIN(AA,BB)
      IF (AZ.GT.AA) GO TO 260
      IF (FN.GT.AA) GO TO 260
      AA = SQRT(AA)
      IF (AZ.GT.AA) IERR=3
      IF (FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
C-----------------------------------------------------------------------
      UFL = D1MACH(1)*1.0D+3
      IF (AZ.LT.UFL) GO TO 230
      IF (FNU.GT.FNUL) GO TO 90
      IF (FN.LE.1.0D0) GO TO 70
      IF (FN.GT.2.0D0) GO TO 60
      IF (AZ.GT.TOL) GO TO 70
      ARG = 0.5D0*AZ
      ALN = -FN*LOG(ARG)
      IF (ALN.GT.ELIM) GO TO 230
      GO TO 70
   60 CONTINUE
      CALL ZUOIK(ZNR, ZNI, FNU, KODE, 2, NN, CYR, CYI, NUF, TOL, ELIM,
     * ALIM)
      IF (NUF.LT.0) GO TO 230
      NZ = NZ + NUF
      NN = NN - NUF
C-----------------------------------------------------------------------
C     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
C     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C-----------------------------------------------------------------------
      IF (NN.EQ.0) GO TO 140
   70 CONTINUE
      IF ((ZNR.LT.0.0D0) .OR. (ZNR.EQ.0.0D0 .AND. ZNI.LT.0.0D0 .AND.
     * M.EQ.2)) GO TO 80
C-----------------------------------------------------------------------
C     RIGHT HALF PLANE COMPUTATION, XN.GE.0. .AND. (XN.NE.0. .OR.
C     YN.GE.0. .OR. M=1)
C-----------------------------------------------------------------------
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, NN, CYR, CYI, NZ, TOL, ELIM, ALIM)
      GO TO 110
C-----------------------------------------------------------------------
C     LEFT HALF PLANE COMPUTATION
C-----------------------------------------------------------------------
   80 CONTINUE
      MR = -MM
      CALL ZACON(ZNR, ZNI, FNU, KODE, MR, NN, CYR, CYI, NW, RL, FNUL,
     * TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 240
      NZ=NW
      GO TO 110
   90 CONTINUE
C-----------------------------------------------------------------------
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
C-----------------------------------------------------------------------
      MR = 0
      IF ((ZNR.GE.0.0D0) .AND. (ZNR.NE.0.0D0 .OR. ZNI.GE.0.0D0 .OR.
     * M.NE.2)) GO TO 100
      MR = -MM
      IF (ZNR.NE.0.0D0 .OR. ZNI.GE.0.0D0) GO TO 100
      ZNR = -ZNR
      ZNI = -ZNI
  100 CONTINUE
      CALL ZBUNK(ZNR, ZNI, FNU, KODE, MR, NN, CYR, CYI, NW, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 240
      NZ = NZ + NW
  110 CONTINUE
C-----------------------------------------------------------------------
C     H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)
C
C     ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
C-----------------------------------------------------------------------
      SGN = DSIGN(HPI,-FMM)
C-----------------------------------------------------------------------
C     CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = FNU
      INUH = INU/2
      IR = INU - 2*INUH
      ARG = (FNU-(INU-IR))*SGN
      RHPI = 1.0D0/SGN
C     ZNI = RHPI*COS(ARG)
C     ZNR = -RHPI*SIN(ARG)
      CSGNI = RHPI*COS(ARG)
      CSGNR = -RHPI*SIN(ARG)
      IF (MOD(INUH,2).EQ.0) GO TO 120
C     ZNR = -ZNR
C     ZNI = -ZNI
      CSGNR = -CSGNR
      CSGNI = -CSGNI
  120 CONTINUE
      ZTI = -FMM
      RTOL = 1.0D0/TOL
      ASCLE = UFL*RTOL
      DO 130 I=1,NN
C       STR = CYR(I)*ZNR - CYI(I)*ZNI
C       CYI(I) = CYR(I)*ZNI + CYI(I)*ZNR
C       CYR(I) = STR
C       STR = -ZNI*ZTI
C       ZNI = ZNR*ZTI
C       ZNR = STR
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0D0
        IF (MAX(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 135
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
  135 CONTINUE
      STR = AA*CSGNR - BB*CSGNI
      STI = AA*CSGNI + BB*CSGNR
      CYR(I) = STR*ATOL
      CYI(I) = STI*ATOL
      STR = -CSGNI*ZTI
      CSGNI = CSGNR*ZTI
      CSGNR = STR
  130 CONTINUE
      RETURN
  140 CONTINUE
      IF (ZNR.LT.0.0D0) GO TO 230
      RETURN
  230 CONTINUE
      NZ=0
      IERR=2
      RETURN
  240 CONTINUE
      IF(NW.EQ.(-1)) GO TO 230
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
*DECK ZBINU



      SUBROUTINE ZSERI (ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  ZSERI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CSERI-A, ZSERI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE POWER SERIES FOR LARGE ABS(Z) IN THE
C     REGION ABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
C     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
C     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE
C     CONDITION ABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE
C     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, DGAMLN, ZABS, ZDIV, ZLOG, ZMLT, ZUCHK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZLOG to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZSERI
C     COMPLEX AK1,CK,COEF,CONE,CRSC,CSCL,CZ,CZERO,HZ,RZ,S1,S2,Y,Z
      DOUBLE PRECISION AA, ACZ, AK, AK1I, AK1R, ALIM, ARM, ASCLE, ATOL,
     * AZ, CKI, CKR, COEFI, COEFR, CONEI, CONER, CRSCR, CZI, CZR, DFNU,
     * ELIM, FNU, FNUP, HZI, HZR, RAZ, RS, RTR1, RZI, RZR, S, SS, STI,
     * STR, S1I, S1R, S2I, S2R, TOL, YI, YR, WI, WR, ZEROI, ZEROR, ZI,
     * ZR, DGAMLN, D1MACH, ZABS
      INTEGER I, IB, IDUM, IFLAG, IL, K, KODE, L, M, N, NN, NZ, NW
      DIMENSION YR(N), YI(N), WR(2), WI(2)
      EXTERNAL ZABS, ZLOG
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZSERI
cSAP091104
c      write(*,*)'in ZSERI ZR, ZI, FNU,N',ZR, ZI, FNU,N

      NZ = 0
      AZ = ZABS(ZR,ZI)
      IF (AZ.EQ.0.0D0) GO TO 160
      ARM = 1.0D+3*D1MACH(1)
      RTR1 = SQRT(ARM)
      CRSCR = 1.0D0
      IFLAG = 0
      IF (AZ.LT.ARM) GO TO 150
      HZR = 0.5D0*ZR
      HZI = 0.5D0*ZI
      CZR = ZEROR
      CZI = ZEROI
      IF (AZ.LE.RTR1) GO TO 10
      CALL ZMLT(HZR, HZI, HZR, HZI, CZR, CZI)
   10 CONTINUE
      ACZ = ZABS(CZR,CZI)
      NN = N
      CALL ZLOG(HZR, HZI, CKR, CKI, IDUM)
   20 CONTINUE
      DFNU = FNU + (NN-1)

cSAP091104
c      write(*,*)'in ZSERI after 20  FNU,NN,DFNU', FNU,NN,DFNU

      FNUP = DFNU + 1.0D0

cSAP091104
c      write(*,*)'in ZSERI after 20  FNUP', FNUP

C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      AK1R = CKR*DFNU
      AK1I = CKI*DFNU

cSAP091104
c      write(*,*) 'in ZSERI FNUP,IDUM',FNUP,IDUM

      AK = DGAMLN(FNUP,IDUM)
      AK1R = AK1R - AK
      IF (KODE.EQ.2) AK1R = AK1R - ZR
      IF (AK1R.GT.(-ELIM)) GO TO 40
   30 CONTINUE
      NZ = NZ + 1
      YR(NN) = ZEROR
      YI(NN) = ZEROI
      IF (ACZ.GT.DFNU) GO TO 190
      NN = NN - 1
      IF (NN.EQ.0) RETURN
      GO TO 20
   40 CONTINUE
      IF (AK1R.GT.(-ALIM)) GO TO 50
      IFLAG = 1
      SS = 1.0D0/TOL
      CRSCR = TOL
      ASCLE = ARM*SS
   50 CONTINUE
      AA = EXP(AK1R)
      IF (IFLAG.EQ.1) AA = AA*SS
      COEFR = AA*COS(AK1I)
      COEFI = AA*SIN(AK1I)
      ATOL = TOL*ACZ/FNUP
      IL = MIN(2,NN)
      DO 90 I=1,IL
        DFNU = FNU + (NN-I)
        FNUP = DFNU + 1.0D0
        S1R = CONER
        S1I = CONEI
        IF (ACZ.LT.TOL*FNUP) GO TO 70
        AK1R = CONER
        AK1I = CONEI
        AK = FNUP + 2.0D0
        S = FNUP
        AA = 2.0D0
   60   CONTINUE
        RS = 1.0D0/S
        STR = AK1R*CZR - AK1I*CZI
        STI = AK1R*CZI + AK1I*CZR
        AK1R = STR*RS
        AK1I = STI*RS
        S1R = S1R + AK1R
        S1I = S1I + AK1I
        S = S + AK
        AK = AK + 2.0D0
        AA = AA*ACZ*RS
        IF (AA.GT.ATOL) GO TO 60
   70   CONTINUE
        S2R = S1R*COEFR - S1I*COEFI
        S2I = S1R*COEFI + S1I*COEFR
        WR(I) = S2R
        WI(I) = S2I
        IF (IFLAG.EQ.0) GO TO 80
        CALL ZUCHK(S2R, S2I, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 30
   80   CONTINUE
        M = NN - I + 1
        YR(M) = S2R*CRSCR
        YI(M) = S2I*CRSCR
        IF (I.EQ.IL) GO TO 90
        CALL ZDIV(COEFR, COEFI, HZR, HZI, STR, STI)
        COEFR = STR*DFNU
        COEFI = STI*DFNU
   90 CONTINUE
      IF (NN.LE.2) RETURN
      K = NN - 2
      AK = K
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      IF (IFLAG.EQ.1) GO TO 120
      IB = 3
  100 CONTINUE
      DO 110 I=IB,NN
        YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2)
        YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2)
        AK = AK - 1.0D0
        K = K - 1
  110 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD WITH SCALED VALUES
C-----------------------------------------------------------------------
  120 CONTINUE
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
C     UNDERFLOW LIMIT = ASCLE = D1MACH(1)*SS*1.0D+3
C-----------------------------------------------------------------------
      S1R = WR(1)
      S1I = WI(1)
      S2R = WR(2)
      S2I = WI(2)
      DO 130 L=3,NN
        CKR = S2R
        CKI = S2I
        S2R = S1R + (AK+FNU)*(RZR*CKR-RZI*CKI)
        S2I = S1I + (AK+FNU)*(RZR*CKI+RZI*CKR)
        S1R = CKR
        S1I = CKI
        CKR = S2R*CRSCR
        CKI = S2I*CRSCR
        YR(K) = CKR
        YI(K) = CKI
        AK = AK - 1.0D0
        K = K - 1
        IF (ZABS(CKR,CKI).GT.ASCLE) GO TO 140
  130 CONTINUE
      RETURN
  140 CONTINUE
      IB = L + 1
      IF (IB.GT.NN) RETURN
      GO TO 100
  150 CONTINUE
      NZ = N
      IF (FNU.EQ.0.0D0) NZ = NZ - 1
  160 CONTINUE
      YR(1) = ZEROR
      YI(1) = ZEROI
      IF (FNU.NE.0.0D0) GO TO 170
      YR(1) = CONER
      YI(1) = CONEI
  170 CONTINUE
      IF (N.EQ.1) RETURN
      DO 180 I=2,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  180 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     RETURN WITH NZ.LT.0 IF ABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE
C     THE CALCULATION IN CBINU WITH N=N-ABS(NZ)
C-----------------------------------------------------------------------
  190 CONTINUE
      NZ = -NZ
      RETURN
      END


*DECK ZASYI
      SUBROUTINE ZASYI (ZR, ZI, FNU, KODE, N, YR, YI, NZ, RL, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  ZASYI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CASYI-A, ZASYI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE ABS(Z) IN THE
C     REGION ABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
C     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1.
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZDIV, ZEXP, ZMLT, ZSQRT
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZEXP and ZSQRT to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZASYI
C     COMPLEX AK1,CK,CONE,CS1,CS2,CZ,CZERO,DK,EZ,P1,RZ,S2,Y,Z
      DOUBLE PRECISION AA, AEZ, AK, AK1I, AK1R, ALIM, ARG, ARM, ATOL,
     * AZ, BB, BK, CKI, CKR, CONEI, CONER, CS1I, CS1R, CS2I, CS2R, CZI,
     * CZR, DFNU, DKI, DKR, DNU2, ELIM, EZI, EZR, FDN, FNU, PI, P1I,
     * P1R, RAZ, RL, RTPI, RTR1, RZI, RZR, S, SGN, SQK, STI, STR, S2I,
     * S2R, TOL, TZI, TZR, YI, YR, ZEROI, ZEROR, ZI, ZR, D1MACH, ZABS
      INTEGER I, IB, IL, INU, J, JL, K, KODE, KODED, M, N, NN, NZ
      DIMENSION YR(N), YI(N)
      EXTERNAL ZABS, ZEXP, ZSQRT
      DATA PI, RTPI  /3.14159265358979324D0 , 0.159154943091895336D0 /
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZASYI
      NZ = 0
      AZ = ZABS(ZR,ZI)
      ARM = 1.0D+3*D1MACH(1)
      RTR1 = SQRT(ARM)
      IL = MIN(2,N)
      DFNU = FNU + (N-IL)
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      AK1R = RTPI*STR*RAZ
      AK1I = RTPI*STI*RAZ
      CALL ZSQRT(AK1R, AK1I, AK1R, AK1I)
      CZR = ZR
      CZI = ZI
      IF (KODE.NE.2) GO TO 10
      CZR = ZEROR
      CZI = ZI
   10 CONTINUE
      IF (ABS(CZR).GT.ELIM) GO TO 100
      DNU2 = DFNU + DFNU
      KODED = 1
      IF ((ABS(CZR).GT.ALIM) .AND. (N.GT.2)) GO TO 20
      KODED = 0
      CALL ZEXP(CZR, CZI, STR, STI)
      CALL ZMLT(AK1R, AK1I, STR, STI, AK1R, AK1I)
   20 CONTINUE
      FDN = 0.0D0
      IF (DNU2.GT.RTR1) FDN = DNU2*DNU2
      EZR = ZR*8.0D0
      EZI = ZI*8.0D0
C-----------------------------------------------------------------------
C     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
C     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
C     EXPANSION FOR THE IMAGINARY PART.
C-----------------------------------------------------------------------
      AEZ = 8.0D0*AZ
      S = TOL/AEZ
      JL = RL+RL + 2
      P1R = ZEROR
      P1I = ZEROI
      IF (ZI.EQ.0.0D0) GO TO 30
C-----------------------------------------------------------------------
C     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
C     SIGNIFICANCE WHEN FNU OR N IS LARGE
C-----------------------------------------------------------------------
      INU = FNU
      ARG = (FNU-INU)*PI
      INU = INU + N - IL
      AK = -SIN(ARG)
      BK = COS(ARG)
      IF (ZI.LT.0.0D0) BK = -BK
      P1R = AK
      P1I = BK
      IF (MOD(INU,2).EQ.0) GO TO 30
      P1R = -P1R
      P1I = -P1I
   30 CONTINUE
      DO 70 K=1,IL
        SQK = FDN - 1.0D0
        ATOL = S*ABS(SQK)
        SGN = 1.0D0
        CS1R = CONER
        CS1I = CONEI
        CS2R = CONER
        CS2I = CONEI
        CKR = CONER
        CKI = CONEI
        AK = 0.0D0
        AA = 1.0D0
        BB = AEZ
        DKR = EZR
        DKI = EZI
        DO 40 J=1,JL
          CALL ZDIV(CKR, CKI, DKR, DKI, STR, STI)
          CKR = STR*SQK
          CKI = STI*SQK
          CS2R = CS2R + CKR
          CS2I = CS2I + CKI
          SGN = -SGN
          CS1R = CS1R + CKR*SGN
          CS1I = CS1I + CKI*SGN
          DKR = DKR + EZR
          DKI = DKI + EZI
          AA = AA*ABS(SQK)/BB
          BB = BB + AEZ
          AK = AK + 8.0D0
          SQK = SQK - AK
          IF (AA.LE.ATOL) GO TO 50
   40   CONTINUE
        GO TO 110
   50   CONTINUE
        S2R = CS1R
        S2I = CS1I
        IF (ZR+ZR.GE.ELIM) GO TO 60
        TZR = ZR + ZR
        TZI = ZI + ZI
        CALL ZEXP(-TZR, -TZI, STR, STI)
        CALL ZMLT(STR, STI, P1R, P1I, STR, STI)
        CALL ZMLT(STR, STI, CS2R, CS2I, STR, STI)
        S2R = S2R + STR
        S2I = S2I + STI
   60   CONTINUE
        FDN = FDN + 8.0D0*DFNU + 4.0D0
        P1R = -P1R
        P1I = -P1I
        M = N - IL + K
        YR(M) = S2R*AK1R - S2I*AK1I
        YI(M) = S2R*AK1I + S2I*AK1R
   70 CONTINUE
      IF (N.LE.2) RETURN
      NN = N
      K = NN - 2
      AK = K
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      IB = 3
      DO 80 I=IB,NN
        YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2)
        YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2)
        AK = AK - 1.0D0
        K = K - 1
   80 CONTINUE
      IF (KODED.EQ.0) RETURN
      CALL ZEXP(CZR, CZI, CKR, CKI)
      DO 90 I=1,NN
        STR = YR(I)*CKR - YI(I)*CKI
        YI(I) = YR(I)*CKI + YI(I)*CKR
        YR(I) = STR
   90 CONTINUE
      RETURN
  100 CONTINUE
      NZ = -1
      RETURN
  110 CONTINUE
      NZ=-2
      RETURN
      END

*DECK ZUOIK
      SUBROUTINE ZUOIK (ZR, ZI, FNU, KODE, IKFLG, N, YR, YI, NUF, TOL,
     +   ELIM, ALIM)
C***BEGIN PROLOGUE  ZUOIK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUOIK-A, ZUOIK-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
C     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
C     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
C     WHERE ALIM.LT.ELIM. IF THE MAGNITUDE, BASED ON THE LEADING
C     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
C     THE RESULT IS ON SCALE. IF NOT, THEN A REFINED TEST USING OTHER
C     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
C     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
C     EXP(-ELIM)/TOL
C
C     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
C          =2 MEANS THE K SEQUENCE IS TESTED
C     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
C         =-1 MEANS AN OVERFLOW WOULD OCCUR
C     IKFLG=1 AND NUF.GT.0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
C             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
C     IKFLG=2 AND NUF.EQ.N MEANS ALL Y VALUES WERE SET TO ZERO
C     IKFLG=2 AND 0.LT.NUF.LT.N NOT CONSIDERED. Y MUST BE SET BY
C             ANOTHER ROUTINE
C
C***SEE ALSO  ZBESH, ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZLOG, ZUCHK, ZUNHJ, ZUNIK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZLOG to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZUOIK
C     COMPLEX ARG,ASUM,BSUM,CWRK,CZ,CZERO,PHI,SUM,Y,Z,ZB,ZETA1,ZETA2,ZN,
C    *ZR
      DOUBLE PRECISION AARG, AIC, ALIM, APHI, ARGI, ARGR, ASUMI, ASUMR,
     * ASCLE, AX, AY, BSUMI, BSUMR, CWRKI, CWRKR, CZI, CZR, ELIM, FNN,
     * FNU, GNN, GNU, PHII, PHIR, RCZ, STR, STI, SUMI, SUMR, TOL, YI,
     * YR, ZBI, ZBR, ZEROI, ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI,
     * ZNI, ZNR, ZR, ZRI, ZRR, D1MACH, ZABS
      INTEGER I, IDUM, IFORM, IKFLG, INIT, KODE, N, NN, NUF, NW
      DIMENSION YR(N), YI(N), CWRKR(16), CWRKI(16)
      EXTERNAL ZABS, ZLOG
      DATA ZEROR,ZEROI / 0.0D0, 0.0D0 /
      DATA AIC / 1.265512123484645396D+00 /
C***FIRST EXECUTABLE STATEMENT  ZUOIK
      NUF = 0
      NN = N
      ZRR = ZR
      ZRI = ZI
      IF (ZR.GE.0.0D0) GO TO 10
      ZRR = -ZR
      ZRI = -ZI
   10 CONTINUE
      ZBR = ZRR
      ZBI = ZRI
      AX = ABS(ZR)*1.7321D0
      AY = ABS(ZI)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      GNU = MAX(FNU,1.0D0)
      IF (IKFLG.EQ.1) GO TO 20
      FNN = NN
      GNN = FNU + FNN - 1.0D0
      GNU = MAX(GNN,FNN)
   20 CONTINUE
C-----------------------------------------------------------------------
C     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
C     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
C     THE SIGN OF THE IMAGINARY PART CORRECT.
C-----------------------------------------------------------------------
      IF (IFORM.EQ.2) GO TO 30
      INIT = 0
      CALL ZUNIK(ZRR, ZRI, GNU, IKFLG, 1, TOL, INIT, PHIR, PHII,
     * ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      GO TO 50
   30 CONTINUE
      ZNR = ZRI
      ZNI = -ZRR
      IF (ZI.GT.0.0D0) GO TO 40
      ZNR = -ZNR
   40 CONTINUE
      CALL ZUNHJ(ZNR, ZNI, GNU, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      AARG = ZABS(ARGR,ARGI)
   50 CONTINUE
      IF (KODE.EQ.1) GO TO 60
      CZR = CZR - ZBR
      CZI = CZI - ZBI
   60 CONTINUE
      IF (IKFLG.EQ.1) GO TO 70
      CZR = -CZR
      CZI = -CZI
   70 CONTINUE
      APHI = ZABS(PHIR,PHII)
      RCZ = CZR
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      IF (RCZ.GT.ELIM) GO TO 210
      IF (RCZ.LT.ALIM) GO TO 80
      RCZ = RCZ + LOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*LOG(AARG) - AIC
      IF (RCZ.GT.ELIM) GO TO 210
      GO TO 130
   80 CONTINUE
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      IF (RCZ.LT.(-ELIM)) GO TO 90
      IF (RCZ.GT.(-ALIM)) GO TO 130
      RCZ = RCZ + LOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*LOG(AARG) - AIC
      IF (RCZ.GT.(-ELIM)) GO TO 110
   90 CONTINUE
      DO 100 I=1,NN
        YR(I) = ZEROR
        YI(I) = ZEROI
  100 CONTINUE
      NUF = NN
      RETURN
  110 CONTINUE
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CALL ZLOG(PHIR, PHII, STR, STI, IDUM)
      CZR = CZR + STR
      CZI = CZI + STI
      IF (IFORM.EQ.1) GO TO 120
      CALL ZLOG(ARGR, ARGI, STR, STI, IDUM)
      CZR = CZR - 0.25D0*STR - AIC
      CZI = CZI - 0.25D0*STI
  120 CONTINUE
      AX = EXP(RCZ)/TOL
      AY = CZI
      CZR = AX*COS(AY)
      CZI = AX*SIN(AY)
      CALL ZUCHK(CZR, CZI, NW, ASCLE, TOL)
      IF (NW.NE.0) GO TO 90
  130 CONTINUE
      IF (IKFLG.EQ.2) RETURN
      IF (N.EQ.1) RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOWS ON I SEQUENCE
C-----------------------------------------------------------------------
  140 CONTINUE
      GNU = FNU + (NN-1)
      IF (IFORM.EQ.2) GO TO 150
      INIT = 0
      CALL ZUNIK(ZRR, ZRI, GNU, IKFLG, 1, TOL, INIT, PHIR, PHII,
     * ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      GO TO 160
  150 CONTINUE
      CALL ZUNHJ(ZNR, ZNI, GNU, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      AARG = ZABS(ARGR,ARGI)
  160 CONTINUE
      IF (KODE.EQ.1) GO TO 170
      CZR = CZR - ZBR
      CZI = CZI - ZBI
  170 CONTINUE
      APHI = ZABS(PHIR,PHII)
      RCZ = CZR
      IF (RCZ.LT.(-ELIM)) GO TO 180
      IF (RCZ.GT.(-ALIM)) RETURN
      RCZ = RCZ + LOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*LOG(AARG) - AIC
      IF (RCZ.GT.(-ELIM)) GO TO 190
  180 CONTINUE
      YR(NN) = ZEROR
      YI(NN) = ZEROI
      NN = NN - 1
      NUF = NUF + 1
      IF (NN.EQ.0) RETURN
      GO TO 140
  190 CONTINUE
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CALL ZLOG(PHIR, PHII, STR, STI, IDUM)
      CZR = CZR + STR
      CZI = CZI + STI
      IF (IFORM.EQ.1) GO TO 200
      CALL ZLOG(ARGR, ARGI, STR, STI, IDUM)
      CZR = CZR - 0.25D0*STR - AIC
      CZI = CZI - 0.25D0*STI
  200 CONTINUE
      AX = EXP(RCZ)/TOL
      AY = CZI
      CZR = AX*COS(AY)
      CZI = AX*SIN(AY)
      CALL ZUCHK(CZR, CZI, NW, ASCLE, TOL)
      IF (NW.NE.0) GO TO 180
      RETURN
  210 CONTINUE
      NUF = -1
      RETURN
      END


*DECK ZMLRI
      SUBROUTINE ZMLRI (ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL)
C***BEGIN PROLOGUE  ZMLRI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CMLRI-A, ZMLRI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE
C     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, DGAMLN, ZABS, ZEXP, ZLOG, ZMLT
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZMLRI
C     COMPLEX CK,CNORM,CONE,CTWO,CZERO,PT,P1,P2,RZ,SUM,Y,Z
      DOUBLE PRECISION ACK, AK, AP, AT, AZ, BK, CKI, CKR, CNORMI,
     * CNORMR, CONEI, CONER, FKAP, FKK, FLAM, FNF, FNU, PTI, PTR, P1I,
     * P1R, P2I, P2R, RAZ, RHO, RHO2, RZI, RZR, SCLE, STI, STR, SUMI,
     * SUMR, TFNF, TOL, TST, YI, YR, ZEROI, ZEROR, ZI, ZR, DGAMLN,
     * D1MACH, ZABS
      INTEGER I, IAZ, IDUM, IFNU, INU, ITIME, K, KK, KM, KODE, M, N, NZ
      DIMENSION YR(N), YI(N)
      EXTERNAL ZABS, ZEXP, ZLOG
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZMLRI
      SCLE = D1MACH(1)/TOL
      NZ=0
      AZ = ZABS(ZR,ZI)
      IAZ = AZ
      IFNU = FNU
      INU = IFNU + N - 1
      AT = IAZ + 1.0D0
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      CKR = STR*AT*RAZ
      CKI = STI*AT*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      P1R = ZEROR
      P1I = ZEROI
      P2R = CONER
      P2I = CONEI
      ACK = (AT+1.0D0)*RAZ
      RHO = ACK + SQRT(ACK*ACK-1.0D0)
      RHO2 = RHO*RHO
      TST = (RHO2+RHO2)/((RHO2-1.0D0)*(RHO-1.0D0))
      TST = TST/TOL
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
C-----------------------------------------------------------------------
      AK = AT
      DO 10 I=1,80
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR*PTR-CKI*PTI)
        P2I = P1I - (CKI*PTR+CKR*PTI)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = ZABS(P2R,P2I)
        IF (AP.GT.TST*AK*AK) GO TO 20
        AK = AK + 1.0D0
   10 CONTINUE
      GO TO 110
   20 CONTINUE
      I = I + 1
      K = 0
      IF (INU.LT.IAZ) GO TO 40
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
C-----------------------------------------------------------------------
      P1R = ZEROR
      P1I = ZEROI
      P2R = CONER
      P2I = CONEI
      AT = INU + 1.0D0
      STR = ZR*RAZ
      STI = -ZI*RAZ
      CKR = STR*AT*RAZ
      CKI = STI*AT*RAZ
      ACK = AT*RAZ
      TST = SQRT(ACK/TOL)
      ITIME = 1
      DO 30 K=1,80
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR*PTR-CKI*PTI)
        P2I = P1I - (CKR*PTI+CKI*PTR)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = ZABS(P2R,P2I)
        IF (AP.LT.TST) GO TO 30
        IF (ITIME.EQ.2) GO TO 40
        ACK = ZABS(CKR,CKI)
        FLAM = ACK + SQRT(ACK*ACK-1.0D0)
        FKAP = AP/ZABS(P1R,P1I)
        RHO = MIN(FLAM,FKAP)
        TST = TST*SQRT(RHO/(RHO*RHO-1.0D0))
        ITIME = 2
   30 CONTINUE
      GO TO 110
   40 CONTINUE
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
C-----------------------------------------------------------------------
      K = K + 1
      KK = MAX(I+IAZ,K+INU)
      FKK = KK
      P1R = ZEROR
      P1I = ZEROI
C-----------------------------------------------------------------------
C     SCALE P2 AND SUM BY SCLE
C-----------------------------------------------------------------------
      P2R = SCLE
      P2I = ZEROI
      FNF = FNU - IFNU
      TFNF = FNF + FNF
cSAP091104
c      write(*,*)'in ZMLRI FKK,TFNF,IDU',FKK,TFNF,IDU

      BK = DGAMLN(FKK+TFNF+1.0D0,IDUM) - DGAMLN(FKK+1.0D0,IDUM) -
     * DGAMLN(TFNF+1.0D0,IDUM)
      BK = EXP(BK)
      SUMR = ZEROR
      SUMI = ZEROI
      KM = KK - INU
      DO 50 I=1,KM
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
   50 CONTINUE
      YR(N) = P2R
      YI(N) = P2I
      IF (N.EQ.1) GO TO 70
      DO 60 I=2,N
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
        M = N - I + 1
        YR(M) = P2R
        YI(M) = P2I
   60 CONTINUE
   70 CONTINUE
      IF (IFNU.LE.0) GO TO 90
      DO 80 I=1,IFNU
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZR*PTI+RZI*PTR)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
   80 CONTINUE
   90 CONTINUE
      PTR = ZR
      PTI = ZI
      IF (KODE.EQ.2) PTR = ZEROR
      CALL ZLOG(RZR, RZI, STR, STI, IDUM)
      P1R = -FNF*STR + PTR
      P1I = -FNF*STI + PTI
cSAP091104
c      write(*,*)'ZMLRI FNF,IDUM',FNF,IDUM

      AP = DGAMLN(1.0D0+FNF,IDUM)
      PTR = P1R - AP
      PTI = P1I
C-----------------------------------------------------------------------
C     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
C     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
C-----------------------------------------------------------------------
      P2R = P2R + SUMR
      P2I = P2I + SUMI
      AP = ZABS(P2R,P2I)
      P1R = 1.0D0/AP
      CALL ZEXP(PTR, PTI, STR, STI)
      CKR = STR*P1R
      CKI = STI*P1R
      PTR = P2R*P1R
      PTI = -P2I*P1R
      CALL ZMLT(CKR, CKI, PTR, PTI, CNORMR, CNORMI)
      DO 100 I=1,N
        STR = YR(I)*CNORMR - YI(I)*CNORMI
        YI(I) = YR(I)*CNORMI + YI(I)*CNORMR
        YR(I) = STR
  100 CONTINUE
      RETURN
  110 CONTINUE
      NZ=-2
      RETURN
      END

*DECK ZWRSK
      SUBROUTINE ZWRSK (ZRR, ZRI, FNU, KODE, N, YR, YI, NZ, CWR, CWI,
     +   TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZWRSK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CWRSK-A, ZWRSK-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY
C     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZBKNU, ZRATI
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZWRSK
C     COMPLEX CINU,CSCL,CT,CW,C1,C2,RCT,ST,Y,ZR
      DOUBLE PRECISION ACT, ACW, ALIM, ASCLE, CINUI, CINUR, CSCLR, CTI,
     * CTR, CWI, CWR, C1I, C1R, C2I, C2R, ELIM, FNU, PTI, PTR, RACT,
     * STI, STR, TOL, YI, YR, ZRI, ZRR, ZABS, D1MACH
      INTEGER I, KODE, N, NW, NZ
      DIMENSION YR(N), YI(N), CWR(2), CWI(2)
      EXTERNAL ZABS
C***FIRST EXECUTABLE STATEMENT  ZWRSK
C-----------------------------------------------------------------------
C     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
C     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
C     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
C-----------------------------------------------------------------------
C
      NZ = 0
      CALL ZBKNU(ZRR, ZRI, FNU, KODE, 2, CWR, CWI, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 50
      CALL ZRATI(ZRR, ZRI, FNU, N, YR, YI, TOL)
C-----------------------------------------------------------------------
C     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
C     R(FNU+J-1,Z)=Y(J),  J=1,...,N
C-----------------------------------------------------------------------
      CINUR = 1.0D0
      CINUI = 0.0D0
      IF (KODE.EQ.1) GO TO 10
      CINUR = COS(ZRI)
      CINUI = SIN(ZRI)
   10 CONTINUE
C-----------------------------------------------------------------------
C     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
C     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
C     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
C     THE RESULT IS ON SCALE.
C-----------------------------------------------------------------------
      ACW = ZABS(CWR(2),CWI(2))
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CSCLR = 1.0D0
      IF (ACW.GT.ASCLE) GO TO 20
      CSCLR = 1.0D0/TOL
      GO TO 30
   20 CONTINUE
      ASCLE = 1.0D0/ASCLE
      IF (ACW.LT.ASCLE) GO TO 30
      CSCLR = TOL
   30 CONTINUE
      C1R = CWR(1)*CSCLR
      C1I = CWI(1)*CSCLR
      C2R = CWR(2)*CSCLR
      C2I = CWI(2)*CSCLR
      STR = YR(1)
      STI = YI(1)
C-----------------------------------------------------------------------
C     CINU=CINU*(CONJG(CT)/ABS(CT))*(1.0D0/ABS(CT) PREVENTS
C     UNDER- OR OVERFLOW PREMATURELY BY SQUARING ABS(CT)
C-----------------------------------------------------------------------
      PTR = STR*C1R - STI*C1I
      PTI = STR*C1I + STI*C1R
      PTR = PTR + C2R
      PTI = PTI + C2I
      CTR = ZRR*PTR - ZRI*PTI
      CTI = ZRR*PTI + ZRI*PTR
      ACT = ZABS(CTR,CTI)
      RACT = 1.0D0/ACT
      CTR = CTR*RACT
      CTI = -CTI*RACT
      PTR = CINUR*RACT
      PTI = CINUI*RACT
      CINUR = PTR*CTR - PTI*CTI
      CINUI = PTR*CTI + PTI*CTR
      YR(1) = CINUR*CSCLR
      YI(1) = CINUI*CSCLR
      IF (N.EQ.1) RETURN
      DO 40 I=2,N
        PTR = STR*CINUR - STI*CINUI
        CINUI = STR*CINUI + STI*CINUR
        CINUR = PTR
        STR = YR(I)
        STI = YI(I)
        YR(I) = CINUR*CSCLR
        YI(I) = CINUI*CSCLR
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END


*DECK ZBUNI
      SUBROUTINE ZBUNI (ZR, ZI, FNU, KODE, N, YR, YI, NZ, NUI, NLAST,
     +   FNUL, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZBUNI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CBUNI-A, ZBUNI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE ABS(Z).GT.
C     FNUL AND FNU+N-1.LT.FNUL. THE ORDER IS INCREASED FROM
C     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
C     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZUNI1, ZUNI2
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZBUNI
C     COMPLEX CSCL,CSCR,CY,RZ,ST,S1,S2,Y,Z
      DOUBLE PRECISION ALIM, AX, AY, CSCLR, CSCRR, CYI, CYR, DFNU,
     * ELIM, FNU, FNUI, FNUL, GNU, RAZ, RZI, RZR, STI, STR, S1I, S1R,
     * S2I, S2R, TOL, YI, YR, ZI, ZR, ZABS, ASCLE, BRY, C1R, C1I, C1M,
     * D1MACH
      INTEGER I, IFLAG, IFORM, K, KODE, N, NL, NLAST, NUI, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2), BRY(3)
      EXTERNAL ZABS
C***FIRST EXECUTABLE STATEMENT  ZBUNI
      NZ = 0
      AX = ABS(ZR)*1.7321D0
      AY = ABS(ZI)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      IF (NUI.EQ.0) GO TO 60
      FNUI = NUI
      DFNU = FNU + (N-1)
      GNU = DFNU + FNUI
      IF (IFORM.EQ.2) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL ZUNI1(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
      GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL ZUNI2(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
   20 CONTINUE
      IF (NW.LT.0) GO TO 50
      IF (NW.NE.0) GO TO 90
      STR = ZABS(CYR(1),CYI(1))
C----------------------------------------------------------------------
C     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
C----------------------------------------------------------------------
      BRY(1)=1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = BRY(2)
      IFLAG = 2
      ASCLE = BRY(2)
      CSCLR = 1.0D0
      IF (STR.GT.BRY(1)) GO TO 21
      IFLAG = 1
      ASCLE = BRY(1)
      CSCLR = 1.0D0/TOL
      GO TO 25
   21 CONTINUE
      IF (STR.LT.BRY(2)) GO TO 25
      IFLAG = 3
      ASCLE=BRY(3)
      CSCLR = TOL
   25 CONTINUE
      CSCRR = 1.0D0/CSCLR
      S1R = CYR(2)*CSCLR
      S1I = CYI(2)*CSCLR
      S2R = CYR(1)*CSCLR
      S2I = CYI(1)*CSCLR
      RAZ = 1.0D0/ZABS(ZR,ZI)
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      DO 30 I=1,NUI
        STR = S2R
        STI = S2I
        S2R = (DFNU+FNUI)*(RZR*STR-RZI*STI) + S1R
        S2I = (DFNU+FNUI)*(RZR*STI+RZI*STR) + S1I
        S1R = STR
        S1I = STI
        FNUI = FNUI - 1.0D0
        IF (IFLAG.GE.3) GO TO 30
        STR = S2R*CSCRR
        STI = S2I*CSCRR
        C1R = ABS(STR)
        C1I = ABS(STI)
        C1M = MAX(C1R,C1I)
        IF (C1M.LE.ASCLE) GO TO 30
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSCRR
        S1I = S1I*CSCRR
        S2R = STR
        S2I = STI
        CSCLR = CSCLR*TOL
        CSCRR = 1.0D0/CSCLR
        S1R = S1R*CSCLR
        S1I = S1I*CSCLR
        S2R = S2R*CSCLR
        S2I = S2I*CSCLR
   30 CONTINUE
      YR(N) = S2R*CSCRR
      YI(N) = S2I*CSCRR
      IF (N.EQ.1) RETURN
      NL = N - 1
      FNUI = NL
      K = NL
      DO 40 I=1,NL
        STR = S2R
        STI = S2I
        S2R = (FNU+FNUI)*(RZR*STR-RZI*STI) + S1R
        S2I = (FNU+FNUI)*(RZR*STI+RZI*STR) + S1I
        S1R = STR
        S1I = STI
        STR = S2R*CSCRR
        STI = S2I*CSCRR
        YR(K) = STR
        YI(K) = STI
        FNUI = FNUI - 1.0D0
        K = K - 1
        IF (IFLAG.GE.3) GO TO 40
        C1R = ABS(STR)
        C1I = ABS(STI)
        C1M = MAX(C1R,C1I)
        IF (C1M.LE.ASCLE) GO TO 40
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSCRR
        S1I = S1I*CSCRR
        S2R = STR
        S2I = STI
        CSCLR = CSCLR*TOL
        CSCRR = 1.0D0/CSCLR
        S1R = S1R*CSCLR
        S1I = S1I*CSCLR
        S2R = S2R*CSCLR
        S2I = S2I*CSCLR
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
   60 CONTINUE
      IF (IFORM.EQ.2) GO TO 70
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
      GO TO 80
   70 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
   80 CONTINUE
      IF (NW.LT.0) GO TO 50
      NZ = NW
      RETURN
   90 CONTINUE
      NLAST = N
      RETURN
      END


*DECK ZBKNU
      SUBROUTINE ZBKNU (ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  ZBKNU
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZAIRY, ZBESH, ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CBKNU-A, ZBKNU-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE.
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, DGAMLN, I1MACH, ZABS, ZDIV, ZEXP, ZKSCL,
C                    ZLOG, ZMLT, ZSHCH, ZSQRT, ZUCHK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZEXP, ZLOG and ZSQRT to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZBKNU
C
      DOUBLE PRECISION AA, AK, ALIM, ASCLE, A1, A2, BB, BK, BRY, CAZ,
     * CBI, CBR, CC, CCHI, CCHR, CKI, CKR, COEFI, COEFR, CONEI, CONER,
     * CRSCR, CSCLR, CSHI, CSHR, CSI, CSR, CSRR, CSSR, CTWOR,
     * CZEROI, CZEROR, CZI, CZR, DNU, DNU2, DPI, ELIM, ETEST, FC, FHS,
     * FI, FK, FKS, FMUI, FMUR, FNU, FPI, FR, G1, G2, HPI, PI, PR, PTI,
     * PTR, P1I, P1R, P2I, P2M, P2R, QI, QR, RAK, RCAZ, RTHPI, RZI,
     * RZR, R1, S, SMUI, SMUR, SPI, STI, STR, S1I, S1R, S2I, S2R, TM,
     * TOL, TTH, T1, T2, YI, YR, ZI, ZR, DGAMLN, D1MACH, ZABS, ELM,
     * CELMR, ZDR, ZDI, AS, ALAS, HELIM, CYR, CYI
      INTEGER I, IFLAG, INU, K, KFLAG, KK, KMAX, KODE, KODED, N, NZ,
     * IDUM, I1MACH, J, IC, INUB, NW
      DIMENSION YR(N), YI(N), CC(8), CSSR(3), CSRR(3), BRY(3), CYR(2),
     * CYI(2)
      EXTERNAL ZABS, ZEXP, ZLOG, ZSQRT
C     COMPLEX Z,Y,A,B,RZ,SMU,FU,FMU,F,FLRZ,CZ,S1,S2,CSH,CCH
C     COMPLEX CK,P,Q,COEF,P1,P2,CBK,PT,CZERO,CONE,CTWO,ST,EZ,CS,DK
C
      DATA KMAX / 30 /
      DATA CZEROR,CZEROI,CONER,CONEI,CTWOR,R1/
     1  0.0D0 , 0.0D0 , 1.0D0 , 0.0D0 , 2.0D0 , 2.0D0 /
      DATA DPI, RTHPI, SPI ,HPI, FPI, TTH /
     1     3.14159265358979324D0,       1.25331413731550025D0,
     2     1.90985931710274403D0,       1.57079632679489662D0,
     3     1.89769999331517738D0,       6.66666666666666666D-01/
      DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)/
     1     5.77215664901532861D-01,    -4.20026350340952355D-02,
     2    -4.21977345555443367D-02,     7.21894324666309954D-03,
     3    -2.15241674114950973D-04,    -2.01348547807882387D-05,
     4     1.13302723198169588D-06,     6.11609510448141582D-09/
C***FIRST EXECUTABLE STATEMENT  ZBKNU
      CAZ = ZABS(ZR,ZI)
      CSCLR = 1.0D0/TOL
      CRSCR = TOL
      CSSR(1) = CSCLR
      CSSR(2) = 1.0D0
      CSSR(3) = CRSCR
      CSRR(1) = CRSCR
      CSRR(2) = 1.0D0
      CSRR(3) = CSCLR
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      NZ = 0
      IFLAG = 0
      KODED = KODE
      RCAZ = 1.0D0/CAZ
      STR = ZR*RCAZ
      STI = -ZI*RCAZ
      RZR = (STR+STR)*RCAZ
      RZI = (STI+STI)*RCAZ
      INU = FNU+0.5D0
      DNU = FNU - INU
      IF (ABS(DNU).EQ.0.5D0) GO TO 110
      DNU2 = 0.0D0
      IF (ABS(DNU).GT.TOL) DNU2 = DNU*DNU
      IF (CAZ.GT.R1) GO TO 110
C-----------------------------------------------------------------------
C     SERIES FOR ABS(Z).LE.R1
C-----------------------------------------------------------------------
      FC = 1.0D0
      CALL ZLOG(RZR, RZI, SMUR, SMUI, IDUM)
      FMUR = SMUR*DNU
      FMUI = SMUI*DNU
      CALL ZSHCH(FMUR, FMUI, CSHR, CSHI, CCHR, CCHI)
      IF (DNU.EQ.0.0D0) GO TO 10
      FC = DNU*DPI
      FC = FC/SIN(FC)
      SMUR = CSHR/DNU
      SMUI = CSHI/DNU
   10 CONTINUE
      A2 = 1.0D0 + DNU
C-----------------------------------------------------------------------
C     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
C-----------------------------------------------------------------------
cSAP091104
c      write(*,*)'in ZBKNU A2,IDUM',A2,IDUM

      T2 = EXP(-DGAMLN(A2,IDUM))
      T1 = 1.0D0/(T2*FC)
      IF (ABS(DNU).GT.0.1D0) GO TO 40
C-----------------------------------------------------------------------
C     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
C-----------------------------------------------------------------------
      AK = 1.0D0
      S = CC(1)
      DO 20 K=2,8
        AK = AK*DNU2
        TM = CC(K)*AK
        S = S + TM
        IF (ABS(TM).LT.TOL) GO TO 30
   20 CONTINUE
   30 G1 = -S
      GO TO 50
   40 CONTINUE
      G1 = (T1-T2)/(DNU+DNU)
   50 CONTINUE
      G2 = (T1+T2)*0.5D0
      FR = FC*(CCHR*G1+SMUR*G2)
      FI = FC*(CCHI*G1+SMUI*G2)
      CALL ZEXP(FMUR, FMUI, STR, STI)
      PR = 0.5D0*STR/T2
      PI = 0.5D0*STI/T2
      CALL ZDIV(0.5D0, 0.0D0, STR, STI, PTR, PTI)
      QR = PTR/T1
      QI = PTI/T1
      S1R = FR
      S1I = FI
      S2R = PR
      S2I = PI
      AK = 1.0D0
      A1 = 1.0D0
      CKR = CONER
      CKI = CONEI
      BK = 1.0D0 - DNU2
      IF (INU.GT.0 .OR. N.GT.1) GO TO 80
C-----------------------------------------------------------------------
C     GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1
C-----------------------------------------------------------------------
      IF (CAZ.LT.TOL) GO TO 70
      CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
      CZR = 0.25D0*CZR
      CZI = 0.25D0*CZI
      T1 = 0.25D0*CAZ*CAZ
   60 CONTINUE
      FR = (FR*AK+PR+QR)/BK
      FI = (FI*AK+PI+QI)/BK
      STR = 1.0D0/(AK-DNU)
      PR = PR*STR
      PI = PI*STR
      STR = 1.0D0/(AK+DNU)
      QR = QR*STR
      QI = QI*STR
      STR = CKR*CZR - CKI*CZI
      RAK = 1.0D0/AK
      CKI = (CKR*CZI+CKI*CZR)*RAK
      CKR = STR*RAK
      S1R = CKR*FR - CKI*FI + S1R
      S1I = CKR*FI + CKI*FR + S1I
      A1 = A1*T1*RAK
      BK = BK + AK + AK + 1.0D0
      AK = AK + 1.0D0
      IF (A1.GT.TOL) GO TO 60
   70 CONTINUE
      YR(1) = S1R
      YI(1) = S1I
      IF (KODED.EQ.1) RETURN
      CALL ZEXP(ZR, ZI, STR, STI)
      CALL ZMLT(S1R, S1I, STR, STI, YR(1), YI(1))
      RETURN
C-----------------------------------------------------------------------
C     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
C-----------------------------------------------------------------------
   80 CONTINUE
      IF (CAZ.LT.TOL) GO TO 100
      CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
      CZR = 0.25D0*CZR
      CZI = 0.25D0*CZI
      T1 = 0.25D0*CAZ*CAZ
   90 CONTINUE
      FR = (FR*AK+PR+QR)/BK
      FI = (FI*AK+PI+QI)/BK
      STR = 1.0D0/(AK-DNU)
      PR = PR*STR
      PI = PI*STR
      STR = 1.0D0/(AK+DNU)
      QR = QR*STR
      QI = QI*STR
      STR = CKR*CZR - CKI*CZI
      RAK = 1.0D0/AK
      CKI = (CKR*CZI+CKI*CZR)*RAK
      CKR = STR*RAK
      S1R = CKR*FR - CKI*FI + S1R
      S1I = CKR*FI + CKI*FR + S1I
      STR = PR - FR*AK
      STI = PI - FI*AK
      S2R = CKR*STR - CKI*STI + S2R
      S2I = CKR*STI + CKI*STR + S2I
      A1 = A1*T1*RAK
      BK = BK + AK + AK + 1.0D0
      AK = AK + 1.0D0
      IF (A1.GT.TOL) GO TO 90
  100 CONTINUE
      KFLAG = 2
      A1 = FNU + 1.0D0
      AK = A1*ABS(SMUR)
      IF (AK.GT.ALIM) KFLAG = 3
      STR = CSSR(KFLAG)
      P2R = S2R*STR
      P2I = S2I*STR
      CALL ZMLT(P2R, P2I, RZR, RZI, S2R, S2I)
      S1R = S1R*STR
      S1I = S1I*STR
      IF (KODED.EQ.1) GO TO 210
      CALL ZEXP(ZR, ZI, FR, FI)
      CALL ZMLT(S1R, S1I, FR, FI, S1R, S1I)
      CALL ZMLT(S2R, S2I, FR, FI, S2R, S2I)
      GO TO 210
C-----------------------------------------------------------------------
C     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
C     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
C     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
C     RECURSION
C-----------------------------------------------------------------------
  110 CONTINUE
      CALL ZSQRT(ZR, ZI, STR, STI)
      CALL ZDIV(RTHPI, CZEROI, STR, STI, COEFR, COEFI)
      KFLAG = 2
      IF (KODED.EQ.2) GO TO 120
      IF (ZR.GT.ALIM) GO TO 290
C     BLANK LINE
      STR = EXP(-ZR)*CSSR(KFLAG)
      STI = -STR*SIN(ZI)
      STR = STR*COS(ZI)
      CALL ZMLT(COEFR, COEFI, STR, STI, COEFR, COEFI)
  120 CONTINUE
      IF (ABS(DNU).EQ.0.5D0) GO TO 300
C-----------------------------------------------------------------------
C     MILLER ALGORITHM FOR ABS(Z).GT.R1
C-----------------------------------------------------------------------
      AK = COS(DPI*DNU)
      AK = ABS(AK)
      IF (AK.EQ.CZEROR) GO TO 300
      FHS = ABS(0.25D0-DNU2)
      IF (FHS.EQ.CZEROR) GO TO 300
C-----------------------------------------------------------------------
C     COMPUTE R2=F(E). IF ABS(Z).GE.R2, USE FORWARD RECURRENCE TO
C     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
C     12.LE.E.LE.60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(14))=
C     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
C-----------------------------------------------------------------------
      T1 = I1MACH(14)-1
      T1 = T1*D1MACH(5)*3.321928094D0
      T1 = MAX(T1,12.0D0)
      T1 = MIN(T1,60.0D0)
      T2 = TTH*T1 - 6.0D0
      IF (ZR.NE.0.0D0) GO TO 130
      T1 = HPI
      GO TO 140
  130 CONTINUE
      T1 = DATAN(ZI/ZR)
      T1 = ABS(T1)
  140 CONTINUE
      IF (T2.GT.CAZ) GO TO 170
C-----------------------------------------------------------------------
C     FORWARD RECURRENCE LOOP WHEN ABS(Z).GE.R2
C-----------------------------------------------------------------------
      ETEST = AK/(DPI*CAZ*TOL)
      FK = CONER
      IF (ETEST.LT.CONER) GO TO 180
      FKS = CTWOR
      CKR = CAZ + CAZ + CTWOR
      P1R = CZEROR
      P2R = CONER
      DO 150 I=1,KMAX
        AK = FHS/FKS
        CBR = CKR/(FK+CONER)
        PTR = P2R
        P2R = CBR*P2R - P1R*AK
        P1R = PTR
        CKR = CKR + CTWOR
        FKS = FKS + FK + FK + CTWOR
        FHS = FHS + FK + FK
        FK = FK + CONER
        STR = ABS(P2R)*FK
        IF (ETEST.LT.STR) GO TO 160
  150 CONTINUE
      GO TO 310
  160 CONTINUE
      FK = FK + SPI*T1*SQRT(T2/CAZ)
      FHS = ABS(0.25D0-DNU2)
      GO TO 180
  170 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE BACKWARD INDEX K FOR ABS(Z).LT.R2
C-----------------------------------------------------------------------
      A2 = SQRT(CAZ)
      AK = FPI*AK/(TOL*SQRT(A2))
      AA = 3.0D0*T1/(1.0D0+CAZ)
      BB = 14.7D0*T1/(28.0D0+CAZ)
      AK = (LOG(AK)+CAZ*COS(AA)/(1.0D0+0.008D0*CAZ))/COS(BB)
      FK = 0.12125D0*AK*AK/CAZ + 1.5D0
  180 CONTINUE
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
      K = FK
      FK = K
      FKS = FK*FK
      P1R = CZEROR
      P1I = CZEROI
      P2R = TOL
      P2I = CZEROI
      CSR = P2R
      CSI = P2I
      DO 190 I=1,K
        A1 = FKS - FK
        AK = (FKS+FK)/(A1+FHS)
        RAK = 2.0D0/(FK+CONER)
        CBR = (FK+ZR)*RAK
        CBI = ZI*RAK
        PTR = P2R
        PTI = P2I
        P2R = (PTR*CBR-PTI*CBI-P1R)*AK
        P2I = (PTI*CBR+PTR*CBI-P1I)*AK
        P1R = PTR
        P1I = PTI
        CSR = CSR + P2R
        CSI = CSI + P2I
        FKS = A1 - FK + CONER
        FK = FK - CONER
  190 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE (P2/CS)=(P2/ABS(CS))*(CONJG(CS)/ABS(CS)) FOR BETTER
C     SCALING
C-----------------------------------------------------------------------
      TM = ZABS(CSR,CSI)
      PTR = 1.0D0/TM
      S1R = P2R*PTR
      S1I = P2I*PTR
      CSR = CSR*PTR
      CSI = -CSI*PTR
      CALL ZMLT(COEFR, COEFI, S1R, S1I, STR, STI)
      CALL ZMLT(STR, STI, CSR, CSI, S1R, S1I)
      IF (INU.GT.0 .OR. N.GT.1) GO TO 200
      ZDR = ZR
      ZDI = ZI
      IF(IFLAG.EQ.1) GO TO 270
      GO TO 240
  200 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE P1/P2=(P1/ABS(P2)*CONJG(P2)/ABS(P2) FOR SCALING
C-----------------------------------------------------------------------
      TM = ZABS(P2R,P2I)
      PTR = 1.0D0/TM
      P1R = P1R*PTR
      P1I = P1I*PTR
      P2R = P2R*PTR
      P2I = -P2I*PTR
      CALL ZMLT(P1R, P1I, P2R, P2I, PTR, PTI)
      STR = DNU + 0.5D0 - PTR
      STI = -PTI
      CALL ZDIV(STR, STI, ZR, ZI, STR, STI)
      STR = STR + 1.0D0
      CALL ZMLT(STR, STI, S1R, S1I, S2R, S2I)
C-----------------------------------------------------------------------
C     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
C     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
C-----------------------------------------------------------------------
  210 CONTINUE
      STR = DNU + 1.0D0
      CKR = STR*RZR
      CKI = STR*RZI
      IF (N.EQ.1) INU = INU - 1
      IF (INU.GT.0) GO TO 220
      IF (N.GT.1) GO TO 215
      S1R = S2R
      S1I = S2I
  215 CONTINUE
      ZDR = ZR
      ZDI = ZI
      IF(IFLAG.EQ.1) GO TO 270
      GO TO 240
  220 CONTINUE
      INUB = 1
      IF(IFLAG.EQ.1) GO TO 261
  225 CONTINUE
      P1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 230 I=INUB,INU
        STR = S2R
        STI = S2I
        S2R = CKR*STR - CKI*STI + S1R
        S2I = CKR*STI + CKI*STR + S1I
        S1R = STR
        S1I = STI
        CKR = CKR + RZR
        CKI = CKI + RZI
        IF (KFLAG.GE.3) GO TO 230
        P2R = S2R*P1R
        P2I = S2I*P1R
        STR = ABS(P2R)
        STI = ABS(P2I)
        P2M = MAX(STR,STI)
        IF (P2M.LE.ASCLE) GO TO 230
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*P1R
        S1I = S1I*P1R
        S2R = P2R
        S2I = P2I
        STR = CSSR(KFLAG)
        S1R = S1R*STR
        S1I = S1I*STR
        S2R = S2R*STR
        S2I = S2I*STR
        P1R = CSRR(KFLAG)
  230 CONTINUE
      IF (N.NE.1) GO TO 240
      S1R = S2R
      S1I = S2I
  240 CONTINUE
      STR = CSRR(KFLAG)
      YR(1) = S1R*STR
      YI(1) = S1I*STR
      IF (N.EQ.1) RETURN
      YR(2) = S2R*STR
      YI(2) = S2I*STR
      IF (N.EQ.2) RETURN
      KK = 2
  250 CONTINUE
      KK = KK + 1
      IF (KK.GT.N) RETURN
      P1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 260 I=KK,N
        P2R = S2R
        P2I = S2I
        S2R = CKR*P2R - CKI*P2I + S1R
        S2I = CKI*P2R + CKR*P2I + S1I
        S1R = P2R
        S1I = P2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        P2R = S2R*P1R
        P2I = S2I*P1R
        YR(I) = P2R
        YI(I) = P2I
        IF (KFLAG.GE.3) GO TO 260
        STR = ABS(P2R)
        STI = ABS(P2I)
        P2M = MAX(STR,STI)
        IF (P2M.LE.ASCLE) GO TO 260
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*P1R
        S1I = S1I*P1R
        S2R = P2R
        S2I = P2I
        STR = CSSR(KFLAG)
        S1R = S1R*STR
        S1I = S1I*STR
        S2R = S2R*STR
        S2I = S2I*STR
        P1R = CSRR(KFLAG)
  260 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
C-----------------------------------------------------------------------
  261 CONTINUE
      HELIM = 0.5D0*ELIM
      ELM = EXP(-ELIM)
      CELMR = ELM
      ASCLE = BRY(1)
      ZDR = ZR
      ZDI = ZI
      IC = -1
      J = 2
      DO 262 I=1,INU
        STR = S2R
        STI = S2I
        S2R = STR*CKR-STI*CKI+S1R
        S2I = STI*CKR+STR*CKI+S1I
        S1R = STR
        S1I = STI
        CKR = CKR+RZR
        CKI = CKI+RZI
        AS = ZABS(S2R,S2I)
        ALAS = LOG(AS)
        P2R = -ZDR+ALAS
        IF(P2R.LT.(-ELIM)) GO TO 263
        CALL ZLOG(S2R,S2I,STR,STI,IDUM)
        P2R = -ZDR+STR
        P2I = -ZDI+STI
        P2M = EXP(P2R)/TOL
        P1R = P2M*COS(P2I)
        P1I = P2M*SIN(P2I)
        CALL ZUCHK(P1R,P1I,NW,ASCLE,TOL)
        IF(NW.NE.0) GO TO 263
        J = 3 - J
        CYR(J) = P1R
        CYI(J) = P1I
        IF(IC.EQ.(I-1)) GO TO 264
        IC = I
        GO TO 262
  263   CONTINUE
        IF(ALAS.LT.HELIM) GO TO 262
        ZDR = ZDR-ELIM
        S1R = S1R*CELMR
        S1I = S1I*CELMR
        S2R = S2R*CELMR
        S2I = S2I*CELMR
  262 CONTINUE
      IF(N.NE.1) GO TO 270
      S1R = S2R
      S1I = S2I
      GO TO 270
  264 CONTINUE
      KFLAG = 1
      INUB = I+1
      S2R = CYR(J)
      S2I = CYI(J)
      J = 3 - J
      S1R = CYR(J)
      S1I = CYI(J)
      IF(INUB.LE.INU) GO TO 225
      IF(N.NE.1) GO TO 240
      S1R = S2R
      S1I = S2I
      GO TO 240
  270 CONTINUE
      YR(1) = S1R
      YI(1) = S1I
      IF(N.EQ.1) GO TO 280
      YR(2) = S2R
      YI(2) = S2I
  280 CONTINUE
      ASCLE = BRY(1)
      CALL ZKSCL(ZDR,ZDI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM)
      INU = N - NZ
      IF (INU.LE.0) RETURN
      KK = NZ + 1
      S1R = YR(KK)
      S1I = YI(KK)
      YR(KK) = S1R*CSRR(1)
      YI(KK) = S1I*CSRR(1)
      IF (INU.EQ.1) RETURN
      KK = NZ + 2
      S2R = YR(KK)
      S2I = YI(KK)
      YR(KK) = S2R*CSRR(1)
      YI(KK) = S2I*CSRR(1)
      IF (INU.EQ.2) RETURN
      T2 = FNU + (KK-1)
      CKR = T2*RZR
      CKI = T2*RZI
      KFLAG = 1
      GO TO 250
  290 CONTINUE
C-----------------------------------------------------------------------
C     SCALE BY EXP(Z), IFLAG = 1 CASES
C-----------------------------------------------------------------------
      KODED = 2
      IFLAG = 1
      KFLAG = 2
      GO TO 120
C-----------------------------------------------------------------------
C     FNU=HALF ODD INTEGER CASE, DNU=-0.5
C-----------------------------------------------------------------------
  300 CONTINUE
      S1R = COEFR
      S1I = COEFI
      S2R = COEFR
      S2I = COEFI
      GO TO 210
C
C
  310 CONTINUE
      NZ=-2
      RETURN
      END




C-ENM added ZACON
*DECK ZACON
      SUBROUTINE ZACON (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, FNUL,
     +   TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZACON
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CACON-A, ZACON-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZACON APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE
C
C***SEE ALSO  ZBESH, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZBINU, ZBKNU, ZMLT, ZS1S2
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZACON
C     COMPLEX CK,CONE,CSCL,CSCR,CSGN,CSPN,CY,CZERO,C1,C2,RZ,SC1,SC2,ST,
C    *S1,S2,Y,Z,ZN
      DOUBLE PRECISION ALIM, ARG, ASCLE, AS2, AZN, BRY, BSCLE, CKI,
     * CKR, CONER, CPN, CSCL, CSCR, CSGNI, CSGNR, CSPNI, CSPNR,
     * CSR, CSRR, CSSR, CYI, CYR, C1I, C1M, C1R, C2I, C2R, ELIM, FMR,
     * FN, FNU, FNUL, PI, PTI, PTR, RAZN, RL, RZI, RZR, SC1I, SC1R,
     * SC2I, SC2R, SGN, SPN, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR,
     * YY, ZEROR, ZI, ZNI, ZNR, ZR, D1MACH, ZABS
      INTEGER I, INU, IUF, KFLAG, KODE, MR, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2), CSSR(3), CSRR(3), BRY(3)
      EXTERNAL ZABS
      DATA PI / 3.14159265358979324D0 /
      DATA ZEROR,CONER / 0.0D0,1.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZACON
      NZ = 0
      ZNR = -ZR
      ZNI = -ZI
      NN = N
cSAP091104
c      write(*,*)'in ZACON before ZBINU(ZNR, ZNI, FNU',ZNR, ZNI, FNU

      CALL ZBINU(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, FNUL, TOL,
     * ELIM, ALIM)
      IF (NW.LT.0) GO TO 90
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
      NN = MIN(2,N)
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 90
      S1R = CYR(1)
      S1I = CYI(1)
      FMR = MR
      SGN = -DSIGN(PI,FMR)
      CSGNR = ZEROR
      CSGNI = SGN
      IF (KODE.EQ.1) GO TO 10
      YY = -ZNI
      CPN = COS(YY)
      SPN = SIN(YY)
      CALL ZMLT(CSGNR, CSGNI, CPN, SPN, CSGNR, CSGNI)
   10 CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = FNU
      ARG = (FNU-INU)*SGN
      CPN = COS(ARG)
      SPN = SIN(ARG)
      CSPNR = CPN
      CSPNI = SPN
      IF (MOD(INU,2).EQ.0) GO TO 20
      CSPNR = -CSPNR
      CSPNI = -CSPNI
   20 CONTINUE
      IUF = 0
      C1R = S1R
      C1I = S1I
      C2R = YR(1)
      C2I = YI(1)
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      IF (KODE.EQ.1) GO TO 30
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
      SC1R = C1R
      SC1I = C1I
   30 CONTINUE
      CALL ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI)
      CALL ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI)
      YR(1) = STR + PTR
      YI(1) = STI + PTI
      IF (N.EQ.1) RETURN
      CSPNR = -CSPNR
      CSPNI = -CSPNI
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = S2R
      C1I = S2I
      C2R = YR(2)
      C2I = YI(2)
      IF (KODE.EQ.1) GO TO 40
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
      SC2R = C1R
      SC2I = C1I
   40 CONTINUE
      CALL ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI)
      CALL ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI)
      YR(2) = STR + PTR
      YI(2) = STI + PTI
      IF (N.EQ.2) RETURN
      CSPNR = -CSPNR
      CSPNI = -CSPNI
      AZN = ZABS(ZNR,ZNI)
      RAZN = 1.0D0/AZN
      STR = ZNR*RAZN
      STI = -ZNI*RAZN
      RZR = (STR+STR)*RAZN
      RZI = (STI+STI)*RAZN
      FN = FNU + 1.0D0
      CKR = FN*RZR
      CKI = FN*RZI
C-----------------------------------------------------------------------
C     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CSCR = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CSCR
      CSRR(1) = CSCR
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = ASCLE
      BRY(2) = 1.0D0/ASCLE
      BRY(3) = D1MACH(2)
      AS2 = ZABS(S2R,S2I)
      KFLAG = 2
      IF (AS2.GT.BRY(1)) GO TO 50
      KFLAG = 1
      GO TO 60
   50 CONTINUE
      IF (AS2.LT.BRY(2)) GO TO 60
      KFLAG = 3
   60 CONTINUE
      BSCLE = BRY(KFLAG)
      S1R = S1R*CSSR(KFLAG)
      S1I = S1I*CSSR(KFLAG)
      S2R = S2R*CSSR(KFLAG)
      S2I = S2I*CSSR(KFLAG)
      CSR = CSRR(KFLAG)
      DO 80 I=3,N
        STR = S2R
        STI = S2I
        S2R = CKR*STR - CKI*STI + S1R
        S2I = CKR*STI + CKI*STR + S1I
        S1R = STR
        S1I = STI
        C1R = S2R*CSR
        C1I = S2I*CSR
        STR = C1R
        STI = C1I
        C2R = YR(I)
        C2I = YI(I)
        IF (KODE.EQ.1) GO TO 70
        IF (IUF.LT.0) GO TO 70
        CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
        NZ = NZ + NW
        SC1R = SC2R
        SC1I = SC2I
        SC2R = C1R
        SC2I = C1I
        IF (IUF.NE.3) GO TO 70
        IUF = -4
        S1R = SC1R*CSSR(KFLAG)
        S1I = SC1I*CSSR(KFLAG)
        S2R = SC2R*CSSR(KFLAG)
        S2I = SC2I*CSSR(KFLAG)
        STR = SC2R
        STI = SC2I
   70   CONTINUE
        PTR = CSPNR*C1R - CSPNI*C1I
        PTI = CSPNR*C1I + CSPNI*C1R
        YR(I) = PTR + CSGNR*C2R - CSGNI*C2I
        YI(I) = PTI + CSGNR*C2I + CSGNI*C2R
        CKR = CKR + RZR
        CKI = CKI + RZI
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        IF (KFLAG.GE.3) GO TO 80
        PTR = ABS(C1R)
        PTI = ABS(C1I)
        C1M = MAX(PTR,PTI)
        IF (C1M.LE.BSCLE) GO TO 80
        KFLAG = KFLAG + 1
        BSCLE = BRY(KFLAG)
        S1R = S1R*CSR
        S1I = S1I*CSR
        S2R = STR
        S2I = STI
        S1R = S1R*CSSR(KFLAG)
        S1I = S1I*CSSR(KFLAG)
        S2R = S2R*CSSR(KFLAG)
        S2I = S2I*CSSR(KFLAG)
        CSR = CSRR(KFLAG)
   80 CONTINUE
      RETURN
   90 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END


C-ENM added here ZBUNK
*DECK ZBUNK
      SUBROUTINE ZBUNK (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  ZBUNK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CBUNI-A, ZBUNI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL.
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
C     IN ZUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN ZUNK2
C
C***SEE ALSO  ZBESH, ZBESK
C***ROUTINES CALLED  ZUNK1, ZUNK2
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZBUNK
C     COMPLEX Y,Z
      DOUBLE PRECISION ALIM, AX, AY, ELIM, FNU, TOL, YI, YR, ZI, ZR
      INTEGER KODE, MR, N, NZ
      DIMENSION YR(N), YI(N)
C***FIRST EXECUTABLE STATEMENT  ZBUNK
      NZ = 0
      AX = ABS(ZR)*1.7321D0
      AY = ABS(ZI)
      IF (AY.GT.AX) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL ZUNK1(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
      GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL ZUNK2(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
   20 CONTINUE
      RETURN
      END

C-ENM added ZUNK1 and ZUNK2 here
*DECK ZUNK1
      SUBROUTINE ZUNK1 (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  ZUNK1
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUNK1-A, ZUNK1-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSION.
C     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***SEE ALSO  ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZS1S2, ZUCHK, ZUNIK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZUNK1
C     COMPLEX CFN,CK,CONE,CRSC,CS,CSCL,CSGN,CSPN,CSR,CSS,CWRK,CY,CZERO,
C    *C1,C2,PHI,PHID,RZ,SUM,SUMD,S1,S2,Y,Z,ZETA1,ZETA1D,ZETA2,ZETA2D,ZR
      DOUBLE PRECISION ALIM, ANG, APHI, ASC, ASCLE, BRY, CKI, CKR,
     * CONER, CRSC, CSCL, CSGNI, CSPNI, CSPNR, CSR, CSRR, CSSR,
     * CWRKI, CWRKR, CYI, CYR, C1I, C1R, C2I, C2M, C2R, ELIM, FMR, FN,
     * FNF, FNU, PHIDI, PHIDR, PHII, PHIR, PI, RAST, RAZR, RS1, RZI,
     * RZR, SGN, STI, STR, SUMDI, SUMDR, SUMI, SUMR, S1I, S1R, S2I,
     * S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R,
     * ZET1DI, ZET1DR, ZET2DI, ZET2DR, ZI, ZR, ZRI, ZRR, D1MACH, ZABS
      INTEGER I, IB, IFLAG, IFN, IL, INIT, INU, IUF, K, KDFLG, KFLAG,
     * KK, KODE, MR, N, NW, NZ, INITD, IC, IPARD, J, M
      DIMENSION BRY(3), INIT(2), YR(N), YI(N), SUMR(2), SUMI(2),
     * ZETA1R(2), ZETA1I(2), ZETA2R(2), ZETA2I(2), CYR(2), CYI(2),
     * CWRKR(16,3), CWRKI(16,3), CSSR(3), CSRR(3), PHIR(2), PHII(2)
      EXTERNAL ZABS
      DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
      DATA PI / 3.14159265358979324D0 /
C***FIRST EXECUTABLE STATEMENT  ZUNK1
      KDFLG = 1
      NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      ZRR = ZR
      ZRI = ZI
      IF (ZR.GE.0.0D0) GO TO 10
      ZRR = -ZR
      ZRI = -ZI
   10 CONTINUE
      J = 2
      DO 70 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
        J = 3 - J
        FN = FNU + (I-1)
        INIT(J) = 0
        CALL ZUNIK(ZRR, ZRI, FN, 2, 0, TOL, INIT(J), PHIR(J), PHII(J),
     *   ZETA1R(J), ZETA1I(J), ZETA2R(J), ZETA2I(J), SUMR(J), SUMI(J),
     *   CWRKR(1,J), CWRKI(1,J))
        IF (KODE.EQ.1) GO TO 20
        STR = ZRR + ZETA2R(J)
        STI = ZRI + ZETA2I(J)
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = ZETA1R(J) - STR
        S1I = ZETA1I(J) - STI
        GO TO 30
   20   CONTINUE
        S1R = ZETA1R(J) - ZETA2R(J)
        S1I = ZETA1I(J) - ZETA2I(J)
   30   CONTINUE
        RS1 = S1R
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        IF (ABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 40
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIR(J),PHII(J))
        RS1 = RS1 + LOG(APHI)
        IF (ABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 40
        IF (KDFLG.EQ.1) KFLAG = 3
   40   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        S2R = PHIR(J)*SUMR(J) - PHII(J)*SUMI(J)
        S2I = PHIR(J)*SUMI(J) + PHII(J)*SUMR(J)
        STR = EXP(S1R)*CSSR(KFLAG)
        S1R = STR*COS(S1I)
        S1I = STR*SIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S1R*S2I + S2R*S1I
        S2R = STR
        IF (KFLAG.NE.1) GO TO 50
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 60
   50   CONTINUE
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        YR(I) = S2R*CSRR(KFLAG)
        YI(I) = S2I*CSRR(KFLAG)
        IF (KDFLG.EQ.2) GO TO 75
        KDFLG = 2
        GO TO 70
   60   CONTINUE
        IF (RS1.GT.0.0D0) GO TO 300
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
        IF (ZR.LT.0.0D0) GO TO 300
        KDFLG = 1
        YR(I)=ZEROR
        YI(I)=ZEROI
        NZ=NZ+1
        IF (I.EQ.1) GO TO 70
        IF ((YR(I-1).EQ.ZEROR).AND.(YI(I-1).EQ.ZEROI)) GO TO 70
        YR(I-1)=ZEROR
        YI(I-1)=ZEROI
        NZ=NZ+1
   70 CONTINUE
      I = N
   75 CONTINUE
      RAZR = 1.0D0/ZABS(ZRR,ZRI)
      STR = ZRR*RAZR
      STI = -ZRI*RAZR
      RZR = (STR+STR)*RAZR
      RZI = (STI+STI)*RAZR
      CKR = FN*RZR
      CKI = FN*RZI
      IB = I + 1
      IF (N.LT.IB) GO TO 160
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
C     ON UNDERFLOW.
C-----------------------------------------------------------------------
      FN = FNU + (N-1)
      IPARD = 1
      IF (MR.NE.0) IPARD = 0
      INITD = 0
      CALL ZUNIK(ZRR, ZRI, FN, 2, IPARD, TOL, INITD, PHIDR, PHIDI,
     * ZET1DR, ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI, CWRKR(1,3),
     * CWRKI(1,3))
      IF (KODE.EQ.1) GO TO 80
      STR = ZRR + ZET2DR
      STI = ZRI + ZET2DI
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = ZET1DR - STR
      S1I = ZET1DI - STI
      GO TO 90
   80 CONTINUE
      S1R = ZET1DR - ZET2DR
      S1I = ZET1DI - ZET2DI
   90 CONTINUE
      RS1 = S1R
      IF (ABS(RS1).GT.ELIM) GO TO 95
      IF (ABS(RS1).LT.ALIM) GO TO 100
C-----------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-----------------------------------------------------------------------
      APHI = ZABS(PHIDR,PHIDI)
      RS1 = RS1+LOG(APHI)
      IF (ABS(RS1).LT.ELIM) GO TO 100
   95 CONTINUE
      IF (ABS(RS1).GT.0.0D0) GO TO 300
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
      IF (ZR.LT.0.0D0) GO TO 300
      NZ = N
      DO 96 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
   96 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     FORWARD RECUR FOR REMAINDER OF THE SEQUENCE
C-----------------------------------------------------------------------
  100 CONTINUE
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 120 I=IB,N
        C2R = S2R
        C2I = S2I
        S2R = CKR*C2R - CKI*C2I + S1R
        S2I = CKR*C2I + CKI*C2R + S1I
        S1R = C2R
        S1I = C2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(I) = C2R
        YI(I) = C2I
        IF (KFLAG.GE.3) GO TO 120
        STR = ABS(C2R)
        STI = ABS(C2I)
        C2M = MAX(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 120
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(KFLAG)
        S1I = S1I*CSSR(KFLAG)
        S2R = S2R*CSSR(KFLAG)
        S2I = S2I*CSSR(KFLAG)
        C1R = CSRR(KFLAG)
  120 CONTINUE
  160 CONTINUE
      IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0
C-----------------------------------------------------------------------
      NZ = 0
      FMR = MR
      SGN = -DSIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
C-----------------------------------------------------------------------
      CSGNI = SGN
      INU = FNU
      FNF = FNU - INU
      IFN = INU + N - 1
      ANG = FNF*SGN
      CSPNR = COS(ANG)
      CSPNI = SIN(ANG)
      IF (MOD(IFN,2).EQ.0) GO TO 170
      CSPNR = -CSPNR
      CSPNI = -CSPNI
  170 CONTINUE
      ASC = BRY(1)
      IUF = 0
      KK = N
      KDFLG = 1
      IB = IB - 1
      IC = IB - 1
      DO 270 K=1,N
        FN = FNU + (KK-1)
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
        M=3
        IF (N.GT.2) GO TO 175
  172   CONTINUE
        INITD = INIT(J)
        PHIDR = PHIR(J)
        PHIDI = PHII(J)
        ZET1DR = ZETA1R(J)
        ZET1DI = ZETA1I(J)
        ZET2DR = ZETA2R(J)
        ZET2DI = ZETA2I(J)
        SUMDR = SUMR(J)
        SUMDI = SUMI(J)
        M = J
        J = 3 - J
        GO TO 180
  175   CONTINUE
        IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 180
        IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 172
        INITD = 0
  180   CONTINUE
        CALL ZUNIK(ZRR, ZRI, FN, 1, 0, TOL, INITD, PHIDR, PHIDI,
     *   ZET1DR, ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI,
     *   CWRKR(1,M), CWRKI(1,M))
        IF (KODE.EQ.1) GO TO 200
        STR = ZRR + ZET2DR
        STI = ZRI + ZET2DI
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZET1DR + STR
        S1I = -ZET1DI + STI
        GO TO 210
  200   CONTINUE
        S1R = -ZET1DR + ZET2DR
        S1I = -ZET1DI + ZET2DI
  210   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (ABS(RS1).GT.ELIM) GO TO 260
        IF (KDFLG.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 220
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIDR,PHIDI)
        RS1 = RS1 + LOG(APHI)
        IF (ABS(RS1).GT.ELIM) GO TO 260
        IF (KDFLG.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 220
        IF (KDFLG.EQ.1) IFLAG = 3
  220   CONTINUE
        STR = PHIDR*SUMDR - PHIDI*SUMDI
        STI = PHIDR*SUMDI + PHIDI*SUMDR
        S2R = -CSGNI*STI
        S2I = CSGNI*STR
        STR = EXP(S1R)*CSSR(IFLAG)
        S1R = STR*COS(S1I)
        S1I = STR*SIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 230
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.EQ.0) GO TO 230
        S2R = ZEROR
        S2I = ZEROI
  230   CONTINUE
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        C2R = S2R
        C2I = S2I
        S2R = S2R*CSRR(IFLAG)
        S2I = S2I*CSRR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
        S1R = YR(KK)
        S1I = YI(KK)
        IF (KODE.EQ.1) GO TO 250
        CALL ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  250   CONTINUE
        YR(KK) = S1R*CSPNR - S1I*CSPNI + S2R
        YI(KK) = CSPNR*S1I + CSPNI*S1R + S2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        IF (C2R.NE.0.0D0 .OR. C2I.NE.0.0D0) GO TO 255
        KDFLG = 1
        GO TO 270
  255   CONTINUE
        IF (KDFLG.EQ.2) GO TO 275
        KDFLG = 2
        GO TO 270
  260   CONTINUE
        IF (RS1.GT.0.0D0) GO TO 300
        S2R = ZEROR
        S2I = ZEROI
        GO TO 230
  270 CONTINUE
      K = N
  275 CONTINUE
      IL = N - K
      IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      CSR = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = INU+IL
      DO 290 I=1,IL
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FN+FNF)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FN+FNF)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        FN = FN - 1.0D0
        C2R = S2R*CSR
        C2I = S2I*CSR
        CKR = C2R
        CKI = C2I
        C1R = YR(KK)
        C1I = YI(KK)
        IF (KODE.EQ.1) GO TO 280
        CALL ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  280   CONTINUE
        YR(KK) = C1R*CSPNR - C1I*CSPNI + C2R
        YI(KK) = C1R*CSPNI + C1I*CSPNR + C2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        IF (IFLAG.GE.3) GO TO 290
        C2R = ABS(CKR)
        C2I = ABS(CKI)
        C2M = MAX(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 290
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSR
        S1I = S1I*CSR
        S2R = CKR
        S2I = CKI
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        CSR = CSRR(IFLAG)
  290 CONTINUE
      RETURN
  300 CONTINUE
      NZ = -1
      RETURN
      END
*DECK ZUNK2
      SUBROUTINE ZUNK2 (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  ZUNK2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUNK2-A, ZUNK2-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
C     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
C     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT
C     HALF PLANE OR ZR=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC-
C     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***SEE ALSO  ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZAIRY, ZS1S2, ZUCHK, ZUNHJ
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZUNK2
C     COMPLEX AI,ARG,ARGD,ASUM,ASUMD,BSUM,BSUMD,CFN,CI,CIP,CK,CONE,CRSC,
C    *CR1,CR2,CS,CSCL,CSGN,CSPN,CSR,CSS,CY,CZERO,C1,C2,DAI,PHI,PHID,RZ,
C    *S1,S2,Y,Z,ZB,ZETA1,ZETA1D,ZETA2,ZETA2D,ZN,ZR
      DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGDI,
     * ARGDR, ARGI, ARGR, ASC, ASCLE, ASUMDI, ASUMDR, ASUMI, ASUMR,
     * BRY, BSUMDI, BSUMDR, BSUMI, BSUMR, CAR, CIPI, CIPR, CKI, CKR,
     * CONER, CRSC, CR1I, CR1R, CR2I, CR2R, CSCL, CSGNI, CSI,
     * CSPNI, CSPNR, CSR, CSRR, CSSR, CYI, CYR, C1I, C1R, C2I, C2M,
     * C2R, DAII, DAIR, ELIM, FMR, FN, FNF, FNU, HPI, PHIDI, PHIDR,
     * PHII, PHIR, PI, PTI, PTR, RAST, RAZR, RS1, RZI, RZR, SAR, SGN,
     * STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, YY, ZBI, ZBR, ZEROI,
     * ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZET1DI, ZET1DR, ZET2DI,
     * ZET2DR, ZI, ZNI, ZNR, ZR, ZRI, ZRR, D1MACH, ZABS
      INTEGER I, IB, IFLAG, IFN, IL, IN, INU, IUF, K, KDFLG, KFLAG, KK,
     * KODE, MR, N, NAI, NDAI, NW, NZ, IDUM, J, IPARD, IC
      DIMENSION BRY(3), YR(N), YI(N), ASUMR(2), ASUMI(2), BSUMR(2),
     * BSUMI(2), PHIR(2), PHII(2), ARGR(2), ARGI(2), ZETA1R(2),
     * ZETA1I(2), ZETA2R(2), ZETA2I(2), CYR(2), CYI(2), CIPR(4),
     * CIPI(4), CSSR(3), CSRR(3)
      EXTERNAL ZABS
      DATA ZEROR,ZEROI,CONER,CR1R,CR1I,CR2R,CR2I /
     1         0.0D0, 0.0D0, 1.0D0,
     1 1.0D0,1.73205080756887729D0 , -0.5D0,-8.66025403784438647D-01 /
      DATA HPI, PI, AIC /
     1     1.57079632679489662D+00,     3.14159265358979324D+00,
     1     1.26551212348464539D+00/
      DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4),
     * CIPI(4) /
     1  1.0D0,0.0D0 ,  0.0D0,-1.0D0 ,  -1.0D0,0.0D0 ,  0.0D0,1.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZUNK2
      KDFLG = 1
      NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      ZRR = ZR
      ZRI = ZI
      IF (ZR.GE.0.0D0) GO TO 10
      ZRR = -ZR
      ZRI = -ZI
   10 CONTINUE
      YY = ZRI
      ZNR = ZRI
      ZNI = -ZRR
      ZBR = ZRR
      ZBI = ZRI
      INU = FNU
      FNF = FNU - INU
      ANG = -HPI*FNF
      CAR = COS(ANG)
      SAR = SIN(ANG)
      C2R = HPI*SAR
      C2I = -HPI*CAR
      KK = MOD(INU,4) + 1
      STR = C2R*CIPR(KK) - C2I*CIPI(KK)
      STI = C2R*CIPI(KK) + C2I*CIPR(KK)
      CSR = CR1R*STR - CR1I*STI
      CSI = CR1R*STI + CR1I*STR
      IF (YY.GT.0.0D0) GO TO 20
      ZNR = -ZNR
      ZBI = -ZBI
   20 CONTINUE
C-----------------------------------------------------------------------
C     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
      J = 2
      DO 80 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
        J = 3 - J
        FN = FNU + (I-1)
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR(J), PHII(J), ARGR(J),
     *   ARGI(J), ZETA1R(J), ZETA1I(J), ZETA2R(J), ZETA2I(J), ASUMR(J),
     *   ASUMI(J), BSUMR(J), BSUMI(J))
        IF (KODE.EQ.1) GO TO 30
        STR = ZBR + ZETA2R(J)
        STI = ZBI + ZETA2I(J)
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = ZETA1R(J) - STR
        S1I = ZETA1I(J) - STI
        GO TO 40
   30   CONTINUE
        S1R = ZETA1R(J) - ZETA2R(J)
        S1I = ZETA1I(J) - ZETA2I(J)
   40   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (ABS(RS1).GT.ELIM) GO TO 70
        IF (KDFLG.EQ.1) KFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 50
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIR(J),PHII(J))
        AARG = ZABS(ARGR(J),ARGI(J))
        RS1 = RS1 + LOG(APHI) - 0.25D0*LOG(AARG) - AIC
        IF (ABS(RS1).GT.ELIM) GO TO 70
        IF (KDFLG.EQ.1) KFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 50
        IF (KDFLG.EQ.1) KFLAG = 3
   50   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        C2R = ARGR(J)*CR2R - ARGI(J)*CR2I
        C2I = ARGR(J)*CR2I + ARGI(J)*CR2R
        CALL ZAIRY(C2R, C2I, 0, 2, AIR, AII, NAI, IDUM)
        CALL ZAIRY(C2R, C2I, 1, 2, DAIR, DAII, NDAI, IDUM)
        STR = DAIR*BSUMR(J) - DAII*BSUMI(J)
        STI = DAIR*BSUMI(J) + DAII*BSUMR(J)
        PTR = STR*CR2R - STI*CR2I
        PTI = STR*CR2I + STI*CR2R
        STR = PTR + (AIR*ASUMR(J)-AII*ASUMI(J))
        STI = PTI + (AIR*ASUMI(J)+AII*ASUMR(J))
        PTR = STR*PHIR(J) - STI*PHII(J)
        PTI = STR*PHII(J) + STI*PHIR(J)
        S2R = PTR*CSR - PTI*CSI
        S2I = PTR*CSI + PTI*CSR
        STR = EXP(S1R)*CSSR(KFLAG)
        S1R = STR*COS(S1I)
        S1I = STR*SIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S1R*S2I + S2R*S1I
        S2R = STR
        IF (KFLAG.NE.1) GO TO 60
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 70
   60   CONTINUE
        IF (YY.LE.0.0D0) S2I = -S2I
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        YR(I) = S2R*CSRR(KFLAG)
        YI(I) = S2I*CSRR(KFLAG)
        STR = CSI
        CSI = -CSR
        CSR = STR
        IF (KDFLG.EQ.2) GO TO 85
        KDFLG = 2
        GO TO 80
   70   CONTINUE
        IF (RS1.GT.0.0D0) GO TO 320
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
        IF (ZR.LT.0.0D0) GO TO 320
        KDFLG = 1
        YR(I)=ZEROR
        YI(I)=ZEROI
        NZ=NZ+1
        STR = CSI
        CSI =-CSR
        CSR = STR
        IF (I.EQ.1) GO TO 80
        IF ((YR(I-1).EQ.ZEROR).AND.(YI(I-1).EQ.ZEROI)) GO TO 80
        YR(I-1)=ZEROR
        YI(I-1)=ZEROI
        NZ=NZ+1
   80 CONTINUE
      I = N
   85 CONTINUE
      RAZR = 1.0D0/ZABS(ZRR,ZRI)
      STR = ZRR*RAZR
      STI = -ZRI*RAZR
      RZR = (STR+STR)*RAZR
      RZI = (STI+STI)*RAZR
      CKR = FN*RZR
      CKI = FN*RZI
      IB = I + 1
      IF (N.LT.IB) GO TO 180
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
C     ON UNDERFLOW.
C-----------------------------------------------------------------------
      FN = FNU + (N-1)
      IPARD = 1
      IF (MR.NE.0) IPARD = 0
      CALL ZUNHJ(ZNR, ZNI, FN, IPARD, TOL, PHIDR, PHIDI, ARGDR, ARGDI,
     * ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR, ASUMDI, BSUMDR, BSUMDI)
      IF (KODE.EQ.1) GO TO 90
      STR = ZBR + ZET2DR
      STI = ZBI + ZET2DI
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = ZET1DR - STR
      S1I = ZET1DI - STI
      GO TO 100
   90 CONTINUE
      S1R = ZET1DR - ZET2DR
      S1I = ZET1DI - ZET2DI
  100 CONTINUE
      RS1 = S1R
      IF (ABS(RS1).GT.ELIM) GO TO 105
      IF (ABS(RS1).LT.ALIM) GO TO 120
C-----------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-----------------------------------------------------------------------
      APHI = ZABS(PHIDR,PHIDI)
      RS1 = RS1+LOG(APHI)
      IF (ABS(RS1).LT.ELIM) GO TO 120
  105 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 320
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
      IF (ZR.LT.0.0D0) GO TO 320
      NZ = N
      DO 106 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  106 CONTINUE
      RETURN
  120 CONTINUE
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 130 I=IB,N
        C2R = S2R
        C2I = S2I
        S2R = CKR*C2R - CKI*C2I + S1R
        S2I = CKR*C2I + CKI*C2R + S1I
        S1R = C2R
        S1I = C2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(I) = C2R
        YI(I) = C2I
        IF (KFLAG.GE.3) GO TO 130
        STR = ABS(C2R)
        STI = ABS(C2I)
        C2M = MAX(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 130
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(KFLAG)
        S1I = S1I*CSSR(KFLAG)
        S2R = S2R*CSSR(KFLAG)
        S2I = S2I*CSSR(KFLAG)
        C1R = CSRR(KFLAG)
  130 CONTINUE
  180 CONTINUE
      IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0
C-----------------------------------------------------------------------
      NZ = 0
      FMR = MR
      SGN = -DSIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
C-----------------------------------------------------------------------
      CSGNI = SGN
      IF (YY.LE.0.0D0) CSGNI = -CSGNI
      IFN = INU + N - 1
      ANG = FNF*SGN
      CSPNR = COS(ANG)
      CSPNI = SIN(ANG)
      IF (MOD(IFN,2).EQ.0) GO TO 190
      CSPNR = -CSPNR
      CSPNI = -CSPNI
  190 CONTINUE
C-----------------------------------------------------------------------
C     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
C     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
      CSR = SAR*CSGNI
      CSI = CAR*CSGNI
      IN = MOD(IFN,4) + 1
      C2R = CIPR(IN)
      C2I = CIPI(IN)
      STR = CSR*C2R + CSI*C2I
      CSI = -CSR*C2I + CSI*C2R
      CSR = STR
      ASC = BRY(1)
      IUF = 0
      KK = N
      KDFLG = 1
      IB = IB - 1
      IC = IB - 1
      DO 290 K=1,N
        FN = FNU + (KK-1)
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
        IF (N.GT.2) GO TO 175
  172   CONTINUE
        PHIDR = PHIR(J)
        PHIDI = PHII(J)
        ARGDR = ARGR(J)
        ARGDI = ARGI(J)
        ZET1DR = ZETA1R(J)
        ZET1DI = ZETA1I(J)
        ZET2DR = ZETA2R(J)
        ZET2DI = ZETA2I(J)
        ASUMDR = ASUMR(J)
        ASUMDI = ASUMI(J)
        BSUMDR = BSUMR(J)
        BSUMDI = BSUMI(J)
        J = 3 - J
        GO TO 210
  175   CONTINUE
        IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 210
        IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 172
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIDR, PHIDI, ARGDR,
     *   ARGDI, ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR,
     *   ASUMDI, BSUMDR, BSUMDI)
  210   CONTINUE
        IF (KODE.EQ.1) GO TO 220
        STR = ZBR + ZET2DR
        STI = ZBI + ZET2DI
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZET1DR + STR
        S1I = -ZET1DI + STI
        GO TO 230
  220   CONTINUE
        S1R = -ZET1DR + ZET2DR
        S1I = -ZET1DI + ZET2DI
  230   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (ABS(RS1).GT.ELIM) GO TO 280
        IF (KDFLG.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 240
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIDR,PHIDI)
        AARG = ZABS(ARGDR,ARGDI)
        RS1 = RS1 + LOG(APHI) - 0.25D0*LOG(AARG) - AIC
        IF (ABS(RS1).GT.ELIM) GO TO 280
        IF (KDFLG.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 240
        IF (KDFLG.EQ.1) IFLAG = 3
  240   CONTINUE
        CALL ZAIRY(ARGDR, ARGDI, 0, 2, AIR, AII, NAI, IDUM)
        CALL ZAIRY(ARGDR, ARGDI, 1, 2, DAIR, DAII, NDAI, IDUM)
        STR = DAIR*BSUMDR - DAII*BSUMDI
        STI = DAIR*BSUMDI + DAII*BSUMDR
        STR = STR + (AIR*ASUMDR-AII*ASUMDI)
        STI = STI + (AIR*ASUMDI+AII*ASUMDR)
        PTR = STR*PHIDR - STI*PHIDI
        PTI = STR*PHIDI + STI*PHIDR
        S2R = PTR*CSR - PTI*CSI
        S2I = PTR*CSI + PTI*CSR
        STR = EXP(S1R)*CSSR(IFLAG)
        S1R = STR*COS(S1I)
        S1I = STR*SIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 250
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.EQ.0) GO TO 250
        S2R = ZEROR
        S2I = ZEROI
  250   CONTINUE
        IF (YY.LE.0.0D0) S2I = -S2I
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        C2R = S2R
        C2I = S2I
        S2R = S2R*CSRR(IFLAG)
        S2I = S2I*CSRR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
        S1R = YR(KK)
        S1I = YI(KK)
        IF (KODE.EQ.1) GO TO 270
        CALL ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  270   CONTINUE
        YR(KK) = S1R*CSPNR - S1I*CSPNI + S2R
        YI(KK) = S1R*CSPNI + S1I*CSPNR + S2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        STR = CSI
        CSI = -CSR
        CSR = STR
        IF (C2R.NE.0.0D0 .OR. C2I.NE.0.0D0) GO TO 255
        KDFLG = 1
        GO TO 290
  255   CONTINUE
        IF (KDFLG.EQ.2) GO TO 295
        KDFLG = 2
        GO TO 290
  280   CONTINUE
        IF (RS1.GT.0.0D0) GO TO 320
        S2R = ZEROR
        S2I = ZEROI
        GO TO 250
  290 CONTINUE
      K = N
  295 CONTINUE
      IL = N - K
      IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      CSR = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = INU+IL
      DO 310 I=1,IL
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FN+FNF)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FN+FNF)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        FN = FN - 1.0D0
        C2R = S2R*CSR
        C2I = S2I*CSR
        CKR = C2R
        CKI = C2I
        C1R = YR(KK)
        C1I = YI(KK)
        IF (KODE.EQ.1) GO TO 300
        CALL ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  300   CONTINUE
        YR(KK) = C1R*CSPNR - C1I*CSPNI + C2R
        YI(KK) = C1R*CSPNI + C1I*CSPNR + C2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        IF (IFLAG.GE.3) GO TO 310
        C2R = ABS(CKR)
        C2I = ABS(CKI)
        C2M = MAX(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 310
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSR
        S1I = S1I*CSR
        S2R = CKR
        S2I = CKI
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        CSR = CSRR(IFLAG)
  310 CONTINUE
      RETURN
  320 CONTINUE
      NZ = -1
      RETURN
      END


*DECK ZMLT
      SUBROUTINE ZMLT (AR, AI, BR, BI, CR, CI)
C***BEGIN PROLOGUE  ZMLT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZMLT-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B.
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZMLT
      DOUBLE PRECISION AR, AI, BR, BI, CR, CI, CA, CB
C***FIRST EXECUTABLE STATEMENT  ZMLT
      CA = AR*BR - AI*BI
      CB = AR*BI + AI*BR
      CR = CA
      CI = CB
      RETURN
      END

*DECK ZUCHK
      SUBROUTINE ZUCHK (YR, YI, NZ, ASCLE, TOL)
C***BEGIN PROLOGUE  ZUCHK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SERI, ZUOIK, ZUNK1, ZUNK2, ZUNI1, ZUNI2 and
C            ZKSCL
C***LIBRARY   SLATEC
C***TYPE      ALL (CUCHK-A, ZUCHK-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
C      EXP(-ALIM)=ASCLE=1.0E+3*D1MACH(1)/TOL. THE TEST IS MADE TO SEE
C      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
C      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
C      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
C      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
C      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
C
C***SEE ALSO  SERI, ZKSCL, ZUNI1, ZUNI2, ZUNK1, ZUNK2, ZUOIK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZUCHK
C
C     COMPLEX Y
      DOUBLE PRECISION ASCLE, SS, ST, TOL, WR, WI, YR, YI
      INTEGER NZ
C***FIRST EXECUTABLE STATEMENT  ZUCHK
      NZ = 0
      WR = ABS(YR)
      WI = ABS(YI)
      ST = MIN(WR,WI)
      IF (ST.GT.ASCLE) RETURN
      SS = MAX(WR,WI)
      ST = ST/TOL
      IF (SS.LT.ST) NZ = 1
      RETURN
      END


*DECK ZLOG
      SUBROUTINE ZLOG (AR, AI, BR, BI, IERR)
C***BEGIN PROLOGUE  ZLOG
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZLOG-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     DOUBLE PRECISION COMPLEX LOGARITHM B=CLOG(A)
C     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0)
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  ZABS
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZLOG
      DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DHPI
      DOUBLE PRECISION ZABS
      INTEGER IERR
      EXTERNAL ZABS
      DATA DPI , DHPI  / 3.141592653589793238462643383D+0,
     1                   1.570796326794896619231321696D+0/
C***FIRST EXECUTABLE STATEMENT  ZLOG
      IERR=0
      IF (AR.EQ.0.0D+0) GO TO 10
      IF (AI.EQ.0.0D+0) GO TO 20
      DTHETA = DATAN(AI/AR)
      IF (DTHETA.LE.0.0D+0) GO TO 40
      IF (AR.LT.0.0D+0) DTHETA = DTHETA - DPI
      GO TO 50
   10 IF (AI.EQ.0.0D+0) GO TO 60
      BI = DHPI
      BR = LOG(ABS(AI))
      IF (AI.LT.0.0D+0) BI = -BI
      RETURN
   20 IF (AR.GT.0.0D+0) GO TO 30
      BR = LOG(ABS(AR))
      BI = DPI
      RETURN
   30 BR = LOG(AR)
      BI = 0.0D+0
      RETURN
   40 IF (AR.LT.0.0D+0) DTHETA = DTHETA + DPI
   50 ZM = ZABS(AR,AI)
      BR = LOG(ZM)
      BI = DTHETA
      RETURN
   60 CONTINUE
      IERR=1
      RETURN
      END

*DECK DGAMLN
      DOUBLE PRECISION FUNCTION DGAMLN (Z, IERR)
C***BEGIN PROLOGUE  DGAMLN
C***SUBSIDIARY
C***PURPOSE  Compute the logarithm of the Gamma function
C***LIBRARY   SLATEC
C***CATEGORY  C7A
C***TYPE      DOUBLE PRECISION (GAMLN-S, DGAMLN-D)
C***KEYWORDS  LOGARITHM OF GAMMA FUNCTION
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C               **** A DOUBLE PRECISION ROUTINE ****
C         DGAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
C         Z.GT.0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
C         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
C         G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN.  THE FUNCTION WAS MADE AS
C         PORTABLE AS POSSIBLE BY COMPUTING ZMIN FROM THE NUMBER OF BASE
C         10 DIGITS IN A WORD, RLN=MAX(-ALOG10(R1MACH(4)),0.5E-18)
C         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.
C
C         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
C         VALUES IS USED FOR SPEED OF EXECUTION.
C
C     DESCRIPTION OF ARGUMENTS
C
C         INPUT      Z IS D0UBLE PRECISION
C           Z      - ARGUMENT, Z.GT.0.0D0
C
C         OUTPUT      DGAMLN IS DOUBLE PRECISION
C           DGAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z.NE.0.0D0
C           IERR    - ERROR FLAG
C                     IERR=0, NORMAL RETURN, COMPUTATION COMPLETED
C                     IERR=1, Z.LE.0.0D0,    NO COMPUTATION
C
C
C***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C***ROUTINES CALLED  D1MACH, I1MACH
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   830501  REVISION DATE from Version 3.2
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   920128  Category corrected.  (WRB)
C   921215  DGAMLN defined for Z negative.  (WRB)
C***END PROLOGUE  DGAMLN
      DOUBLE PRECISION CF, CON, FLN, FZ, GLN, RLN, S, TLG, TRM, TST,
     * T1, WDTOL, Z, ZDMY, ZINC, ZM, ZMIN, ZP, ZSQ, D1MACH
      INTEGER I, IERR, I1M, K, MZ, NZ, I1MACH
      DIMENSION CF(22), GLN(100)
C           LNGAMMA(N), N=1,100
      DATA GLN(1), GLN(2), GLN(3), GLN(4), GLN(5), GLN(6), GLN(7),
     1     GLN(8), GLN(9), GLN(10), GLN(11), GLN(12), GLN(13), GLN(14),
     2     GLN(15), GLN(16), GLN(17), GLN(18), GLN(19), GLN(20),
     3     GLN(21), GLN(22)/
     4     0.00000000000000000D+00,     0.00000000000000000D+00,
     5     6.93147180559945309D-01,     1.79175946922805500D+00,
     6     3.17805383034794562D+00,     4.78749174278204599D+00,
     7     6.57925121201010100D+00,     8.52516136106541430D+00,
     8     1.06046029027452502D+01,     1.28018274800814696D+01,
     9     1.51044125730755153D+01,     1.75023078458738858D+01,
     A     1.99872144956618861D+01,     2.25521638531234229D+01,
     B     2.51912211827386815D+01,     2.78992713838408916D+01,
     C     3.06718601060806728D+01,     3.35050734501368889D+01,
     D     3.63954452080330536D+01,     3.93398841871994940D+01,
     E     4.23356164607534850D+01,     4.53801388984769080D+01/
      DATA GLN(23), GLN(24), GLN(25), GLN(26), GLN(27), GLN(28),
     1     GLN(29), GLN(30), GLN(31), GLN(32), GLN(33), GLN(34),
     2     GLN(35), GLN(36), GLN(37), GLN(38), GLN(39), GLN(40),
     3     GLN(41), GLN(42), GLN(43), GLN(44)/
     4     4.84711813518352239D+01,     5.16066755677643736D+01,
     5     5.47847293981123192D+01,     5.80036052229805199D+01,
     6     6.12617017610020020D+01,     6.45575386270063311D+01,
     7     6.78897431371815350D+01,     7.12570389671680090D+01,
     8     7.46582363488301644D+01,     7.80922235533153106D+01,
     9     8.15579594561150372D+01,     8.50544670175815174D+01,
     A     8.85808275421976788D+01,     9.21361756036870925D+01,
     B     9.57196945421432025D+01,     9.93306124547874269D+01,
     C     1.02968198614513813D+02,     1.06631760260643459D+02,
     D     1.10320639714757395D+02,     1.14034211781461703D+02,
     E     1.17771881399745072D+02,     1.21533081515438634D+02/
      DATA GLN(45), GLN(46), GLN(47), GLN(48), GLN(49), GLN(50),
     1     GLN(51), GLN(52), GLN(53), GLN(54), GLN(55), GLN(56),
     2     GLN(57), GLN(58), GLN(59), GLN(60), GLN(61), GLN(62),
     3     GLN(63), GLN(64), GLN(65), GLN(66)/
     4     1.25317271149356895D+02,     1.29123933639127215D+02,
     5     1.32952575035616310D+02,     1.36802722637326368D+02,
     6     1.40673923648234259D+02,     1.44565743946344886D+02,
     7     1.48477766951773032D+02,     1.52409592584497358D+02,
     8     1.56360836303078785D+02,     1.60331128216630907D+02,
     9     1.64320112263195181D+02,     1.68327445448427652D+02,
     A     1.72352797139162802D+02,     1.76395848406997352D+02,
     B     1.80456291417543771D+02,     1.84533828861449491D+02,
     C     1.88628173423671591D+02,     1.92739047287844902D+02,
     D     1.96866181672889994D+02,     2.01009316399281527D+02,
     E     2.05168199482641199D+02,     2.09342586752536836D+02/
      DATA GLN(67), GLN(68), GLN(69), GLN(70), GLN(71), GLN(72),
     1     GLN(73), GLN(74), GLN(75), GLN(76), GLN(77), GLN(78),
     2     GLN(79), GLN(80), GLN(81), GLN(82), GLN(83), GLN(84),
     3     GLN(85), GLN(86), GLN(87), GLN(88)/
     4     2.13532241494563261D+02,     2.17736934113954227D+02,
     5     2.21956441819130334D+02,     2.26190548323727593D+02,
     6     2.30439043565776952D+02,     2.34701723442818268D+02,
     7     2.38978389561834323D+02,     2.43268849002982714D+02,
     8     2.47572914096186884D+02,     2.51890402209723194D+02,
     9     2.56221135550009525D+02,     2.60564940971863209D+02,
     A     2.64921649798552801D+02,     2.69291097651019823D+02,
     B     2.73673124285693704D+02,     2.78067573440366143D+02,
     C     2.82474292687630396D+02,     2.86893133295426994D+02,
     D     2.91323950094270308D+02,     2.95766601350760624D+02,
     E     3.00220948647014132D+02,     3.04686856765668715D+02/
      DATA GLN(89), GLN(90), GLN(91), GLN(92), GLN(93), GLN(94),
     1     GLN(95), GLN(96), GLN(97), GLN(98), GLN(99), GLN(100)/
     2     3.09164193580146922D+02,     3.13652829949879062D+02,
     3     3.18152639620209327D+02,     3.22663499126726177D+02,
     4     3.27185287703775217D+02,     3.31717887196928473D+02,
     5     3.36261181979198477D+02,     3.40815058870799018D+02,
     6     3.45379407062266854D+02,     3.49954118040770237D+02,
     7     3.54539085519440809D+02,     3.59134205369575399D+02/
C             COEFFICIENTS OF ASYMPTOTIC EXPANSION
      DATA CF(1), CF(2), CF(3), CF(4), CF(5), CF(6), CF(7), CF(8),
     1     CF(9), CF(10), CF(11), CF(12), CF(13), CF(14), CF(15),
     2     CF(16), CF(17), CF(18), CF(19), CF(20), CF(21), CF(22)/
     3     8.33333333333333333D-02,    -2.77777777777777778D-03,
     4     7.93650793650793651D-04,    -5.95238095238095238D-04,
     5     8.41750841750841751D-04,    -1.91752691752691753D-03,
     6     6.41025641025641026D-03,    -2.95506535947712418D-02,
     7     1.79644372368830573D-01,    -1.39243221690590112D+00,
     8     1.34028640441683920D+01,    -1.56848284626002017D+02,
     9     2.19310333333333333D+03,    -3.61087712537249894D+04,
     A     6.91472268851313067D+05,    -1.52382215394074162D+07,
     B     3.82900751391414141D+08,    -1.08822660357843911D+10,
     C     3.47320283765002252D+11,    -1.23696021422692745D+13,
     D     4.88788064793079335D+14,    -2.13203339609193739D+16/
C
C             LN(2*PI)
      DATA CON                    /     1.83787706640934548D+00/
C
C***FIRST EXECUTABLE STATEMENT  DGAMLN
      IERR=0
      IF (Z.LE.0.0D0) GO TO 70
      IF (Z.GT.101.0D0) GO TO 10
      NZ = Z
      FZ = Z - NZ
      IF (FZ.GT.0.0D0) GO TO 10
      IF (NZ.GT.100) GO TO 10

cSAP091104
c      write(*,*)'besse_slatec.f in FUNCTION DGAMLN Z,NZ',Z,NZ

      DGAMLN = GLN(NZ)
      RETURN
   10 CONTINUE
      WDTOL = D1MACH(4)
      WDTOL = MAX(WDTOL,0.5D-18)
      I1M = I1MACH(14)
      RLN = D1MACH(5)*I1M
      FLN = MIN(RLN,20.0D0)
      FLN = MAX(FLN,3.0D0)
      FLN = FLN - 3.0D0
      ZM = 1.8000D0 + 0.3875D0*FLN
      MZ = ZM + 1
      ZMIN = MZ
      ZDMY = Z
      ZINC = 0.0D0
      IF (Z.GE.ZMIN) GO TO 20
      ZINC = ZMIN - NZ
      ZDMY = Z + ZINC
   20 CONTINUE
      ZP = 1.0D0/ZDMY
      T1 = CF(1)*ZP
      S = T1
      IF (ZP.LT.WDTOL) GO TO 40
      ZSQ = ZP*ZP
      TST = T1*WDTOL
      DO 30 K=2,22
        ZP = ZP*ZSQ
        TRM = CF(K)*ZP
        IF (ABS(TRM).LT.TST) GO TO 40
        S = S + TRM
   30 CONTINUE
   40 CONTINUE
      IF (ZINC.NE.0.0D0) GO TO 50
      TLG = LOG(Z)
      DGAMLN = Z*(TLG-1.0D0) + 0.5D0*(CON-TLG) + S
      RETURN
   50 CONTINUE
      ZP = 1.0D0
      NZ = ZINC
      DO 60 I=1,NZ
        ZP = ZP*(Z+(I-1))
   60 CONTINUE
      TLG = LOG(ZDMY)
      DGAMLN = ZDMY*(TLG-1.0D0) - LOG(ZP) + 0.5D0*(CON-TLG) + S
      RETURN
C
C
   70 CONTINUE
      DGAMLN = D1MACH(2)
      IERR=1
      RETURN
      END


*DECK ZDIV
      SUBROUTINE ZDIV (AR, AI, BR, BI, CR, CI)
C***BEGIN PROLOGUE  ZDIV
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZDIV-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     DOUBLE PRECISION COMPLEX DIVIDE C=A/B.
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  ZABS
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZDIV
      DOUBLE PRECISION AR, AI, BR, BI, CR, CI, BM, CA, CB, CC, CD
      DOUBLE PRECISION ZABS
      EXTERNAL ZABS
C***FIRST EXECUTABLE STATEMENT  ZDIV
      BM = 1.0D0/ZABS(BR,BI)
      CC = BR*BM
      CD = BI*BM
      CA = (AR*CC+AI*CD)*BM
      CB = (AI*CC-AR*CD)*BM
      CR = CA
      CI = CB
      RETURN
      END
*DECK ZEXP
      SUBROUTINE ZEXP (AR, AI, BR, BI)
C***BEGIN PROLOGUE  ZEXP
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZEXP-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B=EXP(A)
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZEXP
      DOUBLE PRECISION AR, AI, BR, BI, ZM, CA, CB
C***FIRST EXECUTABLE STATEMENT  ZEXP
      ZM = EXP(AR)
      CA = ZM*COS(AI)
      CB = ZM*SIN(AI)
      BR = CA
      BI = CB
      RETURN
      END

*DECK ZSQRT
      SUBROUTINE ZSQRT (AR, AI, BR, BI)
C***BEGIN PROLOGUE  ZSQRT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZSQRT-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     DOUBLE PRECISION COMPLEX SQUARE ROOT, B=CSQRT(A)
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  ZABS
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZSQRT
      DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DRT
      DOUBLE PRECISION ZABS
      EXTERNAL ZABS
      DATA DRT , DPI / 7.071067811865475244008443621D-1,
     1                 3.141592653589793238462643383D+0/
C***FIRST EXECUTABLE STATEMENT  ZSQRT
      ZM = ZABS(AR,AI)
      ZM = SQRT(ZM)
      IF (AR.EQ.0.0D+0) GO TO 10
      IF (AI.EQ.0.0D+0) GO TO 20
      DTHETA = DATAN(AI/AR)
      IF (DTHETA.LE.0.0D+0) GO TO 40
      IF (AR.LT.0.0D+0) DTHETA = DTHETA - DPI
      GO TO 50
   10 IF (AI.GT.0.0D+0) GO TO 60
      IF (AI.LT.0.0D+0) GO TO 70
      BR = 0.0D+0
      BI = 0.0D+0
      RETURN
   20 IF (AR.GT.0.0D+0) GO TO 30
      BR = 0.0D+0
      BI = SQRT(ABS(AR))
      RETURN
   30 BR = SQRT(AR)
      BI = 0.0D+0
      RETURN
   40 IF (AR.LT.0.0D+0) DTHETA = DTHETA + DPI
   50 DTHETA = DTHETA*0.5D+0
      BR = ZM*COS(DTHETA)
      BI = ZM*SIN(DTHETA)
      RETURN
   60 BR = ZM*DRT
      BI = ZM*DRT
      RETURN
   70 BR = ZM*DRT
      BI = -ZM*DRT
      RETURN
      END
*DECK ZUNHJ
      SUBROUTINE ZUNHJ (ZR, ZI, FNU, IPMTR, TOL, PHIR, PHII, ARGR, ARGI,
     +   ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
C***BEGIN PROLOGUE  ZUNHJ
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUNHJ-A, ZUNHJ-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     REFERENCES
C         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
C         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
C
C         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
C         PRESS, N.Y., 1974, PAGE 420
C
C     ABSTRACT
C         ZUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
C         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
C         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
C
C         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
C
C         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
C         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
C
C               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
C
C         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
C         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
C
C         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
C         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
C         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZDIV, ZLOG, ZSQRT
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZLOG and ZSQRT to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZUNHJ
C     COMPLEX ARG,ASUM,BSUM,CFNU,CONE,CR,CZERO,DR,P,PHI,PRZTH,PTFN,
C    *RFN13,RTZTA,RZTH,SUMA,SUMB,TFN,T2,UP,W,W2,Z,ZA,ZB,ZC,ZETA,ZETA1,
C    *ZETA2,ZTH
      DOUBLE PRECISION ALFA, ANG, AP, AR, ARGI, ARGR, ASUMI, ASUMR,
     * ATOL, AW2, AZTH, BETA, BR, BSUMI, BSUMR, BTOL, C, CONEI, CONER,
     * CRI, CRR, DRI, DRR, EX1, EX2, FNU, FN13, FN23, GAMA, GPI, HPI,
     * PHII, PHIR, PI, PP, PR, PRZTHI, PRZTHR, PTFNI, PTFNR, RAW, RAW2,
     * RAZTH, RFNU, RFNU2, RFN13, RTZTI, RTZTR, RZTHI, RZTHR, STI, STR,
     * SUMAI, SUMAR, SUMBI, SUMBR, TEST, TFNI, TFNR, THPI, TOL, TZAI,
     * TZAR, T2I, T2R, UPI, UPR, WI, WR, W2I, W2R, ZAI, ZAR, ZBI, ZBR,
     * ZCI, ZCR, ZEROI, ZEROR, ZETAI, ZETAR, ZETA1I, ZETA1R, ZETA2I,
     * ZETA2R, ZI, ZR, ZTHI, ZTHR, ZABS, AC, D1MACH
      INTEGER IAS, IBS, IPMTR, IS, J, JR, JU, K, KMAX, KP1, KS, L, LR,
     * LRP1, L1, L2, M, IDUM
      DIMENSION AR(14), BR(14), C(105), ALFA(180), BETA(210), GAMA(30),
     * AP(30), PR(30), PI(30), UPR(14), UPI(14), CRR(14), CRI(14),
     * DRR(14), DRI(14)
      EXTERNAL ZABS, ZLOG, ZSQRT
      DATA AR(1), AR(2), AR(3), AR(4), AR(5), AR(6), AR(7), AR(8),
     1     AR(9), AR(10), AR(11), AR(12), AR(13), AR(14)/
     2     1.00000000000000000D+00,     1.04166666666666667D-01,
     3     8.35503472222222222D-02,     1.28226574556327160D-01,
     4     2.91849026464140464D-01,     8.81627267443757652D-01,
     5     3.32140828186276754D+00,     1.49957629868625547D+01,
     6     7.89230130115865181D+01,     4.74451538868264323D+02,
     7     3.20749009089066193D+03,     2.40865496408740049D+04,
     8     1.98923119169509794D+05,     1.79190200777534383D+06/
      DATA BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),
     1     BR(9), BR(10), BR(11), BR(12), BR(13), BR(14)/
     2     1.00000000000000000D+00,    -1.45833333333333333D-01,
     3    -9.87413194444444444D-02,    -1.43312053915895062D-01,
     4    -3.17227202678413548D-01,    -9.42429147957120249D-01,
     5    -3.51120304082635426D+00,    -1.57272636203680451D+01,
     6    -8.22814390971859444D+01,    -4.92355370523670524D+02,
     7    -3.31621856854797251D+03,    -2.48276742452085896D+04,
     8    -2.04526587315129788D+05,    -1.83844491706820990D+06/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3     1.00000000000000000D+00,    -2.08333333333333333D-01,
     4     1.25000000000000000D-01,     3.34201388888888889D-01,
     5    -4.01041666666666667D-01,     7.03125000000000000D-02,
     6    -1.02581259645061728D+00,     1.84646267361111111D+00,
     7    -8.91210937500000000D-01,     7.32421875000000000D-02,
     8     4.66958442342624743D+00,    -1.12070026162229938D+01,
     9     8.78912353515625000D+00,    -2.36408691406250000D+00,
     A     1.12152099609375000D-01,    -2.82120725582002449D+01,
     B     8.46362176746007346D+01,    -9.18182415432400174D+01,
     C     4.25349987453884549D+01,    -7.36879435947963170D+00,
     D     2.27108001708984375D-01,     2.12570130039217123D+02,
     E    -7.65252468141181642D+02,     1.05999045252799988D+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3    -6.99579627376132541D+02,     2.18190511744211590D+02,
     4    -2.64914304869515555D+01,     5.72501420974731445D-01,
     5    -1.91945766231840700D+03,     8.06172218173730938D+03,
     6    -1.35865500064341374D+04,     1.16553933368645332D+04,
     7    -5.30564697861340311D+03,     1.20090291321635246D+03,
     8    -1.08090919788394656D+02,     1.72772750258445740D+00,
     9     2.02042913309661486D+04,    -9.69805983886375135D+04,
     A     1.92547001232531532D+05,    -2.03400177280415534D+05,
     B     1.22200464983017460D+05,    -4.11926549688975513D+04,
     C     7.10951430248936372D+03,    -4.93915304773088012D+02,
     D     6.07404200127348304D+00,    -2.42919187900551333D+05,
     E     1.31176361466297720D+06,    -2.99801591853810675D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3     3.76327129765640400D+06,    -2.81356322658653411D+06,
     4     1.26836527332162478D+06,    -3.31645172484563578D+05,
     5     4.52187689813627263D+04,    -2.49983048181120962D+03,
     6     2.43805296995560639D+01,     3.28446985307203782D+06,
     7    -1.97068191184322269D+07,     5.09526024926646422D+07,
     8    -7.41051482115326577D+07,     6.63445122747290267D+07,
     9    -3.75671766607633513D+07,     1.32887671664218183D+07,
     A    -2.78561812808645469D+06,     3.08186404612662398D+05,
     B    -1.38860897537170405D+04,     1.10017140269246738D+02,
     C    -4.93292536645099620D+07,     3.25573074185765749D+08,
     D    -9.39462359681578403D+08,     1.55359689957058006D+09,
     E    -1.62108055210833708D+09,     1.10684281682301447D+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3    -4.95889784275030309D+08,     1.42062907797533095D+08,
     4    -2.44740627257387285D+07,     2.24376817792244943D+06,
     5    -8.40054336030240853D+04,     5.51335896122020586D+02,
     6     8.14789096118312115D+08,    -5.86648149205184723D+09,
     7     1.86882075092958249D+10,    -3.46320433881587779D+10,
     8     4.12801855797539740D+10,    -3.30265997498007231D+10,
     9     1.79542137311556001D+10,    -6.56329379261928433D+09,
     A     1.55927986487925751D+09,    -2.25105661889415278D+08,
     B     1.73951075539781645D+07,    -5.49842327572288687D+05,
     C     3.03809051092238427D+03,    -1.46792612476956167D+10,
     D     1.14498237732025810D+11,    -3.99096175224466498D+11,
     E     8.19218669548577329D+11,    -1.09837515608122331D+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1     C(105)/
     2     1.00815810686538209D+12,    -6.45364869245376503D+11,
     3     2.87900649906150589D+11,    -8.78670721780232657D+10,
     4     1.76347306068349694D+10,    -2.16716498322379509D+09,
     5     1.43157876718888981D+08,    -3.87183344257261262D+06,
     6     1.82577554742931747D+04/
      DATA ALFA(1), ALFA(2), ALFA(3), ALFA(4), ALFA(5), ALFA(6),
     1     ALFA(7), ALFA(8), ALFA(9), ALFA(10), ALFA(11), ALFA(12),
     2     ALFA(13), ALFA(14), ALFA(15), ALFA(16), ALFA(17), ALFA(18),
     3     ALFA(19), ALFA(20), ALFA(21), ALFA(22)/
     4    -4.44444444444444444D-03,    -9.22077922077922078D-04,
     5    -8.84892884892884893D-05,     1.65927687832449737D-04,
     6     2.46691372741792910D-04,     2.65995589346254780D-04,
     7     2.61824297061500945D-04,     2.48730437344655609D-04,
     8     2.32721040083232098D-04,     2.16362485712365082D-04,
     9     2.00738858762752355D-04,     1.86267636637545172D-04,
     A     1.73060775917876493D-04,     1.61091705929015752D-04,
     B     1.50274774160908134D-04,     1.40503497391269794D-04,
     C     1.31668816545922806D-04,     1.23667445598253261D-04,
     D     1.16405271474737902D-04,     1.09798298372713369D-04,
     E     1.03772410422992823D-04,     9.82626078369363448D-05/
      DATA ALFA(23), ALFA(24), ALFA(25), ALFA(26), ALFA(27), ALFA(28),
     1     ALFA(29), ALFA(30), ALFA(31), ALFA(32), ALFA(33), ALFA(34),
     2     ALFA(35), ALFA(36), ALFA(37), ALFA(38), ALFA(39), ALFA(40),
     3     ALFA(41), ALFA(42), ALFA(43), ALFA(44)/
     4     9.32120517249503256D-05,     8.85710852478711718D-05,
     5     8.42963105715700223D-05,     8.03497548407791151D-05,
     6     7.66981345359207388D-05,     7.33122157481777809D-05,
     7     7.01662625163141333D-05,     6.72375633790160292D-05,
     8     6.93735541354588974D-04,     2.32241745182921654D-04,
     9    -1.41986273556691197D-05,    -1.16444931672048640D-04,
     A    -1.50803558053048762D-04,    -1.55121924918096223D-04,
     B    -1.46809756646465549D-04,    -1.33815503867491367D-04,
     C    -1.19744975684254051D-04,    -1.06184319207974020D-04,
     D    -9.37699549891194492D-05,    -8.26923045588193274D-05,
     E    -7.29374348155221211D-05,    -6.44042357721016283D-05/
      DATA ALFA(45), ALFA(46), ALFA(47), ALFA(48), ALFA(49), ALFA(50),
     1     ALFA(51), ALFA(52), ALFA(53), ALFA(54), ALFA(55), ALFA(56),
     2     ALFA(57), ALFA(58), ALFA(59), ALFA(60), ALFA(61), ALFA(62),
     3     ALFA(63), ALFA(64), ALFA(65), ALFA(66)/
     4    -5.69611566009369048D-05,    -5.04731044303561628D-05,
     5    -4.48134868008882786D-05,    -3.98688727717598864D-05,
     6    -3.55400532972042498D-05,    -3.17414256609022480D-05,
     7    -2.83996793904174811D-05,    -2.54522720634870566D-05,
     8    -2.28459297164724555D-05,    -2.05352753106480604D-05,
     9    -1.84816217627666085D-05,    -1.66519330021393806D-05,
     A    -1.50179412980119482D-05,    -1.35554031379040526D-05,
     B    -1.22434746473858131D-05,    -1.10641884811308169D-05,
     C    -3.54211971457743841D-04,    -1.56161263945159416D-04,
     D     3.04465503594936410D-05,     1.30198655773242693D-04,
     E     1.67471106699712269D-04,     1.70222587683592569D-04/
      DATA ALFA(67), ALFA(68), ALFA(69), ALFA(70), ALFA(71), ALFA(72),
     1     ALFA(73), ALFA(74), ALFA(75), ALFA(76), ALFA(77), ALFA(78),
     2     ALFA(79), ALFA(80), ALFA(81), ALFA(82), ALFA(83), ALFA(84),
     3     ALFA(85), ALFA(86), ALFA(87), ALFA(88)/
     4     1.56501427608594704D-04,     1.36339170977445120D-04,
     5     1.14886692029825128D-04,     9.45869093034688111D-05,
     6     7.64498419250898258D-05,     6.07570334965197354D-05,
     7     4.74394299290508799D-05,     3.62757512005344297D-05,
     8     2.69939714979224901D-05,     1.93210938247939253D-05,
     9     1.30056674793963203D-05,     7.82620866744496661D-06,
     A     3.59257485819351583D-06,     1.44040049814251817D-07,
     B    -2.65396769697939116D-06,    -4.91346867098485910D-06,
     C    -6.72739296091248287D-06,    -8.17269379678657923D-06,
     D    -9.31304715093561232D-06,    -1.02011418798016441D-05,
     E    -1.08805962510592880D-05,    -1.13875481509603555D-05/
      DATA ALFA(89), ALFA(90), ALFA(91), ALFA(92), ALFA(93), ALFA(94),
     1     ALFA(95), ALFA(96), ALFA(97), ALFA(98), ALFA(99), ALFA(100),
     2     ALFA(101), ALFA(102), ALFA(103), ALFA(104), ALFA(105),
     3     ALFA(106), ALFA(107), ALFA(108), ALFA(109), ALFA(110)/
     4    -1.17519675674556414D-05,    -1.19987364870944141D-05,
     5     3.78194199201772914D-04,     2.02471952761816167D-04,
     6    -6.37938506318862408D-05,    -2.38598230603005903D-04,
     7    -3.10916256027361568D-04,    -3.13680115247576316D-04,
     8    -2.78950273791323387D-04,    -2.28564082619141374D-04,
     9    -1.75245280340846749D-04,    -1.25544063060690348D-04,
     A    -8.22982872820208365D-05,    -4.62860730588116458D-05,
     B    -1.72334302366962267D-05,     5.60690482304602267D-06,
     C     2.31395443148286800D-05,     3.62642745856793957D-05,
     D     4.58006124490188752D-05,     5.24595294959114050D-05,
     E     5.68396208545815266D-05,     5.94349820393104052D-05/
      DATA ALFA(111), ALFA(112), ALFA(113), ALFA(114), ALFA(115),
     1     ALFA(116), ALFA(117), ALFA(118), ALFA(119), ALFA(120),
     2     ALFA(121), ALFA(122), ALFA(123), ALFA(124), ALFA(125),
     3     ALFA(126), ALFA(127), ALFA(128), ALFA(129), ALFA(130)/
     4     6.06478527578421742D-05,     6.08023907788436497D-05,
     5     6.01577894539460388D-05,     5.89199657344698500D-05,
     6     5.72515823777593053D-05,     5.52804375585852577D-05,
     7     5.31063773802880170D-05,     5.08069302012325706D-05,
     8     4.84418647620094842D-05,     4.60568581607475370D-05,
     9    -6.91141397288294174D-04,    -4.29976633058871912D-04,
     A     1.83067735980039018D-04,     6.60088147542014144D-04,
     B     8.75964969951185931D-04,     8.77335235958235514D-04,
     C     7.49369585378990637D-04,     5.63832329756980918D-04,
     D     3.68059319971443156D-04,     1.88464535514455599D-04/
      DATA ALFA(131), ALFA(132), ALFA(133), ALFA(134), ALFA(135),
     1     ALFA(136), ALFA(137), ALFA(138), ALFA(139), ALFA(140),
     2     ALFA(141), ALFA(142), ALFA(143), ALFA(144), ALFA(145),
     3     ALFA(146), ALFA(147), ALFA(148), ALFA(149), ALFA(150)/
     4     3.70663057664904149D-05,    -8.28520220232137023D-05,
     5    -1.72751952869172998D-04,    -2.36314873605872983D-04,
     6    -2.77966150694906658D-04,    -3.02079514155456919D-04,
     7    -3.12594712643820127D-04,    -3.12872558758067163D-04,
     8    -3.05678038466324377D-04,    -2.93226470614557331D-04,
     9    -2.77255655582934777D-04,    -2.59103928467031709D-04,
     A    -2.39784014396480342D-04,    -2.20048260045422848D-04,
     B    -2.00443911094971498D-04,    -1.81358692210970687D-04,
     C    -1.63057674478657464D-04,    -1.45712672175205844D-04,
     D    -1.29425421983924587D-04,    -1.14245691942445952D-04/
      DATA ALFA(151), ALFA(152), ALFA(153), ALFA(154), ALFA(155),
     1     ALFA(156), ALFA(157), ALFA(158), ALFA(159), ALFA(160),
     2     ALFA(161), ALFA(162), ALFA(163), ALFA(164), ALFA(165),
     3     ALFA(166), ALFA(167), ALFA(168), ALFA(169), ALFA(170)/
     4     1.92821964248775885D-03,     1.35592576302022234D-03,
     5    -7.17858090421302995D-04,    -2.58084802575270346D-03,
     6    -3.49271130826168475D-03,    -3.46986299340960628D-03,
     7    -2.82285233351310182D-03,    -1.88103076404891354D-03,
     8    -8.89531718383947600D-04,     3.87912102631035228D-06,
     9     7.28688540119691412D-04,     1.26566373053457758D-03,
     A     1.62518158372674427D-03,     1.83203153216373172D-03,
     B     1.91588388990527909D-03,     1.90588846755546138D-03,
     C     1.82798982421825727D-03,     1.70389506421121530D-03,
     D     1.55097127171097686D-03,     1.38261421852276159D-03/
      DATA ALFA(171), ALFA(172), ALFA(173), ALFA(174), ALFA(175),
     1     ALFA(176), ALFA(177), ALFA(178), ALFA(179), ALFA(180)/
     2     1.20881424230064774D-03,     1.03676532638344962D-03,
     3     8.71437918068619115D-04,     7.16080155297701002D-04,
     4     5.72637002558129372D-04,     4.42089819465802277D-04,
     5     3.24724948503090564D-04,     2.20342042730246599D-04,
     6     1.28412898401353882D-04,     4.82005924552095464D-05/
      DATA BETA(1), BETA(2), BETA(3), BETA(4), BETA(5), BETA(6),
     1     BETA(7), BETA(8), BETA(9), BETA(10), BETA(11), BETA(12),
     2     BETA(13), BETA(14), BETA(15), BETA(16), BETA(17), BETA(18),
     3     BETA(19), BETA(20), BETA(21), BETA(22)/
     4     1.79988721413553309D-02,     5.59964911064388073D-03,
     5     2.88501402231132779D-03,     1.80096606761053941D-03,
     6     1.24753110589199202D-03,     9.22878876572938311D-04,
     7     7.14430421727287357D-04,     5.71787281789704872D-04,
     8     4.69431007606481533D-04,     3.93232835462916638D-04,
     9     3.34818889318297664D-04,     2.88952148495751517D-04,
     A     2.52211615549573284D-04,     2.22280580798883327D-04,
     B     1.97541838033062524D-04,     1.76836855019718004D-04,
     C     1.59316899661821081D-04,     1.44347930197333986D-04,
     D     1.31448068119965379D-04,     1.20245444949302884D-04,
     E     1.10449144504599392D-04,     1.01828770740567258D-04/
      DATA BETA(23), BETA(24), BETA(25), BETA(26), BETA(27), BETA(28),
     1     BETA(29), BETA(30), BETA(31), BETA(32), BETA(33), BETA(34),
     2     BETA(35), BETA(36), BETA(37), BETA(38), BETA(39), BETA(40),
     3     BETA(41), BETA(42), BETA(43), BETA(44)/
     4     9.41998224204237509D-05,     8.74130545753834437D-05,
     5     8.13466262162801467D-05,     7.59002269646219339D-05,
     6     7.09906300634153481D-05,     6.65482874842468183D-05,
     7     6.25146958969275078D-05,     5.88403394426251749D-05,
     8    -1.49282953213429172D-03,    -8.78204709546389328D-04,
     9    -5.02916549572034614D-04,    -2.94822138512746025D-04,
     A    -1.75463996970782828D-04,    -1.04008550460816434D-04,
     B    -5.96141953046457895D-05,    -3.12038929076098340D-05,
     C    -1.26089735980230047D-05,    -2.42892608575730389D-07,
     D     8.05996165414273571D-06,     1.36507009262147391D-05,
     E     1.73964125472926261D-05,     1.98672978842133780D-05/
      DATA BETA(45), BETA(46), BETA(47), BETA(48), BETA(49), BETA(50),
     1     BETA(51), BETA(52), BETA(53), BETA(54), BETA(55), BETA(56),
     2     BETA(57), BETA(58), BETA(59), BETA(60), BETA(61), BETA(62),
     3     BETA(63), BETA(64), BETA(65), BETA(66)/
     4     2.14463263790822639D-05,     2.23954659232456514D-05,
     5     2.28967783814712629D-05,     2.30785389811177817D-05,
     6     2.30321976080909144D-05,     2.28236073720348722D-05,
     7     2.25005881105292418D-05,     2.20981015361991429D-05,
     8     2.16418427448103905D-05,     2.11507649256220843D-05,
     9     2.06388749782170737D-05,     2.01165241997081666D-05,
     A     1.95913450141179244D-05,     1.90689367910436740D-05,
     B     1.85533719641636667D-05,     1.80475722259674218D-05,
     C     5.52213076721292790D-04,     4.47932581552384646D-04,
     D     2.79520653992020589D-04,     1.52468156198446602D-04,
     E     6.93271105657043598D-05,     1.76258683069991397D-05/
      DATA BETA(67), BETA(68), BETA(69), BETA(70), BETA(71), BETA(72),
     1     BETA(73), BETA(74), BETA(75), BETA(76), BETA(77), BETA(78),
     2     BETA(79), BETA(80), BETA(81), BETA(82), BETA(83), BETA(84),
     3     BETA(85), BETA(86), BETA(87), BETA(88)/
     4    -1.35744996343269136D-05,    -3.17972413350427135D-05,
     5    -4.18861861696693365D-05,    -4.69004889379141029D-05,
     6    -4.87665447413787352D-05,    -4.87010031186735069D-05,
     7    -4.74755620890086638D-05,    -4.55813058138628452D-05,
     8    -4.33309644511266036D-05,    -4.09230193157750364D-05,
     9    -3.84822638603221274D-05,    -3.60857167535410501D-05,
     A    -3.37793306123367417D-05,    -3.15888560772109621D-05,
     B    -2.95269561750807315D-05,    -2.75978914828335759D-05,
     C    -2.58006174666883713D-05,    -2.41308356761280200D-05,
     D    -2.25823509518346033D-05,    -2.11479656768912971D-05,
     E    -1.98200638885294927D-05,    -1.85909870801065077D-05/
      DATA BETA(89), BETA(90), BETA(91), BETA(92), BETA(93), BETA(94),
     1     BETA(95), BETA(96), BETA(97), BETA(98), BETA(99), BETA(100),
     2     BETA(101), BETA(102), BETA(103), BETA(104), BETA(105),
     3     BETA(106), BETA(107), BETA(108), BETA(109), BETA(110)/
     4    -1.74532699844210224D-05,    -1.63997823854497997D-05,
     5    -4.74617796559959808D-04,    -4.77864567147321487D-04,
     6    -3.20390228067037603D-04,    -1.61105016119962282D-04,
     7    -4.25778101285435204D-05,     3.44571294294967503D-05,
     8     7.97092684075674924D-05,     1.03138236708272200D-04,
     9     1.12466775262204158D-04,     1.13103642108481389D-04,
     A     1.08651634848774268D-04,     1.01437951597661973D-04,
     B     9.29298396593363896D-05,     8.40293133016089978D-05,
     C     7.52727991349134062D-05,     6.69632521975730872D-05,
     D     5.92564547323194704D-05,     5.22169308826975567D-05,
     E     4.58539485165360646D-05,     4.01445513891486808D-05/
      DATA BETA(111), BETA(112), BETA(113), BETA(114), BETA(115),
     1     BETA(116), BETA(117), BETA(118), BETA(119), BETA(120),
     2     BETA(121), BETA(122), BETA(123), BETA(124), BETA(125),
     3     BETA(126), BETA(127), BETA(128), BETA(129), BETA(130)/
     4     3.50481730031328081D-05,     3.05157995034346659D-05,
     5     2.64956119950516039D-05,     2.29363633690998152D-05,
     6     1.97893056664021636D-05,     1.70091984636412623D-05,
     7     1.45547428261524004D-05,     1.23886640995878413D-05,
     8     1.04775876076583236D-05,     8.79179954978479373D-06,
     9     7.36465810572578444D-04,     8.72790805146193976D-04,
     A     6.22614862573135066D-04,     2.85998154194304147D-04,
     B     3.84737672879366102D-06,    -1.87906003636971558D-04,
     C    -2.97603646594554535D-04,    -3.45998126832656348D-04,
     D    -3.53382470916037712D-04,    -3.35715635775048757D-04/
      DATA BETA(131), BETA(132), BETA(133), BETA(134), BETA(135),
     1     BETA(136), BETA(137), BETA(138), BETA(139), BETA(140),
     2     BETA(141), BETA(142), BETA(143), BETA(144), BETA(145),
     3     BETA(146), BETA(147), BETA(148), BETA(149), BETA(150)/
     4    -3.04321124789039809D-04,    -2.66722723047612821D-04,
     5    -2.27654214122819527D-04,    -1.89922611854562356D-04,
     6    -1.55058918599093870D-04,    -1.23778240761873630D-04,
     7    -9.62926147717644187D-05,    -7.25178327714425337D-05,
     8    -5.22070028895633801D-05,    -3.50347750511900522D-05,
     9    -2.06489761035551757D-05,    -8.70106096849767054D-06,
     A     1.13698686675100290D-06,     9.16426474122778849D-06,
     B     1.56477785428872620D-05,     2.08223629482466847D-05,
     C     2.48923381004595156D-05,     2.80340509574146325D-05,
     D     3.03987774629861915D-05,     3.21156731406700616D-05/
      DATA BETA(151), BETA(152), BETA(153), BETA(154), BETA(155),
     1     BETA(156), BETA(157), BETA(158), BETA(159), BETA(160),
     2     BETA(161), BETA(162), BETA(163), BETA(164), BETA(165),
     3     BETA(166), BETA(167), BETA(168), BETA(169), BETA(170)/
     4    -1.80182191963885708D-03,    -2.43402962938042533D-03,
     5    -1.83422663549856802D-03,    -7.62204596354009765D-04,
     6     2.39079475256927218D-04,     9.49266117176881141D-04,
     7     1.34467449701540359D-03,     1.48457495259449178D-03,
     8     1.44732339830617591D-03,     1.30268261285657186D-03,
     9     1.10351597375642682D-03,     8.86047440419791759D-04,
     A     6.73073208165665473D-04,     4.77603872856582378D-04,
     B     3.05991926358789362D-04,     1.60315694594721630D-04,
     C     4.00749555270613286D-05,    -5.66607461635251611D-05,
     D    -1.32506186772982638D-04,    -1.90296187989614057D-04/
      DATA BETA(171), BETA(172), BETA(173), BETA(174), BETA(175),
     1     BETA(176), BETA(177), BETA(178), BETA(179), BETA(180),
     2     BETA(181), BETA(182), BETA(183), BETA(184), BETA(185),
     3     BETA(186), BETA(187), BETA(188), BETA(189), BETA(190)/
     4    -2.32811450376937408D-04,    -2.62628811464668841D-04,
     5    -2.82050469867598672D-04,    -2.93081563192861167D-04,
     6    -2.97435962176316616D-04,    -2.96557334239348078D-04,
     7    -2.91647363312090861D-04,    -2.83696203837734166D-04,
     8    -2.73512317095673346D-04,    -2.61750155806768580D-04,
     9     6.38585891212050914D-03,     9.62374215806377941D-03,
     A     7.61878061207001043D-03,     2.83219055545628054D-03,
     B    -2.09841352012720090D-03,    -5.73826764216626498D-03,
     C    -7.70804244495414620D-03,    -8.21011692264844401D-03,
     D    -7.65824520346905413D-03,    -6.47209729391045177D-03/
      DATA BETA(191), BETA(192), BETA(193), BETA(194), BETA(195),
     1     BETA(196), BETA(197), BETA(198), BETA(199), BETA(200),
     2     BETA(201), BETA(202), BETA(203), BETA(204), BETA(205),
     3     BETA(206), BETA(207), BETA(208), BETA(209), BETA(210)/
     4    -4.99132412004966473D-03,    -3.45612289713133280D-03,
     5    -2.01785580014170775D-03,    -7.59430686781961401D-04,
     6     2.84173631523859138D-04,     1.10891667586337403D-03,
     7     1.72901493872728771D-03,     2.16812590802684701D-03,
     8     2.45357710494539735D-03,     2.61281821058334862D-03,
     9     2.67141039656276912D-03,     2.65203073395980430D-03,
     A     2.57411652877287315D-03,     2.45389126236094427D-03,
     B     2.30460058071795494D-03,     2.13684837686712662D-03,
     C     1.95896528478870911D-03,     1.77737008679454412D-03,
     D     1.59690280765839059D-03,     1.42111975664438546D-03/
      DATA GAMA(1), GAMA(2), GAMA(3), GAMA(4), GAMA(5), GAMA(6),
     1     GAMA(7), GAMA(8), GAMA(9), GAMA(10), GAMA(11), GAMA(12),
     2     GAMA(13), GAMA(14), GAMA(15), GAMA(16), GAMA(17), GAMA(18),
     3     GAMA(19), GAMA(20), GAMA(21), GAMA(22)/
     4     6.29960524947436582D-01,     2.51984209978974633D-01,
     5     1.54790300415655846D-01,     1.10713062416159013D-01,
     6     8.57309395527394825D-02,     6.97161316958684292D-02,
     7     5.86085671893713576D-02,     5.04698873536310685D-02,
     8     4.42600580689154809D-02,     3.93720661543509966D-02,
     9     3.54283195924455368D-02,     3.21818857502098231D-02,
     A     2.94646240791157679D-02,     2.71581677112934479D-02,
     B     2.51768272973861779D-02,     2.34570755306078891D-02,
     C     2.19508390134907203D-02,     2.06210828235646240D-02,
     D     1.94388240897880846D-02,     1.83810633800683158D-02,
     E     1.74293213231963172D-02,     1.65685837786612353D-02/
      DATA GAMA(23), GAMA(24), GAMA(25), GAMA(26), GAMA(27), GAMA(28),
     1     GAMA(29), GAMA(30)/
     2     1.57865285987918445D-02,     1.50729501494095594D-02,
     3     1.44193250839954639D-02,     1.38184805735341786D-02,
     4     1.32643378994276568D-02,     1.27517121970498651D-02,
     5     1.22761545318762767D-02,     1.18338262398482403D-02/
      DATA EX1, EX2, HPI, GPI, THPI /
     1     3.33333333333333333D-01,     6.66666666666666667D-01,
     2     1.57079632679489662D+00,     3.14159265358979324D+00,
     3     4.71238898038468986D+00/
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZUNHJ
      RFNU = 1.0D0/FNU
C-----------------------------------------------------------------------
C     OVERFLOW TEST (Z/FNU TOO SMALL)
C-----------------------------------------------------------------------
      TEST = D1MACH(1)*1.0D+3
      AC = FNU*TEST
      IF (ABS(ZR).GT.AC .OR. ABS(ZI).GT.AC) GO TO 15
      ZETA1R = 2.0D0*ABS(LOG(TEST))+FNU
      ZETA1I = 0.0D0
      ZETA2R = FNU
      ZETA2I = 0.0D0
      PHIR = 1.0D0
      PHII = 0.0D0
      ARGR = 1.0D0
      ARGI = 0.0D0
      RETURN
   15 CONTINUE
      ZBR = ZR*RFNU
      ZBI = ZI*RFNU
      RFNU2 = RFNU*RFNU
C-----------------------------------------------------------------------
C     COMPUTE IN THE FOURTH QUADRANT
C-----------------------------------------------------------------------
      FN13 = FNU**EX1
      FN23 = FN13*FN13
      RFN13 = 1.0D0/FN13
      W2R = CONER - ZBR*ZBR + ZBI*ZBI
      W2I = CONEI - ZBR*ZBI - ZBR*ZBI
      AW2 = ZABS(W2R,W2I)
      IF (AW2.GT.0.25D0) GO TO 130
C-----------------------------------------------------------------------
C     POWER SERIES FOR ABS(W2).LE.0.25D0
C-----------------------------------------------------------------------
      K = 1
      PR(1) = CONER
      PI(1) = CONEI
      SUMAR = GAMA(1)
      SUMAI = ZEROI
      AP(1) = 1.0D0
      IF (AW2.LT.TOL) GO TO 20
      DO 10 K=2,30
        PR(K) = PR(K-1)*W2R - PI(K-1)*W2I
        PI(K) = PR(K-1)*W2I + PI(K-1)*W2R
        SUMAR = SUMAR + PR(K)*GAMA(K)
        SUMAI = SUMAI + PI(K)*GAMA(K)
        AP(K) = AP(K-1)*AW2
        IF (AP(K).LT.TOL) GO TO 20
   10 CONTINUE
      K = 30
   20 CONTINUE
      KMAX = K
      ZETAR = W2R*SUMAR - W2I*SUMAI
      ZETAI = W2R*SUMAI + W2I*SUMAR
      ARGR = ZETAR*FN23
      ARGI = ZETAI*FN23
      CALL ZSQRT(SUMAR, SUMAI, ZAR, ZAI)
      CALL ZSQRT(W2R, W2I, STR, STI)
      ZETA2R = STR*FNU
      ZETA2I = STI*FNU
      STR = CONER + EX2*(ZETAR*ZAR-ZETAI*ZAI)
      STI = CONEI + EX2*(ZETAR*ZAI+ZETAI*ZAR)
      ZETA1R = STR*ZETA2R - STI*ZETA2I
      ZETA1I = STR*ZETA2I + STI*ZETA2R
      ZAR = ZAR + ZAR
      ZAI = ZAI + ZAI
      CALL ZSQRT(ZAR, ZAI, STR, STI)
      PHIR = STR*RFN13
      PHII = STI*RFN13
      IF (IPMTR.EQ.1) GO TO 120
C-----------------------------------------------------------------------
C     SUM SERIES FOR ASUM AND BSUM
C-----------------------------------------------------------------------
      SUMBR = ZEROR
      SUMBI = ZEROI
      DO 30 K=1,KMAX
        SUMBR = SUMBR + PR(K)*BETA(K)
        SUMBI = SUMBI + PI(K)*BETA(K)
   30 CONTINUE
      ASUMR = ZEROR
      ASUMI = ZEROI
      BSUMR = SUMBR
      BSUMI = SUMBI
      L1 = 0
      L2 = 30
      BTOL = TOL*(ABS(BSUMR)+ABS(BSUMI))
      ATOL = TOL
      PP = 1.0D0
      IAS = 0
      IBS = 0
      IF (RFNU2.LT.TOL) GO TO 110
      DO 100 IS=2,7
        ATOL = ATOL/RFNU2
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 60
        SUMAR = ZEROR
        SUMAI = ZEROI
        DO 40 K=1,KMAX
          M = L1 + K
          SUMAR = SUMAR + PR(K)*ALFA(M)
          SUMAI = SUMAI + PI(K)*ALFA(M)
          IF (AP(K).LT.ATOL) GO TO 50
   40   CONTINUE
   50   CONTINUE
        ASUMR = ASUMR + SUMAR*PP
        ASUMI = ASUMI + SUMAI*PP
        IF (PP.LT.TOL) IAS = 1
   60   CONTINUE
        IF (IBS.EQ.1) GO TO 90
        SUMBR = ZEROR
        SUMBI = ZEROI
        DO 70 K=1,KMAX
          M = L2 + K
          SUMBR = SUMBR + PR(K)*BETA(M)
          SUMBI = SUMBI + PI(K)*BETA(M)
          IF (AP(K).LT.ATOL) GO TO 80
   70   CONTINUE
   80   CONTINUE
        BSUMR = BSUMR + SUMBR*PP
        BSUMI = BSUMI + SUMBI*PP
        IF (PP.LT.BTOL) IBS = 1
   90   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 110
        L1 = L1 + 30
        L2 = L2 + 30
  100 CONTINUE
  110 CONTINUE
      ASUMR = ASUMR + CONER
      PP = RFNU*RFN13
      BSUMR = BSUMR*PP
      BSUMI = BSUMI*PP
  120 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     ABS(W2).GT.0.25D0
C-----------------------------------------------------------------------
  130 CONTINUE
      CALL ZSQRT(W2R, W2I, WR, WI)
      IF (WR.LT.0.0D0) WR = 0.0D0
      IF (WI.LT.0.0D0) WI = 0.0D0
      STR = CONER + WR
      STI = WI
      CALL ZDIV(STR, STI, ZBR, ZBI, ZAR, ZAI)
      CALL ZLOG(ZAR, ZAI, ZCR, ZCI, IDUM)
      IF (ZCI.LT.0.0D0) ZCI = 0.0D0
      IF (ZCI.GT.HPI) ZCI = HPI
      IF (ZCR.LT.0.0D0) ZCR = 0.0D0
      ZTHR = (ZCR-WR)*1.5D0
      ZTHI = (ZCI-WI)*1.5D0
      ZETA1R = ZCR*FNU
      ZETA1I = ZCI*FNU
      ZETA2R = WR*FNU
      ZETA2I = WI*FNU
      AZTH = ZABS(ZTHR,ZTHI)
      ANG = THPI
      IF (ZTHR.GE.0.0D0 .AND. ZTHI.LT.0.0D0) GO TO 140
      ANG = HPI
      IF (ZTHR.EQ.0.0D0) GO TO 140
      ANG = DATAN(ZTHI/ZTHR)
      IF (ZTHR.LT.0.0D0) ANG = ANG + GPI
  140 CONTINUE
      PP = AZTH**EX2
      ANG = ANG*EX2
      ZETAR = PP*COS(ANG)
      ZETAI = PP*SIN(ANG)
      IF (ZETAI.LT.0.0D0) ZETAI = 0.0D0
      ARGR = ZETAR*FN23
      ARGI = ZETAI*FN23
      CALL ZDIV(ZTHR, ZTHI, ZETAR, ZETAI, RTZTR, RTZTI)
      CALL ZDIV(RTZTR, RTZTI, WR, WI, ZAR, ZAI)
      TZAR = ZAR + ZAR
      TZAI = ZAI + ZAI
      CALL ZSQRT(TZAR, TZAI, STR, STI)
      PHIR = STR*RFN13
      PHII = STI*RFN13
      IF (IPMTR.EQ.1) GO TO 120
      RAW = 1.0D0/SQRT(AW2)
      STR = WR*RAW
      STI = -WI*RAW
      TFNR = STR*RFNU*RAW
      TFNI = STI*RFNU*RAW
      RAZTH = 1.0D0/AZTH
      STR = ZTHR*RAZTH
      STI = -ZTHI*RAZTH
      RZTHR = STR*RAZTH*RFNU
      RZTHI = STI*RAZTH*RFNU
      ZCR = RZTHR*AR(2)
      ZCI = RZTHI*AR(2)
      RAW2 = 1.0D0/AW2
      STR = W2R*RAW2
      STI = -W2I*RAW2
      T2R = STR*RAW2
      T2I = STI*RAW2
      STR = T2R*C(2) + C(3)
      STI = T2I*C(2)
      UPR(2) = STR*TFNR - STI*TFNI
      UPI(2) = STR*TFNI + STI*TFNR
      BSUMR = UPR(2) + ZCR
      BSUMI = UPI(2) + ZCI
      ASUMR = ZEROR
      ASUMI = ZEROI
      IF (RFNU.LT.TOL) GO TO 220
      PRZTHR = RZTHR
      PRZTHI = RZTHI
      PTFNR = TFNR
      PTFNI = TFNI
      UPR(1) = CONER
      UPI(1) = CONEI
      PP = 1.0D0
      BTOL = TOL*(ABS(BSUMR)+ABS(BSUMI))
      KS = 0
      KP1 = 2
      L = 3
      IAS = 0
      IBS = 0
      DO 210 LR=2,12,2
        LRP1 = LR + 1
C-----------------------------------------------------------------------
C     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
C     NEXT SUMA AND SUMB
C-----------------------------------------------------------------------
        DO 160 K=LR,LRP1
          KS = KS + 1
          KP1 = KP1 + 1
          L = L + 1
          ZAR = C(L)
          ZAI = ZEROI
          DO 150 J=2,KP1
            L = L + 1
            STR = ZAR*T2R - T2I*ZAI + C(L)
            ZAI = ZAR*T2I + ZAI*T2R
            ZAR = STR
  150     CONTINUE
          STR = PTFNR*TFNR - PTFNI*TFNI
          PTFNI = PTFNR*TFNI + PTFNI*TFNR
          PTFNR = STR
          UPR(KP1) = PTFNR*ZAR - PTFNI*ZAI
          UPI(KP1) = PTFNI*ZAR + PTFNR*ZAI
          CRR(KS) = PRZTHR*BR(KS+1)
          CRI(KS) = PRZTHI*BR(KS+1)
          STR = PRZTHR*RZTHR - PRZTHI*RZTHI
          PRZTHI = PRZTHR*RZTHI + PRZTHI*RZTHR
          PRZTHR = STR
          DRR(KS) = PRZTHR*AR(KS+2)
          DRI(KS) = PRZTHI*AR(KS+2)
  160   CONTINUE
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 180
        SUMAR = UPR(LRP1)
        SUMAI = UPI(LRP1)
        JU = LRP1
        DO 170 JR=1,LR
          JU = JU - 1
          SUMAR = SUMAR + CRR(JR)*UPR(JU) - CRI(JR)*UPI(JU)
          SUMAI = SUMAI + CRR(JR)*UPI(JU) + CRI(JR)*UPR(JU)
  170   CONTINUE
        ASUMR = ASUMR + SUMAR
        ASUMI = ASUMI + SUMAI
        TEST = ABS(SUMAR) + ABS(SUMAI)
        IF (PP.LT.TOL .AND. TEST.LT.TOL) IAS = 1
  180   CONTINUE
        IF (IBS.EQ.1) GO TO 200
        SUMBR = UPR(LR+2) + UPR(LRP1)*ZCR - UPI(LRP1)*ZCI
        SUMBI = UPI(LR+2) + UPR(LRP1)*ZCI + UPI(LRP1)*ZCR
        JU = LRP1
        DO 190 JR=1,LR
          JU = JU - 1
          SUMBR = SUMBR + DRR(JR)*UPR(JU) - DRI(JR)*UPI(JU)
          SUMBI = SUMBI + DRR(JR)*UPI(JU) + DRI(JR)*UPR(JU)
  190   CONTINUE
        BSUMR = BSUMR + SUMBR
        BSUMI = BSUMI + SUMBI
        TEST = ABS(SUMBR) + ABS(SUMBI)
        IF (PP.LT.BTOL .AND. TEST.LT.BTOL) IBS = 1
  200   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 220
  210 CONTINUE
  220 CONTINUE
      ASUMR = ASUMR + CONER
      STR = -BSUMR*RFN13
      STI = -BSUMI*RFN13
      CALL ZDIV(STR, STI, RTZTR, RTZTI, BSUMR, BSUMI)
      GO TO 120
      END



*DECK ZUNIK
      SUBROUTINE ZUNIK (ZRR, ZRI, FNU, IKFLG, IPMTR, TOL, INIT, PHIR,
     +   PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
C***BEGIN PROLOGUE  ZUNIK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUNIK-A, ZUNIK-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C        ZUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
C        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
C        RESPECTIVELY BY
C
C        W(FNU,ZR) = PHI*EXP(ZETA)*SUM
C
C        WHERE       ZETA=-ZETA1 + ZETA2       OR
C                          ZETA1 - ZETA2
C
C        THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
C        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG=
C        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
C        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
C        ZETA1,ZETA2.
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZDIV, ZLOG, ZSQRT
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added EXTERNAL statement with ZLOG and ZSQRT.  (RWC)
C***END PROLOGUE  ZUNIK
C     COMPLEX CFN,CON,CONE,CRFN,CWRK,CZERO,PHI,S,SR,SUM,T,T2,ZETA1,
C    *ZETA2,ZN,ZR
      DOUBLE PRECISION AC, C, CON, CONEI, CONER, CRFNI, CRFNR, CWRKI,
     * CWRKR, FNU, PHII, PHIR, RFN, SI, SR, SRI, SRR, STI, STR, SUMI,
     * SUMR, TEST, TI, TOL, TR, T2I, T2R, ZEROI, ZEROR, ZETA1I, ZETA1R,
     * ZETA2I, ZETA2R, ZNI, ZNR, ZRI, ZRR, D1MACH
      INTEGER I, IDUM, IKFLG, INIT, IPMTR, J, K, L
      DIMENSION C(120), CWRKR(16), CWRKI(16), CON(2)
      EXTERNAL ZLOG, ZSQRT
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
      DATA CON(1), CON(2)  /
     1 3.98942280401432678D-01,  1.25331413731550025D+00 /
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3     1.00000000000000000D+00,    -2.08333333333333333D-01,
     4     1.25000000000000000D-01,     3.34201388888888889D-01,
     5    -4.01041666666666667D-01,     7.03125000000000000D-02,
     6    -1.02581259645061728D+00,     1.84646267361111111D+00,
     7    -8.91210937500000000D-01,     7.32421875000000000D-02,
     8     4.66958442342624743D+00,    -1.12070026162229938D+01,
     9     8.78912353515625000D+00,    -2.36408691406250000D+00,
     A     1.12152099609375000D-01,    -2.82120725582002449D+01,
     B     8.46362176746007346D+01,    -9.18182415432400174D+01,
     C     4.25349987453884549D+01,    -7.36879435947963170D+00,
     D     2.27108001708984375D-01,     2.12570130039217123D+02,
     E    -7.65252468141181642D+02,     1.05999045252799988D+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3    -6.99579627376132541D+02,     2.18190511744211590D+02,
     4    -2.64914304869515555D+01,     5.72501420974731445D-01,
     5    -1.91945766231840700D+03,     8.06172218173730938D+03,
     6    -1.35865500064341374D+04,     1.16553933368645332D+04,
     7    -5.30564697861340311D+03,     1.20090291321635246D+03,
     8    -1.08090919788394656D+02,     1.72772750258445740D+00,
     9     2.02042913309661486D+04,    -9.69805983886375135D+04,
     A     1.92547001232531532D+05,    -2.03400177280415534D+05,
     B     1.22200464983017460D+05,    -4.11926549688975513D+04,
     C     7.10951430248936372D+03,    -4.93915304773088012D+02,
     D     6.07404200127348304D+00,    -2.42919187900551333D+05,
     E     1.31176361466297720D+06,    -2.99801591853810675D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3     3.76327129765640400D+06,    -2.81356322658653411D+06,
     4     1.26836527332162478D+06,    -3.31645172484563578D+05,
     5     4.52187689813627263D+04,    -2.49983048181120962D+03,
     6     2.43805296995560639D+01,     3.28446985307203782D+06,
     7    -1.97068191184322269D+07,     5.09526024926646422D+07,
     8    -7.41051482115326577D+07,     6.63445122747290267D+07,
     9    -3.75671766607633513D+07,     1.32887671664218183D+07,
     A    -2.78561812808645469D+06,     3.08186404612662398D+05,
     B    -1.38860897537170405D+04,     1.10017140269246738D+02,
     C    -4.93292536645099620D+07,     3.25573074185765749D+08,
     D    -9.39462359681578403D+08,     1.55359689957058006D+09,
     E    -1.62108055210833708D+09,     1.10684281682301447D+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3    -4.95889784275030309D+08,     1.42062907797533095D+08,
     4    -2.44740627257387285D+07,     2.24376817792244943D+06,
     5    -8.40054336030240853D+04,     5.51335896122020586D+02,
     6     8.14789096118312115D+08,    -5.86648149205184723D+09,
     7     1.86882075092958249D+10,    -3.46320433881587779D+10,
     8     4.12801855797539740D+10,    -3.30265997498007231D+10,
     9     1.79542137311556001D+10,    -6.56329379261928433D+09,
     A     1.55927986487925751D+09,    -2.25105661889415278D+08,
     B     1.73951075539781645D+07,    -5.49842327572288687D+05,
     C     3.03809051092238427D+03,    -1.46792612476956167D+10,
     D     1.14498237732025810D+11,    -3.99096175224466498D+11,
     E     8.19218669548577329D+11,    -1.09837515608122331D+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1     C(105), C(106), C(107), C(108), C(109), C(110), C(111),
     2     C(112), C(113), C(114), C(115), C(116), C(117), C(118)/
     3     1.00815810686538209D+12,    -6.45364869245376503D+11,
     4     2.87900649906150589D+11,    -8.78670721780232657D+10,
     5     1.76347306068349694D+10,    -2.16716498322379509D+09,
     6     1.43157876718888981D+08,    -3.87183344257261262D+06,
     7     1.82577554742931747D+04,     2.86464035717679043D+11,
     8    -2.40629790002850396D+12,     9.10934118523989896D+12,
     9    -2.05168994109344374D+13,     3.05651255199353206D+13,
     A    -3.16670885847851584D+13,     2.33483640445818409D+13,
     B    -1.23204913055982872D+13,     4.61272578084913197D+12,
     C    -1.19655288019618160D+12,     2.05914503232410016D+11,
     D    -2.18229277575292237D+10,     1.24700929351271032D+09/
      DATA C(119), C(120)/
     1    -2.91883881222208134D+07,     1.18838426256783253D+05/
C***FIRST EXECUTABLE STATEMENT  ZUNIK
      IF (INIT.NE.0) GO TO 40
C-----------------------------------------------------------------------
C     INITIALIZE ALL VARIABLES
C-----------------------------------------------------------------------
      RFN = 1.0D0/FNU
C-----------------------------------------------------------------------
C     OVERFLOW TEST (ZR/FNU TOO SMALL)
C-----------------------------------------------------------------------
      TEST = D1MACH(1)*1.0D+3
      AC = FNU*TEST
      IF (ABS(ZRR).GT.AC .OR. ABS(ZRI).GT.AC) GO TO 15
      ZETA1R = 2.0D0*ABS(LOG(TEST))+FNU
      ZETA1I = 0.0D0
      ZETA2R = FNU
      ZETA2I = 0.0D0
      PHIR = 1.0D0
      PHII = 0.0D0
      RETURN
   15 CONTINUE
      TR = ZRR*RFN
      TI = ZRI*RFN
      SR = CONER + (TR*TR-TI*TI)
      SI = CONEI + (TR*TI+TI*TR)
      CALL ZSQRT(SR, SI, SRR, SRI)
      STR = CONER + SRR
      STI = CONEI + SRI
      CALL ZDIV(STR, STI, TR, TI, ZNR, ZNI)
      CALL ZLOG(ZNR, ZNI, STR, STI, IDUM)
      ZETA1R = FNU*STR
      ZETA1I = FNU*STI
      ZETA2R = FNU*SRR
      ZETA2I = FNU*SRI
      CALL ZDIV(CONER, CONEI, SRR, SRI, TR, TI)
      SRR = TR*RFN
      SRI = TI*RFN
      CALL ZSQRT(SRR, SRI, CWRKR(16), CWRKI(16))
      PHIR = CWRKR(16)*CON(IKFLG)
      PHII = CWRKI(16)*CON(IKFLG)
      IF (IPMTR.NE.0) RETURN
      CALL ZDIV(CONER, CONEI, SR, SI, T2R, T2I)
      CWRKR(1) = CONER
      CWRKI(1) = CONEI
      CRFNR = CONER
      CRFNI = CONEI
      AC = 1.0D0
      L = 1
      DO 20 K=2,15
        SR = ZEROR
        SI = ZEROI
        DO 10 J=1,K
          L = L + 1
          STR = SR*T2R - SI*T2I + C(L)
          SI = SR*T2I + SI*T2R
          SR = STR
   10   CONTINUE
        STR = CRFNR*SRR - CRFNI*SRI
        CRFNI = CRFNR*SRI + CRFNI*SRR
        CRFNR = STR
        CWRKR(K) = CRFNR*SR - CRFNI*SI
        CWRKI(K) = CRFNR*SI + CRFNI*SR
        AC = AC*RFN
        TEST = ABS(CWRKR(K)) + ABS(CWRKI(K))
        IF (AC.LT.TOL .AND. TEST.LT.TOL) GO TO 30
   20 CONTINUE
      K = 15
   30 CONTINUE
      INIT = K
   40 CONTINUE
      IF (IKFLG.EQ.2) GO TO 60
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE I FUNCTION
C-----------------------------------------------------------------------
      SR = ZEROR
      SI = ZEROI
      DO 50 I=1,INIT
        SR = SR + CWRKR(I)
        SI = SI + CWRKI(I)
   50 CONTINUE
      SUMR = SR
      SUMI = SI
      PHIR = CWRKR(16)*CON(1)
      PHII = CWRKI(16)*CON(1)
      RETURN
   60 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE K FUNCTION
C-----------------------------------------------------------------------
      SR = ZEROR
      SI = ZEROI
      TR = CONER
      DO 70 I=1,INIT
        SR = SR + TR*CWRKR(I)
        SI = SI + TR*CWRKI(I)
        TR = -TR
   70 CONTINUE
      SUMR = SR
      SUMI = SI
      PHIR = CWRKR(16)*CON(2)
      PHII = CWRKI(16)*CON(2)
      RETURN
      END

*DECK ZRATI
      SUBROUTINE ZRATI (ZR, ZI, FNU, N, CYR, CYI, TOL)
C***BEGIN PROLOGUE  ZRATI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CRATI-A, ZRATI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
C     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
C     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
C     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
C     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
C     BY D. J. SOOKNE.
C
C***SEE ALSO  ZBESH, ZBESI, ZBESK
C***ROUTINES CALLED  ZABS, ZDIV
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZRATI
      DOUBLE PRECISION AK, AMAGZ, AP1, AP2, ARG, AZ, CDFNUI, CDFNUR,
     * CONEI, CONER, CYI, CYR, CZEROI, CZEROR, DFNU, FDNU, FLAM, FNU,
     * FNUP, PTI, PTR, P1I, P1R, P2I, P2R, RAK, RAP1, RHO, RT2, RZI,
     * RZR, TEST, TEST1, TOL, TTI, TTR, T1I, T1R, ZI, ZR, ZABS
      INTEGER I, ID, IDNU, INU, ITIME, K, KK, MAGZ, N
      DIMENSION CYR(N), CYI(N)
      EXTERNAL ZABS
      DATA CZEROR,CZEROI,CONER,CONEI,RT2/
     1 0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.41421356237309505D0 /
C***FIRST EXECUTABLE STATEMENT  ZRATI
      AZ = ZABS(ZR,ZI)
      INU = FNU
      IDNU = INU + N - 1
      MAGZ = AZ
      AMAGZ = MAGZ+1
      FDNU = IDNU
      FNUP = MAX(AMAGZ,FDNU)
      ID = IDNU - MAGZ - 1
      ITIME = 1
      K = 1
      PTR = 1.0D0/AZ
      RZR = PTR*(ZR+ZR)*PTR
      RZI = -PTR*(ZI+ZI)*PTR
      T1R = RZR*FNUP
      T1I = RZI*FNUP
      P2R = -T1R
      P2I = -T1I
      P1R = CONER
      P1I = CONEI
      T1R = T1R + RZR
      T1I = T1I + RZI
      IF (ID.GT.0) ID = 0
      AP2 = ZABS(P2R,P2I)
      AP1 = ZABS(P1R,P1I)
C-----------------------------------------------------------------------
C     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
C     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
C     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
C     PREMATURELY.
C-----------------------------------------------------------------------
      ARG = (AP2+AP2)/(AP1*TOL)
      TEST1 = SQRT(ARG)
      TEST = TEST1
      RAP1 = 1.0D0/AP1
      P1R = P1R*RAP1
      P1I = P1I*RAP1
      P2R = P2R*RAP1
      P2I = P2I*RAP1
      AP2 = AP2*RAP1
   10 CONTINUE
      K = K + 1
      AP1 = AP2
      PTR = P2R
      PTI = P2I
      P2R = P1R - (T1R*PTR-T1I*PTI)
      P2I = P1I - (T1R*PTI+T1I*PTR)
      P1R = PTR
      P1I = PTI
      T1R = T1R + RZR
      T1I = T1I + RZI
      AP2 = ZABS(P2R,P2I)
      IF (AP1.LE.TEST) GO TO 10
      IF (ITIME.EQ.2) GO TO 20
      AK = ZABS(T1R,T1I)*0.5D0
      FLAM = AK + SQRT(AK*AK-1.0D0)
      RHO = MIN(AP2/AP1,FLAM)
      TEST = TEST1*SQRT(RHO/(RHO*RHO-1.0D0))
      ITIME = 2
      GO TO 10
   20 CONTINUE
      KK = K + 1 - ID
      AK = KK
      T1R = AK
      T1I = CZEROI
      DFNU = FNU + (N-1)
      P1R = 1.0D0/AP2
      P1I = CZEROI
      P2R = CZEROR
      P2I = CZEROI
      DO 30 I=1,KK
        PTR = P1R
        PTI = P1I
        RAP1 = DFNU + T1R
        TTR = RZR*RAP1
        TTI = RZI*RAP1
        P1R = (PTR*TTR-PTI*TTI) + P2R
        P1I = (PTR*TTI+PTI*TTR) + P2I
        P2R = PTR
        P2I = PTI
        T1R = T1R - CONER
   30 CONTINUE
      IF (P1R.NE.CZEROR .OR. P1I.NE.CZEROI) GO TO 40
      P1R = TOL
      P1I = TOL
   40 CONTINUE
      CALL ZDIV(P2R, P2I, P1R, P1I, CYR(N), CYI(N))
      IF (N.EQ.1) RETURN
      K = N - 1
      AK = K
      T1R = AK
      T1I = CZEROI
      CDFNUR = FNU*RZR
      CDFNUI = FNU*RZI
      DO 60 I=2,N
        PTR = CDFNUR + (T1R*RZR-T1I*RZI) + CYR(K+1)
        PTI = CDFNUI + (T1R*RZI+T1I*RZR) + CYI(K+1)
        AK = ZABS(PTR,PTI)
        IF (AK.NE.CZEROR) GO TO 50
        PTR = TOL
        PTI = TOL
        AK = TOL*RT2
   50   CONTINUE
        RAK = CONER/AK
        CYR(K) = RAK*PTR*RAK
        CYI(K) = -RAK*PTI*RAK
        T1R = T1R - CONER
        K = K - 1
   60 CONTINUE
      RETURN
      END

*DECK ZUNI1
      SUBROUTINE ZUNI1 (ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL,
     +   TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZUNI1
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUNI1-A, ZUNI1-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
C     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZUCHK, ZUNIK, ZUOIK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZUNI1
C     COMPLEX CFN,CONE,CRSC,CSCL,CSR,CSS,CWRK,CZERO,C1,C2,PHI,RZ,SUM,S1,
C    *S2,Y,Z,ZETA1,ZETA2
      DOUBLE PRECISION ALIM, APHI, ASCLE, BRY, CONER, CRSC,
     * CSCL, CSRR, CSSR, CWRKI, CWRKR, C1R, C2I, C2M, C2R, ELIM, FN,
     * FNU, FNUL, PHII, PHIR, RAST, RS1, RZI, RZR, STI, STR, SUMI,
     * SUMR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I,
     * ZETA1R, ZETA2I, ZETA2R, ZI, ZR, CYR, CYI, D1MACH, ZABS
      INTEGER I, IFLAG, INIT, K, KODE, M, N, ND, NLAST, NN, NUF, NW, NZ
      DIMENSION BRY(3), YR(N), YI(N), CWRKR(16), CWRKI(16), CSSR(3),
     * CSRR(3), CYR(2), CYI(2)
      EXTERNAL ZABS
      DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZUNI1
      NZ = 0
      ND = N
      NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
      FN = MAX(FNU,1.0D0)
      INIT = 0
      CALL ZUNIK(ZR, ZI, FN, 1, 1, TOL, INIT, PHIR, PHII, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
      IF (KODE.EQ.1) GO TO 10
      STR = ZR + ZETA2R
      STI = ZI + ZETA2I
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = -ZETA1R + STR
      S1I = -ZETA1I + STI
      GO TO 20
   10 CONTINUE
      S1R = -ZETA1R + ZETA2R
      S1I = -ZETA1I + ZETA2I
   20 CONTINUE
      RS1 = S1R
      IF (ABS(RS1).GT.ELIM) GO TO 130
   30 CONTINUE
      NN = MIN(2,ND)
      DO 80 I=1,NN
        FN = FNU + (ND-I)
        INIT = 0
        CALL ZUNIK(ZR, ZI, FN, 1, 0, TOL, INIT, PHIR, PHII, ZETA1R,
     *   ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
        IF (KODE.EQ.1) GO TO 40
        STR = ZR + ZETA2R
        STI = ZI + ZETA2I
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZETA1R + STR
        S1I = -ZETA1I + STI + ZI
        GO TO 50
   40   CONTINUE
        S1R = -ZETA1R + ZETA2R
        S1I = -ZETA1I + ZETA2I
   50   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (ABS(RS1).GT.ELIM) GO TO 110
        IF (I.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 60
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIR,PHII)
        RS1 = RS1 + LOG(APHI)
        IF (ABS(RS1).GT.ELIM) GO TO 110
        IF (I.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 60
        IF (I.EQ.1) IFLAG = 3
   60   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 IF ABS(S1).LT.ASCLE
C-----------------------------------------------------------------------
        S2R = PHIR*SUMR - PHII*SUMI
        S2I = PHIR*SUMI + PHII*SUMR
        STR = EXP(S1R)*CSSR(IFLAG)
        S1R = STR*COS(S1I)
        S1I = STR*SIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 70
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 110
   70   CONTINUE
        CYR(I) = S2R
        CYI(I) = S2I
        M = ND - I + 1
        YR(M) = S2R*CSRR(IFLAG)
        YI(M) = S2I*CSRR(IFLAG)
   80 CONTINUE
      IF (ND.LE.2) GO TO 100
      RAST = 1.0D0/ZABS(ZR,ZI)
      STR = ZR*RAST
      STI = -ZI*RAST
      RZR = (STR+STR)*RAST
      RZI = (STI+STI)*RAST
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = K
      DO 90 I=3,ND
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(K) = C2R
        YI(K) = C2I
        K = K - 1
        FN = FN - 1.0D0
        IF (IFLAG.GE.3) GO TO 90
        STR = ABS(C2R)
        STI = ABS(C2I)
        C2M = MAX(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 90
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        C1R = CSRR(IFLAG)
   90 CONTINUE
  100 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
  110 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 120
      YR(ND) = ZEROR
      YI(ND) = ZEROI
      NZ = NZ + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 100
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 120
      ND = ND - NUF
      NZ = NZ + NUF
      IF (ND.EQ.0) GO TO 100
      FN = FNU + (ND-1)
      IF (FN.GE.FNUL) GO TO 30
      NLAST = ND
      RETURN
  120 CONTINUE
      NZ = -1
      RETURN
  130 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 120
      NZ = N
      DO 140 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  140 CONTINUE
      RETURN
      END
*DECK ZUNI2
      SUBROUTINE ZUNI2 (ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL,
     +   TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZUNI2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUNI2-A, ZUNI2-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
C     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
C     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZAIRY, ZUCHK, ZUNHJ, ZUOIK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZUNI2
C     COMPLEX AI,ARG,ASUM,BSUM,CFN,CI,CID,CIP,CONE,CRSC,CSCL,CSR,CSS,
C    *CZERO,C1,C2,DAI,PHI,RZ,S1,S2,Y,Z,ZB,ZETA1,ZETA2,ZN
      DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGI,
     * ARGR, ASCLE, ASUMI, ASUMR, BRY, BSUMI, BSUMR, CIDI, CIPI, CIPR,
     * CONER, CRSC, CSCL, CSRR, CSSR, C1R, C2I, C2M, C2R, DAII,
     * DAIR, ELIM, FN, FNU, FNUL, HPI, PHII, PHIR, RAST, RAZ, RS1, RZI,
     * RZR, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZBI, ZBR, ZEROI,
     * ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI, ZNI, ZNR, ZR, CYR,
     * CYI, D1MACH, ZABS, CAR, SAR
      INTEGER I, IFLAG, IN, INU, J, K, KODE, N, NAI, ND, NDAI, NLAST,
     * NN, NUF, NW, NZ, IDUM
      DIMENSION BRY(3), YR(N), YI(N), CIPR(4), CIPI(4), CSSR(3),
     * CSRR(3), CYR(2), CYI(2)
      EXTERNAL ZABS
      DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
      DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4),
     * CIPI(4)/ 1.0D0,0.0D0, 0.0D0,1.0D0, -1.0D0,0.0D0, 0.0D0,-1.0D0/
      DATA HPI, AIC  /
     1      1.57079632679489662D+00,     1.265512123484645396D+00/
C***FIRST EXECUTABLE STATEMENT  ZUNI2
      NZ = 0
      ND = N
      NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
C-----------------------------------------------------------------------
C     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
C-----------------------------------------------------------------------
      ZNR = ZI
      ZNI = -ZR
      ZBR = ZR
      ZBI = ZI
      CIDI = -CONER
      INU = FNU
      ANG = HPI*(FNU-INU)
      C2R = COS(ANG)
      C2I = SIN(ANG)
      CAR = C2R
      SAR = C2I
      IN = INU + N - 1
      IN = MOD(IN,4) + 1
      STR = C2R*CIPR(IN) - C2I*CIPI(IN)
      C2I = C2R*CIPI(IN) + C2I*CIPR(IN)
      C2R = STR
      IF (ZI.GT.0.0D0) GO TO 10
      ZNR = -ZNR
      ZBI = -ZBI
      CIDI = -CIDI
      C2I = -C2I
   10 CONTINUE
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
      FN = MAX(FNU,1.0D0)
      CALL ZUNHJ(ZNR, ZNI, FN, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
      IF (KODE.EQ.1) GO TO 20
      STR = ZBR + ZETA2R
      STI = ZBI + ZETA2I
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = -ZETA1R + STR
      S1I = -ZETA1I + STI
      GO TO 30
   20 CONTINUE
      S1R = -ZETA1R + ZETA2R
      S1I = -ZETA1I + ZETA2I
   30 CONTINUE
      RS1 = S1R
      IF (ABS(RS1).GT.ELIM) GO TO 150
   40 CONTINUE
      NN = MIN(2,ND)
      DO 90 I=1,NN
        FN = FNU + (ND-I)
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR, PHII, ARGR, ARGI,
     *   ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
        IF (KODE.EQ.1) GO TO 50
        STR = ZBR + ZETA2R
        STI = ZBI + ZETA2I
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZETA1R + STR
        S1I = -ZETA1I + STI + ABS(ZI)
        GO TO 60
   50   CONTINUE
        S1R = -ZETA1R + ZETA2R
        S1I = -ZETA1I + ZETA2I
   60   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (ABS(RS1).GT.ELIM) GO TO 120
        IF (I.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 70
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        APHI = ZABS(PHIR,PHII)
        AARG = ZABS(ARGR,ARGI)
        RS1 = RS1 + LOG(APHI) - 0.25D0*LOG(AARG) - AIC
        IF (ABS(RS1).GT.ELIM) GO TO 120
        IF (I.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 70
        IF (I.EQ.1) IFLAG = 3
   70   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        CALL ZAIRY(ARGR, ARGI, 0, 2, AIR, AII, NAI, IDUM)
        CALL ZAIRY(ARGR, ARGI, 1, 2, DAIR, DAII, NDAI, IDUM)
        STR = DAIR*BSUMR - DAII*BSUMI
        STI = DAIR*BSUMI + DAII*BSUMR
        STR = STR + (AIR*ASUMR-AII*ASUMI)
        STI = STI + (AIR*ASUMI+AII*ASUMR)
        S2R = PHIR*STR - PHII*STI
        S2I = PHIR*STI + PHII*STR
        STR = EXP(S1R)*CSSR(IFLAG)
        S1R = STR*COS(S1I)
        S1I = STR*SIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 80
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 120
   80   CONTINUE
        IF (ZI.LE.0.0D0) S2I = -S2I
        STR = S2R*C2R - S2I*C2I
        S2I = S2R*C2I + S2I*C2R
        S2R = STR
        CYR(I) = S2R
        CYI(I) = S2I
        J = ND - I + 1
        YR(J) = S2R*CSRR(IFLAG)
        YI(J) = S2I*CSRR(IFLAG)
        STR = -C2I*CIDI
        C2I = C2R*CIDI
        C2R = STR
   90 CONTINUE
      IF (ND.LE.2) GO TO 110
      RAZ = 1.0D0/ZABS(ZR,ZI)
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = K
      DO 100 I=3,ND
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(K) = C2R
        YI(K) = C2I
        K = K - 1
        FN = FN - 1.0D0
        IF (IFLAG.GE.3) GO TO 100
        STR = ABS(C2R)
        STI = ABS(C2I)
        C2M = MAX(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 100
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        C1R = CSRR(IFLAG)
  100 CONTINUE
  110 CONTINUE
      RETURN
  120 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 140
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
      YR(ND) = ZEROR
      YI(ND) = ZEROI
      NZ = NZ + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 110
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 140
      ND = ND - NUF
      NZ = NZ + NUF
      IF (ND.EQ.0) GO TO 110
      FN = FNU + (ND-1)
      IF (FN.LT.FNUL) GO TO 130
C      FN = CIDI
C      J = NUF + 1
C      K = MOD(J,4) + 1
C      S1R = CIPR(K)
C      S1I = CIPI(K)
C      IF (FN.LT.0.0D0) S1I = -S1I
C      STR = C2R*S1R - C2I*S1I
C      C2I = C2R*S1I + C2I*S1R
C      C2R = STR
      IN = INU + ND - 1
      IN = MOD(IN,4) + 1
      C2R = CAR*CIPR(IN) - SAR*CIPI(IN)
      C2I = CAR*CIPI(IN) + SAR*CIPR(IN)
      IF (ZI.LE.0.0D0) C2I = -C2I
      GO TO 40
  130 CONTINUE
      NLAST = ND
      RETURN
  140 CONTINUE
      NZ = -1
      RETURN
  150 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 140
      NZ = N
      DO 160 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  160 CONTINUE
      RETURN
      END

*DECK ZSHCH
      SUBROUTINE ZSHCH (ZR, ZI, CSHR, CSHI, CCHR, CCHI)
C***BEGIN PROLOGUE  ZSHCH
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CSHCH-A, ZSHCH-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
C     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
C
C***SEE ALSO  ZBESH, ZBESK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZSHCH
C
      DOUBLE PRECISION CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, ZI, ZR
C***FIRST EXECUTABLE STATEMENT  ZSHCH
      SH = SINH(ZR)
      CH = COSH(ZR)
      SN = SIN(ZI)
      CN = COS(ZI)
      CSHR = SH*CN
      CSHI = CH*SN
      CCHR = CH*CN
      CCHI = SH*SN
      RETURN
      END


*DECK ZKSCL
      SUBROUTINE ZKSCL (ZRR, ZRI, FNU, N, YR, YI, NZ, RZR, RZI, ASCLE,
     +   TOL, ELIM)
C***BEGIN PROLOGUE  ZKSCL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CKSCL-A, ZKSCL-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
C     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
C     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
C
C***SEE ALSO  ZBESK
C***ROUTINES CALLED  ZABS, ZLOG, ZUCHK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZLOG to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZKSCL
C     COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM
      DOUBLE PRECISION ACS, AS, ASCLE, CKI, CKR, CSI, CSR, CYI,
     * CYR, ELIM, FN, FNU, RZI, RZR, STR, S1I, S1R, S2I,
     * S2R, TOL, YI, YR, ZEROI, ZEROR, ZRI, ZRR, ZABS,
     * ZDR, ZDI, CELMR, ELM, HELIM, ALAS
      INTEGER I, IC, IDUM, KK, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2)
      EXTERNAL ZABS, ZLOG
      DATA ZEROR,ZEROI / 0.0D0 , 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZKSCL
      NZ = 0
      IC = 0
      NN = MIN(2,N)
      DO 10 I=1,NN
        S1R = YR(I)
        S1I = YI(I)
        CYR(I) = S1R
        CYI(I) = S1I
        AS = ZABS(S1R,S1I)
        ACS = -ZRR + LOG(AS)
        NZ = NZ + 1
        YR(I) = ZEROR
        YI(I) = ZEROI
        IF (ACS.LT.(-ELIM)) GO TO 10
        CALL ZLOG(S1R, S1I, CSR, CSI, IDUM)
        CSR = CSR - ZRR
        CSI = CSI - ZRI
        STR = EXP(CSR)/TOL
        CSR = STR*COS(CSI)
        CSI = STR*SIN(CSI)
        CALL ZUCHK(CSR, CSI, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 10
        YR(I) = CSR
        YI(I) = CSI
        IC = I
        NZ = NZ - 1
   10 CONTINUE
      IF (N.EQ.1) RETURN
      IF (IC.GT.1) GO TO 20
      YR(1) = ZEROR
      YI(1) = ZEROI
      NZ = 2
   20 CONTINUE
      IF (N.EQ.2) RETURN
      IF (NZ.EQ.0) RETURN
      FN = FNU + 1.0D0
      CKR = FN*RZR
      CKI = FN*RZI
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      HELIM = 0.5D0*ELIM
      ELM = EXP(-ELIM)
      CELMR = ELM
      ZDR = ZRR
      ZDI = ZRI
C
C     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
C     S2 GETS LARGER THAN EXP(ELIM/2)
C
      DO 30 I=3,N
        KK = I
        CSR = S2R
        CSI = S2I
        S2R = CKR*CSR - CKI*CSI + S1R
        S2I = CKI*CSR + CKR*CSI + S1I
        S1R = CSR
        S1I = CSI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AS = ZABS(S2R,S2I)
        ALAS = LOG(AS)
        ACS = -ZDR + ALAS
        NZ = NZ + 1
        YR(I) = ZEROR
        YI(I) = ZEROI
        IF (ACS.LT.(-ELIM)) GO TO 25
        CALL ZLOG(S2R, S2I, CSR, CSI, IDUM)
        CSR = CSR - ZDR
        CSI = CSI - ZDI
        STR = EXP(CSR)/TOL
        CSR = STR*COS(CSI)
        CSI = STR*SIN(CSI)
        CALL ZUCHK(CSR, CSI, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 25
        YR(I) = CSR
        YI(I) = CSI
        NZ = NZ - 1
        IF (IC.EQ.KK-1) GO TO 40
        IC = KK
        GO TO 30
   25   CONTINUE
        IF(ALAS.LT.HELIM) GO TO 30
        ZDR = ZDR - ELIM
        S1R = S1R*CELMR
        S1I = S1I*CELMR
        S2R = S2R*CELMR
        S2I = S2I*CELMR
   30 CONTINUE
      NZ = N
      IF(IC.EQ.N) NZ=N-1
      GO TO 45
   40 CONTINUE
      NZ = KK - 2
   45 CONTINUE
      DO 50 I=1,NZ
        YR(I) = ZEROR
        YI(I) = ZEROI
   50 CONTINUE
      RETURN
      END

*DECK ZS1S2
      SUBROUTINE ZS1S2 (ZRR, ZRI, S1R, S1I, S2R, S2I, NZ, ASCLE, ALIM,
     +   IUF)
C***BEGIN PROLOGUE  ZS1S2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZAIRY and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CS1S2-A, ZS1S2-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
C     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
C     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
C     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
C     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
C     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
C     PRECISION ABOVE THE UNDERFLOW LIMIT.
C
C***SEE ALSO  ZAIRY, ZBESK
C***ROUTINES CALLED  ZABS, ZEXP, ZLOG
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZS1S2
C     COMPLEX CZERO,C1,S1,S1D,S2,ZR
      DOUBLE PRECISION AA, ALIM, ALN, ASCLE, AS1, AS2, C1I, C1R, S1DI,
     * S1DR, S1I, S1R, S2I, S2R, ZEROI, ZEROR, ZRI, ZRR, ZABS
      INTEGER IUF, IDUM, NZ
      EXTERNAL ZABS, ZEXP, ZLOG
      DATA ZEROR,ZEROI  / 0.0D0 , 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZS1S2
      NZ = 0
      AS1 = ZABS(S1R,S1I)
      AS2 = ZABS(S2R,S2I)
      IF (S1R.EQ.0.0D0 .AND. S1I.EQ.0.0D0) GO TO 10
      IF (AS1.EQ.0.0D0) GO TO 10
      ALN = -ZRR - ZRR + LOG(AS1)
      S1DR = S1R
      S1DI = S1I
      S1R = ZEROR
      S1I = ZEROI
      AS1 = ZEROR
      IF (ALN.LT.(-ALIM)) GO TO 10
      CALL ZLOG(S1DR, S1DI, C1R, C1I, IDUM)
      C1R = C1R - ZRR - ZRR
      C1I = C1I - ZRI - ZRI
      CALL ZEXP(C1R, C1I, S1R, S1I)
      AS1 = ZABS(S1R,S1I)
      IUF = IUF + 1
   10 CONTINUE
      AA = MAX(AS1,AS2)
      IF (AA.GT.ASCLE) RETURN
      S1R = ZEROR
      S1I = ZEROI
      S2R = ZEROR
      S2I = ZEROI
      NZ = 1
      IUF = 0
      RETURN
      END

*DECK ZAIRY
      SUBROUTINE ZAIRY (ZR, ZI, ID, KODE, AIR, AII, NZ, IERR)
C***BEGIN PROLOGUE  ZAIRY
C***PURPOSE  Compute the Airy function Ai(z) or its derivative dAi/dz
C            for complex argument z.  A scaling option is available
C            to help avoid underflow and overflow.
C***LIBRARY   SLATEC
C***CATEGORY  C10D
C***TYPE      COMPLEX (CAIRY-C, ZAIRY-C)
C***KEYWORDS  AIRY FUNCTION, BESSEL FUNCTION OF ORDER ONE THIRD,
C             BESSEL FUNCTION OF ORDER TWO THIRDS
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         On KODE=1, ZAIRY computes the complex Airy function Ai(z)
C         or its derivative dAi/dz on ID=0 or ID=1 respectively. On
C         KODE=2, a scaling option exp(zeta)*Ai(z) or exp(zeta)*dAi/dz
C         is provided to remove the exponential decay in -pi/3<arg(z)
C         <pi/3 and the exponential growth in pi/3<abs(arg(z))<pi where
C         zeta=(2/3)*z**(3/2).
C
C         While the Airy functions Ai(z) and dAi/dz are analytic in
C         the whole z-plane, the corresponding scaled functions defined
C         for KODE=2 have a cut along the negative real axis.
C
C         Input
C           ZR     - DOUBLE PRECISION real part of argument Z
C           ZI     - DOUBLE PRECISION imag part of argument Z
C           ID     - Order of derivative, ID=0 or ID=1
C           KODE   - A parameter to indicate the scaling option
C                    KODE=1  returns
C                            AI=Ai(z)  on ID=0
C                            AI=dAi/dz on ID=1
C                            at z=Z
C                        =2  returns
C                            AI=exp(zeta)*Ai(z)  on ID=0
C                            AI=exp(zeta)*dAi/dz on ID=1
C                            at z=Z where zeta=(2/3)*z**(3/2)
C
C         Output
C           AIR    - DOUBLE PRECISION real part of result
C           AII    - DOUBLE PRECISION imag part of result
C           NZ     - Underflow indicator
C                    NZ=0    Normal return
C                    NZ=1    AI=0 due to underflow in
C                            -pi/3<arg(Z)<pi/3 on KODE=1
C           IERR   - Error flag
C                    IERR=0  Normal return     - COMPUTATION COMPLETED
C                    IERR=1  Input error       - NO COMPUTATION
C                    IERR=2  Overflow          - NO COMPUTATION
C                            (Re(Z) too large with KODE=1)
C                    IERR=3  Precision warning - COMPUTATION COMPLETED
C                            (Result has less than half precision)
C                    IERR=4  Precision error   - NO COMPUTATION
C                            (Result has no precision)
C                    IERR=5  Algorithmic error - NO COMPUTATION
C                            (Termination condition not met)
C
C *Long Description:
C
C         Ai(z) and dAi/dz are computed from K Bessel functions by
C
C                Ai(z) =  c*sqrt(z)*K(1/3,zeta)
C               dAi/dz = -c*   z   *K(2/3,zeta)
C                    c =  1/(pi*sqrt(3))
C                 zeta =  (2/3)*z**(3/2)
C
C         when abs(z)>1 and from power series when abs(z)<=1.
C
C         In most complex variable computation, one must evaluate ele-
C         mentary functions.  When the magnitude of Z is large, losses
C         of significance by argument reduction occur.  Consequently, if
C         the magnitude of ZETA=(2/3)*Z**(3/2) exceeds U1=SQRT(0.5/UR),
C         then losses exceeding half precision are likely and an error
C         flag IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is
C         double precision unit roundoff limited to 18 digits precision.
C         Also, if the magnitude of ZETA is larger than U2=0.5/UR, then
C         all significance is lost and IERR=4.  In order to use the INT
C         function, ZETA must be further restricted not to exceed
C         U3=I1MACH(9)=LARGEST INTEGER.  Thus, the magnitude of ZETA
C         must be restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2,
C         and U3 are approximately 2.0E+3, 4.2E+6, 2.1E+9 in single
C         precision and 4.7E+7, 2.3E+15, 2.1E+9 in double precision.
C         This makes U2 limiting is single precision and U3 limiting
C         in double precision.  This means that the magnitude of Z
C         cannot exceed approximately 3.4E+4 in single precision and
C         2.1E+6 in double precision.  This also means that one can
C         expect to retain, in the worst cases on 32-bit machines,
C         no digits in single precision and only 6 digits in double
C         precision.
C
C         The approximate relative error in the magnitude of a complex
C         Bessel function can be expressed as P*10**S where P=MAX(UNIT
C         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre-
C         sents the increase in error due to argument reduction in the
C         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))),
C         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF
C         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may
C         have only absolute accuracy.  This is most likely to occur
C         when one component (in magnitude) is larger than the other by
C         several orders of magnitude.  If one component is 10**K larger
C         than the other, then one can expect only MAX(ABS(LOG10(P))-K,
C         0) significant digits; or, stated another way, when K exceeds
C         the exponent of P, no significant digits remain in the smaller
C         component.  However, the phase angle retains absolute accuracy
C         because, in complex arithmetic with precision P, the smaller
C         component will not (as a rule) decrease below P times the
C         magnitude of the larger component. In these extreme cases,
C         the principal phase angle is on the order of +P, -P, PI/2-P,
C         or -PI/2+P.
C
C***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
C                 matical Functions, National Bureau of Standards
C                 Applied Mathematics Series 55, U. S. Department
C                 of Commerce, Tenth Printing (1972) or later.
C               2. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument and Large Order, Report SAND83-0643,
C                 Sandia National Laboratories, Albuquerque, NM, May
C                 1983.
C               3. D. E. Amos, A Subroutine Package for Bessel Functions
C                 of a Complex Argument and Nonnegative Order, Report
C                 SAND85-1018, Sandia National Laboratory, Albuquerque,
C                 NM, May 1985.
C               4. D. E. Amos, A portable package for Bessel functions
C                 of a complex argument and nonnegative order, ACM
C                 Transactions on Mathematical Software, 12 (September
C                 1986), pp. 265-273.
C
C***ROUTINES CALLED  D1MACH, I1MACH, ZABS, ZACAI, ZBKNU, ZEXP, ZSQRT
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   890801  REVISION DATE from Version 3.2
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   920128  Category corrected.  (WRB)
C   920811  Prologue revised.  (DWL)
C   930122  Added ZEXP and ZSQRT to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZAIRY
C     COMPLEX AI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3
      DOUBLE PRECISION AA, AD, AII, AIR, AK, ALIM, ATRM, AZ, AZ3, BK,
     * CC, CK, COEF, CONEI, CONER, CSQI, CSQR, CYI, CYR, C1, C2, DIG,
     * DK, D1, D2, ELIM, FID, FNU, PTR, RL, R1M5, SFAC, STI, STR,
     * S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I, TRM2R, TTH, ZEROI,
     * ZEROR, ZI, ZR, ZTAI, ZTAR, Z3I, Z3R, D1MACH, ZABS, ALAZ, BB
      INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, I1MACH
      DIMENSION CYR(1), CYI(1)
      EXTERNAL ZABS, ZEXP, ZSQRT
      DATA TTH, C1, C2, COEF /6.66666666666666667D-01,
     * 3.55028053887817240D-01,2.58819403792806799D-01,
     * 1.83776298473930683D-01/
      DATA ZEROR, ZEROI, CONER, CONEI /0.0D0,0.0D0,1.0D0,0.0D0/
C***FIRST EXECUTABLE STATEMENT  ZAIRY
      IERR = 0
      NZ=0
      IF (ID.LT.0 .OR. ID.GT.1) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (IERR.NE.0) RETURN
      AZ = ZABS(ZR,ZI)
      TOL = MAX(D1MACH(4),1.0D-18)
      FID = ID
      IF (AZ.GT.1.0D0) GO TO 70
C-----------------------------------------------------------------------
C     POWER SERIES FOR ABS(Z).LE.1.
C-----------------------------------------------------------------------
      S1R = CONER
      S1I = CONEI
      S2R = CONER
      S2I = CONEI
      IF (AZ.LT.TOL) GO TO 170
      AA = AZ*AZ
      IF (AA.LT.TOL/AZ) GO TO 40
      TRM1R = CONER
      TRM1I = CONEI
      TRM2R = CONER
      TRM2I = CONEI
      ATRM = 1.0D0
      STR = ZR*ZR - ZI*ZI
      STI = ZR*ZI + ZI*ZR
      Z3R = STR*ZR - STI*ZI
      Z3I = STR*ZI + STI*ZR
      AZ3 = AZ*AA
      AK = 2.0D0 + FID
      BK = 3.0D0 - FID - FID
      CK = 4.0D0 - FID
      DK = 3.0D0 + FID + FID
      D1 = AK*DK
      D2 = BK*CK
      AD = MIN(D1,D2)
      AK = 24.0D0 + 9.0D0*FID
      BK = 30.0D0 - 9.0D0*FID
      DO 30 K=1,25
        STR = (TRM1R*Z3R-TRM1I*Z3I)/D1
        TRM1I = (TRM1R*Z3I+TRM1I*Z3R)/D1
        TRM1R = STR
        S1R = S1R + TRM1R
        S1I = S1I + TRM1I
        STR = (TRM2R*Z3R-TRM2I*Z3I)/D2
        TRM2I = (TRM2R*Z3I+TRM2I*Z3R)/D2
        TRM2R = STR
        S2R = S2R + TRM2R
        S2I = S2I + TRM2I
        ATRM = ATRM*AZ3/AD
        D1 = D1 + AK
        D2 = D2 + BK
        AD = MIN(D1,D2)
        IF (ATRM.LT.TOL*AD) GO TO 40
        AK = AK + 18.0D0
        BK = BK + 18.0D0
   30 CONTINUE
   40 CONTINUE
      IF (ID.EQ.1) GO TO 50
      AIR = S1R*C1 - C2*(ZR*S2R-ZI*S2I)
      AII = S1I*C1 - C2*(ZR*S2I+ZI*S2R)
      IF (KODE.EQ.1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      CALL ZEXP(ZTAR, ZTAI, STR, STI)
      PTR = AIR*STR - AII*STI
      AII = AIR*STI + AII*STR
      AIR = PTR
      RETURN
   50 CONTINUE
      AIR = -S2R*C2
      AII = -S2I*C2
      IF (AZ.LE.TOL) GO TO 60
      STR = ZR*S1R - ZI*S1I
      STI = ZR*S1I + ZI*S1R
      CC = C1/(1.0D0+FID)
      AIR = AIR + CC*(STR*ZR-STI*ZI)
      AII = AII + CC*(STR*ZI+STI*ZR)
   60 CONTINUE
      IF (KODE.EQ.1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      CALL ZEXP(ZTAR, ZTAI, STR, STI)
      PTR = STR*AIR - STI*AII
      AII = STR*AII + STI*AIR
      AIR = PTR
      RETURN
C-----------------------------------------------------------------------
C     CASE FOR ABS(Z).GT.1.0
C-----------------------------------------------------------------------
   70 CONTINUE
      FNU = (1.0D0+FID)/3.0D0
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0D-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C-----------------------------------------------------------------------
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN(ABS(K1),ABS(K2))
      ELIM = 2.303D0*(K*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*K1
      DIG = MIN(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + MAX(-AA,-41.45D0)
      RL = 1.2D0*DIG + 3.0D0
      ALAZ = LOG(AZ)
C-----------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
      AA=0.5D0/TOL
      BB=I1MACH(9)*0.5D0
      AA=MIN(AA,BB)
      AA=AA**TTH
      IF (AZ.GT.AA) GO TO 260
      AA=SQRT(AA)
      IF (AZ.GT.AA) IERR=3
      CALL ZSQRT(ZR, ZI, CSQR, CSQI)
      ZTAR = TTH*(ZR*CSQR-ZI*CSQI)
      ZTAI = TTH*(ZR*CSQI+ZI*CSQR)
C-----------------------------------------------------------------------
C     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
C-----------------------------------------------------------------------
      IFLAG = 0
      SFAC = 1.0D0
      AK = ZTAI
      IF (ZR.GE.0.0D0) GO TO 80
      BK = ZTAR
      CK = -ABS(BK)
      ZTAR = CK
      ZTAI = AK
   80 CONTINUE
      IF (ZI.NE.0.0D0) GO TO 90
      IF (ZR.GT.0.0D0) GO TO 90
      ZTAR = 0.0D0
      ZTAI = AK
   90 CONTINUE
      AA = ZTAR
      IF (AA.GE.0.0D0 .AND. ZR.GT.0.0D0) GO TO 110
      IF (KODE.EQ.2) GO TO 100
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.GT.(-ALIM)) GO TO 100
      AA = -AA + 0.25D0*ALAZ
      IFLAG = 1
      SFAC = TOL
      IF (AA.GT.ELIM) GO TO 270
  100 CONTINUE
C-----------------------------------------------------------------------
C     CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
C-----------------------------------------------------------------------
      MR = 1
      IF (ZI.LT.0.0D0) MR = -1
      CALL ZACAI(ZTAR, ZTAI, FNU, KODE, MR, 1, CYR, CYI, NN, RL, TOL,
     * ELIM, ALIM)
      IF (NN.LT.0) GO TO 280
      NZ = NZ + NN
      GO TO 130
  110 CONTINUE
      IF (KODE.EQ.2) GO TO 120
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.LT.ALIM) GO TO 120
      AA = -AA - 0.25D0*ALAZ
      IFLAG = 2
      SFAC = 1.0D0/TOL
      IF (AA.LT.(-ELIM)) GO TO 210
  120 CONTINUE
      CALL ZBKNU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, TOL, ELIM,
     * ALIM)
  130 CONTINUE
      S1R = CYR(1)*COEF
      S1I = CYI(1)*COEF
      IF (IFLAG.NE.0) GO TO 150
      IF (ID.EQ.1) GO TO 140
      AIR = CSQR*S1R - CSQI*S1I
      AII = CSQR*S1I + CSQI*S1R
      RETURN
  140 CONTINUE
      AIR = -(ZR*S1R-ZI*S1I)
      AII = -(ZR*S1I+ZI*S1R)
      RETURN
  150 CONTINUE
      S1R = S1R*SFAC
      S1I = S1I*SFAC
      IF (ID.EQ.1) GO TO 160
      STR = S1R*CSQR - S1I*CSQI
      S1I = S1R*CSQI + S1I*CSQR
      S1R = STR
      AIR = S1R/SFAC
      AII = S1I/SFAC
      RETURN
  160 CONTINUE
      STR = -(S1R*ZR-S1I*ZI)
      S1I = -(S1R*ZI+S1I*ZR)
      S1R = STR
      AIR = S1R/SFAC
      AII = S1I/SFAC
      RETURN
  170 CONTINUE
      AA = 1.0D+3*D1MACH(1)
      S1R = ZEROR
      S1I = ZEROI
      IF (ID.EQ.1) GO TO 190
      IF (AZ.LE.AA) GO TO 180
      S1R = C2*ZR
      S1I = C2*ZI
  180 CONTINUE
      AIR = C1 - S1R
      AII = -S1I
      RETURN
  190 CONTINUE
      AIR = -C2
      AII = 0.0D0
      AA = SQRT(AA)
      IF (AZ.LE.AA) GO TO 200
      S1R = 0.5D0*(ZR*ZR-ZI*ZI)
      S1I = ZR*ZI
  200 CONTINUE
      AIR = AIR + C1*S1R
      AII = AII + C1*S1I
      RETURN
  210 CONTINUE
      NZ = 1
      AIR = ZEROR
      AII = ZEROI
      RETURN
  270 CONTINUE
      NZ = 0
      IERR=2
      RETURN
  280 CONTINUE
      IF(NN.EQ.(-1)) GO TO 270
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      IERR=4
      NZ=0
      RETURN
      END


*DECK ZACAI
      SUBROUTINE ZACAI (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, TOL,
     +   ELIM, ALIM)
C***BEGIN PROLOGUE  ZACAI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZAIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (CACAI-A, ZACAI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
C     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
C     RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT IF ZACON
C     IS CALLED FROM ZAIRY.
C
C***SEE ALSO  ZAIRY
C***ROUTINES CALLED  D1MACH, ZABS, ZASYI, ZBKNU, ZMLRI, ZS1S2, ZSERI
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZACAI
C     COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY
      DOUBLE PRECISION ALIM, ARG, ASCLE, AZ, CSGNR, CSGNI, CSPNR,
     * CSPNI, C1R, C1I, C2R, C2I, CYR, CYI, DFNU, ELIM, FMR, FNU, PI,
     * RL, SGN, TOL, YY, YR, YI, ZR, ZI, ZNR, ZNI, D1MACH, ZABS
      INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2)
      EXTERNAL ZABS
      DATA PI / 3.14159265358979324D0 /
C***FIRST EXECUTABLE STATEMENT  ZACAI
      NZ = 0
      ZNR = -ZR
      ZNI = -ZI
      AZ = ZABS(ZR,ZI)
      NN = N
      DFNU = FNU + (N-1)
      IF (AZ.LE.2.0D0) GO TO 10
      IF (AZ*AZ*0.25D0.GT.DFNU+1.0D0) GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
cSAP091104
c      write(*,*)'in ZACAI ZNR, ZNI, FNU',ZNR, ZNI, FNU

      CALL ZSERI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL, ELIM, ALIM)
      GO TO 40
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 30
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL ZASYI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 80
      GO TO 40
   30 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL ZMLRI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL)
      IF(NW.LT.0) GO TO 80
   40 CONTINUE
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, 1, CYR, CYI, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 80
      FMR = MR
      SGN = -DSIGN(PI,FMR)
      CSGNR = 0.0D0
      CSGNI = SGN
      IF (KODE.EQ.1) GO TO 50
      YY = -ZNI
      CSGNR = -CSGNI*SIN(YY)
      CSGNI = CSGNI*COS(YY)
   50 CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = FNU
      ARG = (FNU-INU)*SGN
      CSPNR = COS(ARG)
      CSPNI = SIN(ARG)
      IF (MOD(INU,2).EQ.0) GO TO 60
      CSPNR = -CSPNR
      CSPNI = -CSPNI
   60 CONTINUE
      C1R = CYR(1)
      C1I = CYI(1)
      C2R = YR(1)
      C2I = YI(1)
      IF (KODE.EQ.1) GO TO 70
      IUF = 0
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
   70 CONTINUE
      YR(1) = CSPNR*C1R - CSPNI*C1I + CSGNR*C2R - CSGNI*C2I
      YI(1) = CSPNR*C1I + CSPNI*C1R + CSGNR*C2I + CSGNI*C2R
      RETURN
   80 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END







*DECK DBESK0
      DOUBLE PRECISION FUNCTION DBESK0 (X)
C***BEGIN PROLOGUE  DBESK0
C***PURPOSE  Compute the modified (hyperbolic) Bessel function of the
C            third kind of order zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (BESK0-S, DBESK0-D)
C***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESK0(X) calculates the double precision modified (hyperbolic)
C Bessel function of the third kind of order zero for double
C precision argument X.  The argument must be greater than zero
C but not so large that the result underflows.
C
C Series for BK0        on the interval  0.          to  4.00000E+00
C                                        with weighted error   3.08E-33
C                                         log weighted error  32.51
C                               significant figures required  32.05
C                                    decimal places required  33.11
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DBESI0, DBSK0E, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBESK0
      DOUBLE PRECISION X, BK0CS(16), XMAX, XMAXT, XSML, Y,
     1  D1MACH, DCSEVL, DBESI0, DBSK0E
      LOGICAL FIRST
      SAVE BK0CS, NTK0, XSML, XMAX, FIRST
      DATA BK0CS(  1) / -.3532739323 3902768720 1140060063 153 D-1    /
      DATA BK0CS(  2) / +.3442898999 2462848688 6344927529 213 D+0    /
      DATA BK0CS(  3) / +.3597993651 5361501626 5721303687 231 D-1    /
      DATA BK0CS(  4) / +.1264615411 4469259233 8479508673 447 D-2    /
      DATA BK0CS(  5) / +.2286212103 1194517860 8269830297 585 D-4    /
      DATA BK0CS(  6) / +.2534791079 0261494573 0790013428 354 D-6    /
      DATA BK0CS(  7) / +.1904516377 2202088589 7214059381 366 D-8    /
      DATA BK0CS(  8) / +.1034969525 7633624585 1008317853 089 D-10   /
      DATA BK0CS(  9) / +.4259816142 7910825765 2445327170 133 D-13   /
      DATA BK0CS( 10) / +.1374465435 8807508969 4238325440 000 D-15   /
      DATA BK0CS( 11) / +.3570896528 5083735909 9688597333 333 D-18   /
      DATA BK0CS( 12) / +.7631643660 1164373766 7498666666 666 D-21   /
      DATA BK0CS( 13) / +.1365424988 4407818590 8053333333 333 D-23   /
      DATA BK0CS( 14) / +.2075275266 9066680831 9999999999 999 D-26   /
      DATA BK0CS( 15) / +.2712814218 0729856000 0000000000 000 D-29   /
      DATA BK0CS( 16) / +.3082593887 9146666666 6666666666 666 D-32   /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESK0
      IF (FIRST) THEN
         NTK0 = INITDS (BK0CS, 16, 0.1*REAL(D1MACH(3)))
         XSML = SQRT(4.0D0*D1MACH(3))
         XMAXT = -LOG(D1MACH(1))
         XMAX = XMAXT - 0.5D0*XMAXT*LOG(XMAXT)/(XMAXT+0.5D0)
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0.D0) CALL XERMSG ('SLATEC', 'DBESK0',
     +   'X IS ZERO OR NEGATIVE', 2, 2)
      IF (X.GT.2.0D0) GO TO 20
C
      Y = 0.D0
      IF (X.GT.XSML) Y = X*X
      DBESK0 = -LOG(0.5D0*X)*DBESI0(X) - 0.25D0 + DCSEVL (.5D0*Y-1.D0,
     1  BK0CS, NTK0)
      RETURN
C
 20   DBESK0 = 0.D0
      IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'DBESK0',
     +   'X SO BIG K0 UNDERFLOWS', 1, 1)
      IF (X.GT.XMAX) RETURN
C
      DBESK0 = EXP(-X) * DBSK0E(X)
C
      RETURN
      END
*DECK DBSK0E
      DOUBLE PRECISION FUNCTION DBSK0E (X)
C***BEGIN PROLOGUE  DBSK0E
C***PURPOSE  Compute the exponentially scaled modified (hyperbolic)
C            Bessel function of the third kind of order zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (BESK0E-S, DBSK0E-D)
C***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBSK0E(X) computes the double precision exponentially scaled
C modified (hyperbolic) Bessel function of the third kind of
C order zero for positive double precision argument X.
C
C Series for BK0        on the interval  0.          to  4.00000E+00
C                                        with weighted error   3.08E-33
C                                         log weighted error  32.51
C                               significant figures required  32.05
C                                    decimal places required  33.11
C
C Series for AK0        on the interval  1.25000E-01 to  5.00000E-01
C                                        with weighted error   2.85E-32
C                                         log weighted error  31.54
C                               significant figures required  30.19
C                                    decimal places required  32.33
C
C Series for AK02       on the interval  0.          to  1.25000E-01
C                                        with weighted error   2.30E-32
C                                         log weighted error  31.64
C                               significant figures required  29.68
C                                    decimal places required  32.40
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DBESI0, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBSK0E
      DOUBLE PRECISION X, BK0CS(16), AK0CS(38), AK02CS(33),
     1  XSML, Y, D1MACH, DCSEVL, DBESI0
      LOGICAL FIRST
      SAVE BK0CS, AK0CS, AK02CS, NTK0, NTAK0, NTAK02, XSML, FIRST
      DATA BK0CS(  1) / -.3532739323 3902768720 1140060063 153 D-1    /
      DATA BK0CS(  2) / +.3442898999 2462848688 6344927529 213 D+0    /
      DATA BK0CS(  3) / +.3597993651 5361501626 5721303687 231 D-1    /
      DATA BK0CS(  4) / +.1264615411 4469259233 8479508673 447 D-2    /
      DATA BK0CS(  5) / +.2286212103 1194517860 8269830297 585 D-4    /
      DATA BK0CS(  6) / +.2534791079 0261494573 0790013428 354 D-6    /
      DATA BK0CS(  7) / +.1904516377 2202088589 7214059381 366 D-8    /
      DATA BK0CS(  8) / +.1034969525 7633624585 1008317853 089 D-10   /
      DATA BK0CS(  9) / +.4259816142 7910825765 2445327170 133 D-13   /
      DATA BK0CS( 10) / +.1374465435 8807508969 4238325440 000 D-15   /
      DATA BK0CS( 11) / +.3570896528 5083735909 9688597333 333 D-18   /
      DATA BK0CS( 12) / +.7631643660 1164373766 7498666666 666 D-21   /
      DATA BK0CS( 13) / +.1365424988 4407818590 8053333333 333 D-23   /
      DATA BK0CS( 14) / +.2075275266 9066680831 9999999999 999 D-26   /
      DATA BK0CS( 15) / +.2712814218 0729856000 0000000000 000 D-29   /
      DATA BK0CS( 16) / +.3082593887 9146666666 6666666666 666 D-32   /
      DATA AK0CS(  1) / -.7643947903 3279414240 8297827008 8 D-1      /
      DATA AK0CS(  2) / -.2235652605 6998190520 2309555079 1 D-1      /
      DATA AK0CS(  3) / +.7734181154 6938582353 0061817404 7 D-3      /
      DATA AK0CS(  4) / -.4281006688 8860994644 5214643541 6 D-4      /
      DATA AK0CS(  5) / +.3081700173 8629747436 5001482666 0 D-5      /
      DATA AK0CS(  6) / -.2639367222 0096649740 6744889272 3 D-6      /
      DATA AK0CS(  7) / +.2563713036 4034692062 9408826574 2 D-7      /
      DATA AK0CS(  8) / -.2742705549 9002012638 5721191524 4 D-8      /
      DATA AK0CS(  9) / +.3169429658 0974995920 8083287340 3 D-9      /
      DATA AK0CS( 10) / -.3902353286 9621841416 0106571796 2 D-10     /
      DATA AK0CS( 11) / +.5068040698 1885754020 5009212728 6 D-11     /
      DATA AK0CS( 12) / -.6889574741 0078706795 4171355798 4 D-12     /
      DATA AK0CS( 13) / +.9744978497 8259176913 8820133683 1 D-13     /
      DATA AK0CS( 14) / -.1427332841 8845485053 8985534012 2 D-13     /
      DATA AK0CS( 15) / +.2156412571 0214630395 5806297652 7 D-14     /
      DATA AK0CS( 16) / -.3349654255 1495627721 8878205853 0 D-15     /
      DATA AK0CS( 17) / +.5335260216 9529116921 4528039260 1 D-16     /
      DATA AK0CS( 18) / -.8693669980 8907538076 3962237883 7 D-17     /
      DATA AK0CS( 19) / +.1446404347 8622122278 8776344234 6 D-17     /
      DATA AK0CS( 20) / -.2452889825 5001296824 0467875157 3 D-18     /
      DATA AK0CS( 21) / +.4233754526 2321715728 2170634240 0 D-19     /
      DATA AK0CS( 22) / -.7427946526 4544641956 9534129493 3 D-20     /
      DATA AK0CS( 23) / +.1323150529 3926668662 7796746240 0 D-20     /
      DATA AK0CS( 24) / -.2390587164 7396494513 3598146559 9 D-21     /
      DATA AK0CS( 25) / +.4376827585 9232261401 6571255466 6 D-22     /
      DATA AK0CS( 26) / -.8113700607 3451180593 3901141333 3 D-23     /
      DATA AK0CS( 27) / +.1521819913 8321729583 1037815466 6 D-23     /
      DATA AK0CS( 28) / -.2886041941 4833977702 3595861333 3 D-24     /
      DATA AK0CS( 29) / +.5530620667 0547179799 9261013333 3 D-25     /
      DATA AK0CS( 30) / -.1070377329 2498987285 9163306666 6 D-25     /
      DATA AK0CS( 31) / +.2091086893 1423843002 9632853333 3 D-26     /
      DATA AK0CS( 32) / -.4121713723 6462038274 1026133333 3 D-27     /
      DATA AK0CS( 33) / +.8193483971 1213076401 3568000000 0 D-28     /
      DATA AK0CS( 34) / -.1642000275 4592977267 8075733333 3 D-28     /
      DATA AK0CS( 35) / +.3316143281 4802271958 9034666666 6 D-29     /
      DATA AK0CS( 36) / -.6746863644 1452959410 8586666666 6 D-30     /
      DATA AK0CS( 37) / +.1382429146 3184246776 3541333333 3 D-30     /
      DATA AK0CS( 38) / -.2851874167 3598325708 1173333333 3 D-31     /
      DATA AK02CS(  1) / -.1201869826 3075922398 3934621245 2 D-1      /
      DATA AK02CS(  2) / -.9174852691 0256953106 5256107571 3 D-2      /
      DATA AK02CS(  3) / +.1444550931 7750058210 4884387805 7 D-3      /
      DATA AK02CS(  4) / -.4013614175 4357097286 7102107787 9 D-5      /
      DATA AK02CS(  5) / +.1567831810 8523106725 9034899033 3 D-6      /
      DATA AK02CS(  6) / -.7770110438 5217377103 1579975446 0 D-8      /
      DATA AK02CS(  7) / +.4611182576 1797178825 3313052958 6 D-9      /
      DATA AK02CS(  8) / -.3158592997 8605657705 2666580330 9 D-10     /
      DATA AK02CS(  9) / +.2435018039 3650411278 3588781432 9 D-11     /
      DATA AK02CS( 10) / -.2074331387 3983478977 0985337350 6 D-12     /
      DATA AK02CS( 11) / +.1925787280 5899170847 4273650469 3 D-13     /
      DATA AK02CS( 12) / -.1927554805 8389561036 0034718221 8 D-14     /
      DATA AK02CS( 13) / +.2062198029 1978182782 8523786964 4 D-15     /
      DATA AK02CS( 14) / -.2341685117 5792424026 0364019507 1 D-16     /
      DATA AK02CS( 15) / +.2805902810 6430422468 1517882845 8 D-17     /
      DATA AK02CS( 16) / -.3530507631 1618079458 1548246357 3 D-18     /
      DATA AK02CS( 17) / +.4645295422 9351082674 2421633706 6 D-19     /
      DATA AK02CS( 18) / -.6368625941 3442664739 2205346133 3 D-20     /
      DATA AK02CS( 19) / +.9069521310 9865155676 2234880000 0 D-21     /
      DATA AK02CS( 20) / -.1337974785 4236907398 4500531199 9 D-21     /
      DATA AK02CS( 21) / +.2039836021 8599523155 2208896000 0 D-22     /
      DATA AK02CS( 22) / -.3207027481 3678405000 6086997333 3 D-23     /
      DATA AK02CS( 23) / +.5189744413 6623099636 2635946666 6 D-24     /
      DATA AK02CS( 24) / -.8629501497 5405721929 6460799999 9 D-25     /
      DATA AK02CS( 25) / +.1472161183 1025598552 0803840000 0 D-25     /
      DATA AK02CS( 26) / -.2573069023 8670112838 1235199999 9 D-26     /
      DATA AK02CS( 27) / +.4601774086 6435165873 7664000000 0 D-27     /
      DATA AK02CS( 28) / -.8411555324 2010937371 3066666666 6 D-28     /
      DATA AK02CS( 29) / +.1569806306 6353689393 0154666666 6 D-28     /
      DATA AK02CS( 30) / -.2988226453 0057577889 7919999999 9 D-29     /
      DATA AK02CS( 31) / +.5796831375 2168365206 1866666666 6 D-30     /
      DATA AK02CS( 32) / -.1145035994 3476813321 5573333333 3 D-30     /
      DATA AK02CS( 33) / +.2301266594 2496828020 0533333333 3 D-31     /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBSK0E
      IF (FIRST) THEN
         ETA = 0.1*REAL(D1MACH(3))
         NTK0 = INITDS (BK0CS, 16, ETA)
         NTAK0 = INITDS (AK0CS, 38, ETA)
         NTAK02 = INITDS (AK02CS, 33, ETA)
         XSML = SQRT(4.0D0*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0.D0) CALL XERMSG ('SLATEC', 'DBSK0E',
     +   'X IS ZERO OR NEGATIVE', 2, 2)
      IF (X.GT.2.0D0) GO TO 20
C
      Y = 0.D0
      IF (X.GT.XSML) Y = X*X
      DBSK0E = EXP(X)*(-LOG(0.5D0*X)*DBESI0(X) - 0.25D0 +
     1  DCSEVL (.5D0*Y-1.D0, BK0CS, NTK0))
      RETURN
C
 20   IF (X.LE.8.D0) DBSK0E = (1.25D0 + DCSEVL ((16.D0/X-5.D0)/3.D0,
     1  AK0CS, NTAK0))/SQRT(X)
      IF (X.GT.8.D0) DBSK0E = (1.25D0 +
     1  DCSEVL (16.D0/X-1.D0, AK02CS, NTAK02))/SQRT(X)
C
      RETURN
      END

*DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
C***BEGIN PROLOGUE  XERMSG
C***PURPOSE  Process error messages for SLATEC and other libraries.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERMSG-A)
C***KEYWORDS  ERROR MESSAGE, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C   XERMSG processes a diagnostic message in a manner determined by the
C   value of LEVEL and the current value of the library error control
C   flag, KONTRL.  See subroutine XSETF for details.
C
C    LIBRAR   A character constant (or character variable) with the name
C             of the library.  This will be 'SLATEC' for the SLATEC
C             Common Math Library.  The error handling package is
C             general enough to be used by many libraries
C             simultaneously, so it is desirable for the routine that
C             detects and reports an error to identify the library name
C             as well as the routine name.
C
C    SUBROU   A character constant (or character variable) with the name
C             of the routine that detected the error.  Usually it is the
C             name of the routine that is calling XERMSG.  There are
C             some instances where a user callable library routine calls
C             lower level subsidiary routines where the error is
C             detected.  In such cases it may be more informative to
C             supply the name of the routine the user called rather than
C             the name of the subsidiary routine that detected the
C             error.
C
C    MESSG    A character constant (or character variable) with the text
C             of the error or warning message.  In the example below,
C             the message is a character constant that contains a
C             generic message.
C
C                   CALL XERMSG ('SLATEC', 'MMPY',
C                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
C                  *3, 1)
C
C             It is possible (and is sometimes desirable) to generate a
C             specific message--e.g., one that contains actual numeric
C             values.  Specific numeric values can be converted into
C             character strings using formatted WRITE statements into
C             character variables.  This is called standard Fortran
C             internal file I/O and is exemplified in the first three
C             lines of the following example.  You can also catenate
C             substrings of characters to construct the error message.
C             Here is an example showing the use of both writing to
C             an internal file and catenating character strings.
C
C                   CHARACTER*5 CHARN, CHARL
C                   WRITE (CHARN,10) N
C                   WRITE (CHARL,10) LDA
C                10 FORMAT(I5)
C                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
C                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
C                  *   CHARL, 3, 1)
C
C             There are two subtleties worth mentioning.  One is that
C             the // for character catenation is used to construct the
C             error message so that no single character constant is
C             continued to the next line.  This avoids confusion as to
C             whether there are trailing blanks at the end of the line.
C             The second is that by catenating the parts of the message
C             as an actual argument rather than encoding the entire
C             message into one large character variable, we avoid
C             having to know how long the message will be in order to
C             declare an adequate length for that large character
C             variable.  XERMSG calls XERPRN to print the message using
C             multiple lines if necessary.  If the message is very long,
C             XERPRN will break it into pieces of 72 characters (as
C             requested by XERMSG) for printing on multiple lines.
C             Also, XERMSG asks XERPRN to prefix each line with ' *  '
C             so that the total line length could be 76 characters.
C             Note also that XERPRN scans the error message backwards
C             to ignore trailing blanks.  Another feature is that
C             the substring '$$' is treated as a new line sentinel
C             by XERPRN.  If you want to construct a multiline
C             message without having to count out multiples of 72
C             characters, just use '$$' as a separator.  '$$'
C             obviously must occur within 72 characters of the
C             start of each line to have its intended effect since
C             XERPRN is asked to wrap around at 72 characters in
C             addition to looking for '$$'.
C
C    NERR     An integer value that is chosen by the library routine's
C             author.  It must be in the range -99 to 999 (three
C             printable digits).  Each distinct error should have its
C             own error number.  These error numbers should be described
C             in the machine readable documentation for the routine.
C             The error numbers need be unique only within each routine,
C             so it is reasonable for each routine to start enumerating
C             errors from 1 and proceeding to the next integer.
C
C    LEVEL    An integer value in the range 0 to 2 that indicates the
C             level (severity) of the error.  Their meanings are
C
C            -1  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.  An attempt is made to only print this
C                message once.
C
C             0  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.
C
C             1  A recoverable error.  This is used even if the error is
C                so serious that the routine cannot return any useful
C                answer.  If the user has told the error package to
C                return after recoverable errors, then XERMSG will
C                return to the Library routine which can then return to
C                the user's routine.  The user may also permit the error
C                package to terminate the program upon encountering a
C                recoverable error.
C
C             2  A fatal error.  XERMSG will not return to its caller
C                after it receives a fatal error.  This level should
C                hardly ever be used; it is much better to allow the
C                user a chance to recover.  An example of one of the few
C                cases in which it is permissible to declare a level 2
C                error is a reverse communication Library routine that
C                is likely to be called repeatedly until it integrates
C                across some interval.  If there is a serious error in
C                the input such that another step cannot be taken and
C                the Library routine is called again without the input
C                error having been corrected by the caller, the Library
C                routine will probably be called forever with improper
C                input.  In this case, it is reasonable to declare the
C                error to be fatal.
C
C    Each of the arguments to XERMSG is input; none will be modified by
C    XERMSG.  A routine may make multiple calls to XERMSG with warning
C    level messages; however, after a call to XERMSG with a recoverable
C    error, the routine should return to the user.  Do not try to call
C    XERMSG with a second recoverable error after the first recoverable
C    error because the error package saves the error number.  The user
C    can retrieve this error number by calling another entry point in
C    the error handling package and then clear the error number when
C    recovering from the error.  Calling XERMSG in succession causes the
C    old error number to be overwritten by the latest error number.
C    This is considered harmless for error numbers associated with
C    warning messages but must not be done for error numbers of serious
C    errors.  After a call to XERMSG with a recoverable error, the user
C    must be given a chance to call NUMXER or XERCLR to retrieve or
C    clear the error number.
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
C***REVISION HISTORY  (YYMMDD)
C   880101  DATE WRITTEN
C   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
C           THERE ARE TWO BASIC CHANGES.
C           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
C               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
C               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
C               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
C               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
C               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
C               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
C               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
C           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
C               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
C               OF LOWER CASE.
C   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
C           THE PRINCIPAL CHANGES ARE
C           1.  CLARIFY COMMENTS IN THE PROLOGUES
C           2.  RENAME XRPRNT TO XERPRN
C           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
C               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
C               CHARACTER FOR NEW RECORDS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           CLEAN UP THE CODING.
C   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
C           PREFIX.
C   891013  REVISED TO CORRECT COMMENTS.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
C           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
C           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
C           XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERMSG
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
C***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
C
C       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
C       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
C          SHOULD BE PRINTED.
C
C       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
C          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
C          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
C
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.
     *   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //
     *      'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//
     *      'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
C
C       RECORD THE MESSAGE.
C
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
C
C       HANDLE PRINT-ONCE WARNING MESSAGES.
C
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
C
C       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
C
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
C
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
C
C       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
C       ZERO AND THE ERROR IS NOT FATAL.
C
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
C
C       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
C       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
C       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
C       IS NOT ZERO.
C
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
C       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
C       FROM EACH OF THE FOLLOWING THREE OPTIONS.
C       1.  LEVEL OF THE MESSAGE
C              'INFORMATIVE MESSAGE'
C              'POTENTIALLY RECOVERABLE ERROR'
C              'FATAL ERROR'
C       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
C              'PROG CONTINUES'
C              'PROG ABORTED'
C       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
C           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
C           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
C              'TRACEBACK REQUESTED'
C              'TRACEBACK NOT REQUESTED'
C       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
C       EXCEED 74 CHARACTERS.
C       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
C
      IF (LKNTRL .GT. 0) THEN
C
C       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
C
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
C
C       THEN WHETHER THE PROGRAM WILL CONTINUE.
C
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.
     *       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
C
C       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
C
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       NOW SEND OUT THE MESSAGE.
C
      CALL XERPRN (' *  ', -1, MESSG, 72)
C
C       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
C          TRACEBACK.
C
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
C
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
C
C       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
C
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
C
C       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
C       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
C
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
C
C       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
C       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
C       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
C
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN
     *         (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
      END

*DECK DBSI0E
      DOUBLE PRECISION FUNCTION DBSI0E (X)
C***BEGIN PROLOGUE  DBSI0E
C***PURPOSE  Compute the exponentially scaled modified (hyperbolic)
C            Bessel function of the first kind of order zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (BESI0E-S, DBSI0E-D)
C***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
C             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
C             ORDER ZERO, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBSI0E(X) calculates the double precision exponentially scaled
C modified (hyperbolic) Bessel function of the first kind of order
C zero for double precision argument X.  The result is the Bessel
C function I0(X) multiplied by EXP(-ABS(X)).
C
C Series for BI0        on the interval  0.          to  9.00000E+00
C                                        with weighted error   9.51E-34
C                                         log weighted error  33.02
C                               significant figures required  33.31
C                                    decimal places required  33.65
C
C Series for AI0        on the interval  1.25000E-01 to  3.33333E-01
C                                        with weighted error   2.74E-32
C                                         log weighted error  31.56
C                               significant figures required  30.15
C                                    decimal places required  32.39
C
C Series for AI02       on the interval  0.          to  1.25000E-01
C                                        with weighted error   1.97E-32
C                                         log weighted error  31.71
C                               significant figures required  30.15
C                                    decimal places required  32.63
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, INITDS
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DBSI0E
      DOUBLE PRECISION X, BI0CS(18), AI0CS(46), AI02CS(69),
     1  XSML, Y, D1MACH, DCSEVL
      LOGICAL FIRST
      SAVE BI0CS, AI0CS, AI02CS, NTI0, NTAI0, NTAI02, XSML, FIRST
      DATA BI0CS(  1) / -.7660547252 8391449510 8189497624 3285 D-1   /
      DATA BI0CS(  2) / +.1927337953 9938082699 5240875088 1196 D+1   /
      DATA BI0CS(  3) / +.2282644586 9203013389 3702929233 0415 D+0   /
      DATA BI0CS(  4) / +.1304891466 7072904280 7933421069 1888 D-1   /
      DATA BI0CS(  5) / +.4344270900 8164874513 7868268102 6107 D-3   /
      DATA BI0CS(  6) / +.9422657686 0019346639 2317174411 8766 D-5   /
      DATA BI0CS(  7) / +.1434006289 5106910799 6209187817 9957 D-6   /
      DATA BI0CS(  8) / +.1613849069 6617490699 1541971999 4611 D-8   /
      DATA BI0CS(  9) / +.1396650044 5356696994 9509270814 2522 D-10  /
      DATA BI0CS( 10) / +.9579451725 5054453446 2752317189 3333 D-13  /
      DATA BI0CS( 11) / +.5333981859 8625021310 1510774400 0000 D-15  /
      DATA BI0CS( 12) / +.2458716088 4374707746 9678591999 9999 D-17  /
      DATA BI0CS( 13) / +.9535680890 2487700269 4434133333 3333 D-20  /
      DATA BI0CS( 14) / +.3154382039 7214273367 8933333333 3333 D-22  /
      DATA BI0CS( 15) / +.9004564101 0946374314 6666666666 6666 D-25  /
      DATA BI0CS( 16) / +.2240647369 1236700160 0000000000 0000 D-27  /
      DATA BI0CS( 17) / +.4903034603 2428373333 3333333333 3333 D-30  /
      DATA BI0CS( 18) / +.9508172606 1226666666 6666666666 6666 D-33  /
      DATA AI0CS(  1) / +.7575994494 0237959427 2987203743 8 D-1      /
      DATA AI0CS(  2) / +.7591380810 8233455072 9297873320 4 D-2      /
      DATA AI0CS(  3) / +.4153131338 9237505018 6319749138 2 D-3      /
      DATA AI0CS(  4) / +.1070076463 4390730735 8242970217 0 D-4      /
      DATA AI0CS(  5) / -.7901179979 2128946607 5031948573 0 D-5      /
      DATA AI0CS(  6) / -.7826143501 4387522697 8898980690 9 D-6      /
      DATA AI0CS(  7) / +.2783849942 9488708063 8118538985 7 D-6      /
      DATA AI0CS(  8) / +.8252472600 6120271919 6682913319 8 D-8      /
      DATA AI0CS(  9) / -.1204463945 5201991790 5496089110 3 D-7      /
      DATA AI0CS( 10) / +.1559648598 5060764436 1228752792 8 D-8      /
      DATA AI0CS( 11) / +.2292556367 1033165434 7725480285 7 D-9      /
      DATA AI0CS( 12) / -.1191622884 2790646036 7777423447 8 D-9      /
      DATA AI0CS( 13) / +.1757854916 0324098302 1833124774 3 D-10     /
      DATA AI0CS( 14) / +.1128224463 2189005171 4441135682 4 D-11     /
      DATA AI0CS( 15) / -.1146848625 9272988777 2963387698 2 D-11     /
      DATA AI0CS( 16) / +.2715592054 8036628726 4365192160 6 D-12     /
      DATA AI0CS( 17) / -.2415874666 5626878384 4247572028 1 D-13     /
      DATA AI0CS( 18) / -.6084469888 2551250646 0609963922 4 D-14     /
      DATA AI0CS( 19) / +.3145705077 1754772937 0836026730 3 D-14     /
      DATA AI0CS( 20) / -.7172212924 8711877179 6217505917 6 D-15     /
      DATA AI0CS( 21) / +.7874493403 4541033960 8390960332 7 D-16     /
      DATA AI0CS( 22) / +.1004802753 0094624023 4524457183 9 D-16     /
      DATA AI0CS( 23) / -.7566895365 3505348534 2843588881 0 D-17     /
      DATA AI0CS( 24) / +.2150380106 8761198878 1205128784 5 D-17     /
      DATA AI0CS( 25) / -.3754858341 8308744291 5158445260 8 D-18     /
      DATA AI0CS( 26) / +.2354065842 2269925769 0075710532 2 D-19     /
      DATA AI0CS( 27) / +.1114667612 0479285302 2637335511 0 D-19     /
      DATA AI0CS( 28) / -.5398891884 3969903786 9677932270 9 D-20     /
      DATA AI0CS( 29) / +.1439598792 2407526770 4285840452 2 D-20     /
      DATA AI0CS( 30) / -.2591916360 1110934064 6081840196 2 D-21     /
      DATA AI0CS( 31) / +.2238133183 9985839074 3409229824 0 D-22     /
      DATA AI0CS( 32) / +.5250672575 3647711727 7221683199 9 D-23     /
      DATA AI0CS( 33) / -.3249904138 5332307841 7343228586 6 D-23     /
      DATA AI0CS( 34) / +.9924214103 2050379278 5728471040 0 D-24     /
      DATA AI0CS( 35) / -.2164992254 2446695231 4655429973 3 D-24     /
      DATA AI0CS( 36) / +.3233609471 9435940839 7333299199 9 D-25     /
      DATA AI0CS( 37) / -.1184620207 3967424898 2473386666 6 D-26     /
      DATA AI0CS( 38) / -.1281671853 9504986505 4833868799 9 D-26     /
      DATA AI0CS( 39) / +.5827015182 2793905116 0556885333 3 D-27     /
      DATA AI0CS( 40) / -.1668222326 0261097193 6450150399 9 D-27     /
      DATA AI0CS( 41) / +.3625309510 5415699757 0068480000 0 D-28     /
      DATA AI0CS( 42) / -.5733627999 0557135899 4595839999 9 D-29     /
      DATA AI0CS( 43) / +.3736796722 0630982296 4258133333 3 D-30     /
      DATA AI0CS( 44) / +.1602073983 1568519633 6551253333 3 D-30     /
      DATA AI0CS( 45) / -.8700424864 0572298845 2249599999 9 D-31     /
      DATA AI0CS( 46) / +.2741320937 9374811456 0341333333 3 D-31     /
      DATA AI02CS(  1) / +.5449041101 4108831607 8960962268 0 D-1      /
      DATA AI02CS(  2) / +.3369116478 2556940898 9785662979 9 D-2      /
      DATA AI02CS(  3) / +.6889758346 9168239842 6263914301 1 D-4      /
      DATA AI02CS(  4) / +.2891370520 8347564829 6692402323 2 D-5      /
      DATA AI02CS(  5) / +.2048918589 4690637418 2760534093 1 D-6      /
      DATA AI02CS(  6) / +.2266668990 4981780645 9327743136 1 D-7      /
      DATA AI02CS(  7) / +.3396232025 7083863451 5084396952 3 D-8      /
      DATA AI02CS(  8) / +.4940602388 2249695891 0482449783 5 D-9      /
      DATA AI02CS(  9) / +.1188914710 7846438342 4084525196 3 D-10     /
      DATA AI02CS( 10) / -.3149916527 9632413645 3864862961 9 D-10     /
      DATA AI02CS( 11) / -.1321581184 0447713118 7540739926 7 D-10     /
      DATA AI02CS( 12) / -.1794178531 5068061177 7943574026 9 D-11     /
      DATA AI02CS( 13) / +.7180124451 3836662336 7106429346 9 D-12     /
      DATA AI02CS( 14) / +.3852778382 7421427011 4089801777 6 D-12     /
      DATA AI02CS( 15) / +.1540086217 5214098269 1325823339 7 D-13     /
      DATA AI02CS( 16) / -.4150569347 2872220866 2689972015 6 D-13     /
      DATA AI02CS( 17) / -.9554846698 8283076487 0214494312 5 D-14     /
      DATA AI02CS( 18) / +.3811680669 3526224207 4605535511 8 D-14     /
      DATA AI02CS( 19) / +.1772560133 0565263836 0493266675 8 D-14     /
      DATA AI02CS( 20) / -.3425485619 6772191346 1924790328 2 D-15     /
      DATA AI02CS( 21) / -.2827623980 5165834849 4205593759 4 D-15     /
      DATA AI02CS( 22) / +.3461222867 6974610930 9706250813 4 D-16     /
      DATA AI02CS( 23) / +.4465621420 2967599990 1042054284 3 D-16     /
      DATA AI02CS( 24) / -.4830504485 9441820712 5525403795 4 D-17     /
      DATA AI02CS( 25) / -.7233180487 8747539545 6227240924 5 D-17     /
      DATA AI02CS( 26) / +.9921475412 1736985988 8046093981 0 D-18     /
      DATA AI02CS( 27) / +.1193650890 8459820855 0439949924 2 D-17     /
      DATA AI02CS( 28) / -.2488709837 1508072357 2054491660 2 D-18     /
      DATA AI02CS( 29) / -.1938426454 1609059289 8469781132 6 D-18     /
      DATA AI02CS( 30) / +.6444656697 3734438687 8301949394 9 D-19     /
      DATA AI02CS( 31) / +.2886051596 2892243264 8171383073 4 D-19     /
      DATA AI02CS( 32) / -.1601954907 1749718070 6167156200 7 D-19     /
      DATA AI02CS( 33) / -.3270815010 5923147208 9193567485 9 D-20     /
      DATA AI02CS( 34) / +.3686932283 8264091811 4600723939 3 D-20     /
      DATA AI02CS( 35) / +.1268297648 0309501530 1359529710 9 D-22     /
      DATA AI02CS( 36) / -.7549825019 3772739076 9636664410 1 D-21     /
      DATA AI02CS( 37) / +.1502133571 3778353496 3712789053 4 D-21     /
      DATA AI02CS( 38) / +.1265195883 5096485349 3208799248 3 D-21     /
      DATA AI02CS( 39) / -.6100998370 0836807086 2940891600 2 D-22     /
      DATA AI02CS( 40) / -.1268809629 2601282643 6872095924 2 D-22     /
      DATA AI02CS( 41) / +.1661016099 8907414578 4038487490 5 D-22     /
      DATA AI02CS( 42) / -.1585194335 7658855793 7970504881 4 D-23     /
      DATA AI02CS( 43) / -.3302645405 9682178009 5381766755 6 D-23     /
      DATA AI02CS( 44) / +.1313580902 8392397817 4039623117 4 D-23     /
      DATA AI02CS( 45) / +.3689040246 6711567933 1425637280 4 D-24     /
      DATA AI02CS( 46) / -.4210141910 4616891492 1978247249 9 D-24     /
      DATA AI02CS( 47) / +.4791954591 0828657806 3171401373 0 D-25     /
      DATA AI02CS( 48) / +.8459470390 2218217952 9971707412 4 D-25     /
      DATA AI02CS( 49) / -.4039800940 8728324931 4607937181 0 D-25     /
      DATA AI02CS( 50) / -.6434714653 6504313473 0100850469 5 D-26     /
      DATA AI02CS( 51) / +.1225743398 8756659903 4464736990 5 D-25     /
      DATA AI02CS( 52) / -.2934391316 0257089231 9879821175 4 D-26     /
      DATA AI02CS( 53) / -.1961311309 1949829262 0371205728 9 D-26     /
      DATA AI02CS( 54) / +.1503520374 8221934241 6229900309 8 D-26     /
      DATA AI02CS( 55) / -.9588720515 7448265520 3386388206 9 D-28     /
      DATA AI02CS( 56) / -.3483339380 8170454863 9441108511 4 D-27     /
      DATA AI02CS( 57) / +.1690903610 2630436730 6244960725 6 D-27     /
      DATA AI02CS( 58) / +.1982866538 7356030438 9400115718 8 D-28     /
      DATA AI02CS( 59) / -.5317498081 4918162145 7583002528 4 D-28     /
      DATA AI02CS( 60) / +.1803306629 8883929462 3501450390 1 D-28     /
      DATA AI02CS( 61) / +.6213093341 4548931758 8405311242 2 D-29     /
      DATA AI02CS( 62) / -.7692189292 7721618632 0072806673 0 D-29     /
      DATA AI02CS( 63) / +.1858252826 1117025426 2556016596 3 D-29     /
      DATA AI02CS( 64) / +.1237585142 2813957248 9927154554 1 D-29     /
      DATA AI02CS( 65) / -.1102259120 4092238032 1779478779 2 D-29     /
      DATA AI02CS( 66) / +.1886287118 0397044900 7787447943 1 D-30     /
      DATA AI02CS( 67) / +.2160196872 2436589131 4903141406 0 D-30     /
      DATA AI02CS( 68) / -.1605454124 9197432005 8446594965 5 D-30     /
      DATA AI02CS( 69) / +.1965352984 5942906039 3884807331 8 D-31     /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBSI0E
      IF (FIRST) THEN
         ETA = 0.1*REAL(D1MACH(3))
         NTI0 = INITDS (BI0CS, 18, ETA)
         NTAI0 = INITDS (AI0CS, 46, ETA)
         NTAI02 = INITDS (AI02CS, 69, ETA)
         XSML = SQRT(4.5D0*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.3.0D0) GO TO 20
C
      DBSI0E = 1.0D0 - X
      IF (Y.GT.XSML) DBSI0E = EXP(-Y) * (2.75D0 +
     1  DCSEVL (Y*Y/4.5D0-1.D0, BI0CS, NTI0) )
      RETURN
C
 20   IF (Y.LE.8.D0) DBSI0E = (0.375D0 + DCSEVL ((48.D0/Y-11.D0)/5.D0,
     1  AI0CS, NTAI0))/SQRT(Y)
      IF (Y.GT.8.D0) DBSI0E = (0.375D0 + DCSEVL (16.D0/Y-1.D0, AI02CS,
     1  NTAI02))/SQRT(Y)
C
      RETURN
      END
*DECK DBESI0
      DOUBLE PRECISION FUNCTION DBESI0 (X)
C***BEGIN PROLOGUE  DBESI0
C***PURPOSE  Compute the hyperbolic Bessel function of the first kind
C            of order zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (BESI0-S, DBESI0-D)
C***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESI0(X) calculates the double precision modified (hyperbolic)
C Bessel function of the first kind of order zero and double
C precision argument X.
C
C Series for BI0        on the interval  0.          to  9.00000E+00
C                                        with weighted error   9.51E-34
C                                         log weighted error  33.02
C                               significant figures required  33.31
C                                    decimal places required  33.65
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DBSI0E, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBESI0
      DOUBLE PRECISION X, BI0CS(18), XMAX, XSML, Y, D1MACH,
     1  DCSEVL, DBSI0E
      LOGICAL FIRST
      SAVE BI0CS, NTI0, XSML, XMAX, FIRST
      DATA BI0CS(  1) / -.7660547252 8391449510 8189497624 3285 D-1   /
      DATA BI0CS(  2) / +.1927337953 9938082699 5240875088 1196 D+1   /
      DATA BI0CS(  3) / +.2282644586 9203013389 3702929233 0415 D+0   /
      DATA BI0CS(  4) / +.1304891466 7072904280 7933421069 1888 D-1   /
      DATA BI0CS(  5) / +.4344270900 8164874513 7868268102 6107 D-3   /
      DATA BI0CS(  6) / +.9422657686 0019346639 2317174411 8766 D-5   /
      DATA BI0CS(  7) / +.1434006289 5106910799 6209187817 9957 D-6   /
      DATA BI0CS(  8) / +.1613849069 6617490699 1541971999 4611 D-8   /
      DATA BI0CS(  9) / +.1396650044 5356696994 9509270814 2522 D-10  /
      DATA BI0CS( 10) / +.9579451725 5054453446 2752317189 3333 D-13  /
      DATA BI0CS( 11) / +.5333981859 8625021310 1510774400 0000 D-15  /
      DATA BI0CS( 12) / +.2458716088 4374707746 9678591999 9999 D-17  /
      DATA BI0CS( 13) / +.9535680890 2487700269 4434133333 3333 D-20  /
      DATA BI0CS( 14) / +.3154382039 7214273367 8933333333 3333 D-22  /
      DATA BI0CS( 15) / +.9004564101 0946374314 6666666666 6666 D-25  /
      DATA BI0CS( 16) / +.2240647369 1236700160 0000000000 0000 D-27  /
      DATA BI0CS( 17) / +.4903034603 2428373333 3333333333 3333 D-30  /
      DATA BI0CS( 18) / +.9508172606 1226666666 6666666666 6666 D-33  /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESI0
      IF (FIRST) THEN
         NTI0 = INITDS (BI0CS, 18, 0.1*REAL(D1MACH(3)))
         XSML = SQRT(4.5D0*D1MACH(3))
         XMAX = LOG (D1MACH(2))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.3.0D0) GO TO 20
C
      DBESI0 = 1.0D0
      IF (Y.GT.XSML) DBESI0 = 2.75D0 + DCSEVL (Y*Y/4.5D0-1.D0, BI0CS,
     1  NTI0)
      RETURN
C
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'DBESI0',
     +   'ABS(X) SO BIG I0 OVERFLOWS', 2, 2)
C
      DBESI0 = EXP(Y) * DBSI0E(X)
C
      RETURN
      END

*DECK XERPRN
      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
C***BEGIN PROLOGUE  XERPRN
C***SUBSIDIARY
C***PURPOSE  Print error messages processed by XERMSG.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERPRN-A)
C***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C This routine sends one or more lines to each of the (up to five)
C logical units to which error messages are to be sent.  This routine
C is called several times by XERMSG, sometimes with a single line to
C print and sometimes with a (potentially very long) message that may
C wrap around into multiple lines.
C
C PREFIX  Input argument of type CHARACTER.  This argument contains
C         characters to be put at the beginning of each line before
C         the body of the message.  No more than 16 characters of
C         PREFIX will be used.
C
C NPREF   Input argument of type INTEGER.  This argument is the number
C         of characters to use from PREFIX.  If it is negative, the
C         intrinsic function LEN is used to determine its length.  If
C         it is zero, PREFIX is not used.  If it exceeds 16 or if
C         LEN(PREFIX) exceeds 16, only the first 16 characters will be
C         used.  If NPREF is positive and the length of PREFIX is less
C         than NPREF, a copy of PREFIX extended with blanks to length
C         NPREF will be used.
C
C MESSG   Input argument of type CHARACTER.  This is the text of a
C         message to be printed.  If it is a long message, it will be
C         broken into pieces for printing on multiple lines.  Each line
C         will start with the appropriate prefix and be followed by a
C         piece of the message.  NWRAP is the number of characters per
C         piece; that is, after each NWRAP characters, we break and
C         start a new line.  In addition the characters '$$' embedded
C         in MESSG are a sentinel for a new line.  The counting of
C         characters up to NWRAP starts over for each new line.  The
C         value of NWRAP typically used by XERMSG is 72 since many
C         older error messages in the SLATEC Library are laid out to
C         rely on wrap-around every 72 characters.
C
C NWRAP   Input argument of type INTEGER.  This gives the maximum size
C         piece into which to break MESSG for printing on multiple
C         lines.  An embedded '$$' ends a line, and the count restarts
C         at the following character.  If a line break does not occur
C         on a blank (it would split a word) that word is moved to the
C         next line.  Values of NWRAP less than 16 will be treated as
C         16.  Values of NWRAP greater than 132 will be treated as 132.
C         The actual line length will be NPREF + NWRAP after NPREF has
C         been adjusted to fall between 0 and 16 and NWRAP has been
C         adjusted to fall between 16 and 132.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   880621  DATE WRITTEN
C   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
C           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
C           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
C           SLASH CHARACTER IN FORMAT STATEMENTS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
C           LINES TO BE PRINTED.
C   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
C           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
C   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Added code to break messages between words.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERPRN
      CHARACTER*(*) PREFIX, MESSG
      INTEGER NPREF, NWRAP
      CHARACTER*148 CBUFF
      INTEGER IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
C***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
C
C       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
C       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
C       ERROR MESSAGE UNIT.
C
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) .EQ. 0) IU(I) = N
   10 CONTINUE
C
C       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
C       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
C       THE REST OF THIS ROUTINE.
C
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
C
C       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
C       TIME FROM MESSG TO PRINT ON ONE LINE.
C
      LWRAP = MAX(16, MIN(132, NWRAP))
C
C       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
C
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
C
C       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
C
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
C
C       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
C       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
C       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
C       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
C
C       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
C       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
C       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
C       OF THE SECOND ARGUMENT.
C
C       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
C       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
C       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
C       POSITION NEXTC.
C
C       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
C                       REMAINDER OF THE CHARACTER STRING.  LPIECE
C                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
C                       WHICHEVER IS LESS.
C
C       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
C                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
C                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
C                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
C                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
C                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
C                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
C                       SHOULD BE INCREMENTED BY 2.
C
C       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
C
C       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
C                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
C                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
C                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
C                       AT THE END OF A LINE.
C
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
C
C       THERE WAS NO NEW LINE SENTINEL FOUND.
C
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
C
C       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
C       DON'T PRINT A BLANK LINE.
C
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
C
C       LPIECE SHOULD BE SET DOWN TO LWRAP.
C
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
C
C       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
C       WE SHOULD DECREMENT LPIECE BY ONE.
C
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
C
C       PRINT
C
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
C
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN
      END
*DECK XERSVE
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,
     +   ICOUNT)
C***BEGIN PROLOGUE  XERSVE
C***SUBSIDIARY
C***PURPOSE  Record that an error has occurred.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (XERSVE-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C *Usage:
C
C        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
C        CHARACTER * (len) LIBRAR, SUBROU, MESSG
C
C        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
C
C *Arguments:
C
C        LIBRAR :IN    is the library that the message is from.
C        SUBROU :IN    is the subroutine that the message is from.
C        MESSG  :IN    is the message to be saved.
C        KFLAG  :IN    indicates the action to be performed.
C                      when KFLAG > 0, the message in MESSG is saved.
C                      when KFLAG=0 the tables will be dumped and
C                      cleared.
C                      when KFLAG < 0, the tables will be dumped and
C                      not cleared.
C        NERR   :IN    is the error number.
C        LEVEL  :IN    is the error severity.
C        ICOUNT :OUT   the number of times this message has been seen,
C                      or zero if the table has overflowed and does not
C                      contain this message specifically.  When KFLAG=0,
C                      ICOUNT will not be altered.
C
C *Description:
C
C   Record that this error occurred and possibly dump and clear the
C   tables.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   800319  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900413  Routine modified to remove reference to KFLAG.  (WRB)
C   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
C           sequence, use IF-THEN-ELSE, make number of saved entries
C           easily changeable, changed routine name from XERSAV to
C           XERSVE.  (RWC)
C   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERSVE
      PARAMETER (LENTAB=10)
      INTEGER LUN(5)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
C***FIRST EXECUTABLE STATEMENT  XERSVE
C
      IF (KFLAG.LE.0) THEN
C
C        Dump the table.
C
         IF (NMSG.EQ.0) RETURN
C
C        Print to each unit.
C
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C
C           Print the table header.
C
            WRITE (IUNIT,9000)
C
C           Print body of table.
C
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),
     *            NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
C
C           Print number of other errors.
C
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
C
C        Clear the error tables.
C
         IF (KFLAG.EQ.0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
C
C        PROCESS A MESSAGE...
C        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
C
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.
     *         MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.
     *         LEVEL.EQ.LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
C
         IF (NMSG.LT.LENTAB) THEN
C
C           Empty slot found for new message.
C
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
C
C           Table is full.
C
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
C
C     Formats.
C
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /
     +   ' LIBRARY    SUBROUTINE MESSAGE START             NERR',
     +   '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END




*DECK XERHLT
      SUBROUTINE XERHLT (MESSG)
C***BEGIN PROLOGUE  XERHLT
C***SUBSIDIARY
C***PURPOSE  Abort program execution and print error message.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERHLT-A)
C***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        ***Note*** machine dependent routine
C        XERHLT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG is as in XERMSG.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to delete length of character
C           and changed routine name from XERABT to XERHLT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERHLT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERHLT
      STOP
      END

*DECK XERCNT
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
C***BEGIN PROLOGUE  XERCNT
C***SUBSIDIARY
C***PURPOSE  Allow user control over handling of errors.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERCNT-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCNT.
C        If the user has provided his own version of XERCNT, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        LIBRAR - the library that the routine is in.
C        SUBROU - the subroutine that XERMSG is being called from
C        MESSG  - the first 20 characters of the error message.
C        NERR   - same as in the call to XERMSG.
C        LEVEL  - same as in the call to XERMSG.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
C           names, changed routine name from XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERCNT
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
C***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
      END

*DECK DBESK1
      DOUBLE PRECISION FUNCTION DBESK1 (X)
C***BEGIN PROLOGUE  DBESK1
C***PURPOSE  Compute the modified (hyperbolic) Bessel function of the
C            third kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (BESK1-S, DBESK1-D)
C***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESK1(X) calculates the double precision modified (hyperbolic)
C Bessel function of the third kind of order one for double precision
C argument X.  The argument must be large enough that the result does
C not overflow and small enough that the result does not underflow.
C
C Series for BK1        on the interval  0.          to  4.00000E+00
C                                        with weighted error   9.16E-32
C                                         log weighted error  31.04
C                               significant figures required  30.61
C                                    decimal places required  31.64
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DBESI1, DBSK1E, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBESK1
      DOUBLE PRECISION X, BK1CS(16), XMAX, XMAXT, XMIN, XSML, Y,
     1  D1MACH, DCSEVL, DBESI1, DBSK1E
      LOGICAL FIRST
      SAVE BK1CS, NTK1, XMIN, XSML, XMAX, FIRST
      DATA BK1CS(  1) / +.2530022733 8947770532 5311208685 33 D-1     /
      DATA BK1CS(  2) / -.3531559607 7654487566 7238316918 01 D+0     /
      DATA BK1CS(  3) / -.1226111808 2265714823 4790679300 42 D+0     /
      DATA BK1CS(  4) / -.6975723859 6398643501 8129202960 83 D-2     /
      DATA BK1CS(  5) / -.1730288957 5130520630 1765073689 79 D-3     /
      DATA BK1CS(  6) / -.2433406141 5659682349 6007350301 64 D-5     /
      DATA BK1CS(  7) / -.2213387630 7347258558 3152525451 26 D-7     /
      DATA BK1CS(  8) / -.1411488392 6335277610 9583302126 08 D-9     /
      DATA BK1CS(  9) / -.6666901694 1993290060 8537512643 73 D-12    /
      DATA BK1CS( 10) / -.2427449850 5193659339 2631968648 53 D-14    /
      DATA BK1CS( 11) / -.7023863479 3862875971 7837971200 00 D-17    /
      DATA BK1CS( 12) / -.1654327515 5100994675 4910293333 33 D-19    /
      DATA BK1CS( 13) / -.3233834745 9944491991 8933333333 33 D-22    /
      DATA BK1CS( 14) / -.5331275052 9265274999 4666666666 66 D-25    /
      DATA BK1CS( 15) / -.7513040716 2157226666 6666666666 66 D-28    /
      DATA BK1CS( 16) / -.9155085717 6541866666 6666666666 66 D-31    /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESK1
      IF (FIRST) THEN
         NTK1 = INITDS (BK1CS, 16, 0.1*REAL(D1MACH(3)))
         XMIN = EXP(MAX(LOG(D1MACH(1)), -LOG(D1MACH(2))) + 0.01D0)
         XSML = SQRT(4.0D0*D1MACH(3))
         XMAXT = -LOG(D1MACH(1))
         XMAX = XMAXT - 0.5D0*XMAXT*LOG(XMAXT)/(XMAXT+0.5D0)
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0.D0) CALL XERMSG ('SLATEC', 'DBESK1',
     +   'X IS ZERO OR NEGATIVE', 2, 2)
      IF (X.GT.2.0D0) GO TO 20
C
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'DBESK1',
     +   'X SO SMALL K1 OVERFLOWS', 3, 2)
      Y = 0.D0
      IF (X.GT.XSML) Y = X*X
      DBESK1 = LOG(0.5D0*X)*DBESI1(X) + (0.75D0 + DCSEVL (.5D0*Y-1.D0,
     1  BK1CS, NTK1))/X
      RETURN
C
 20   DBESK1 = 0.D0
      IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'DBESK1',
     +   'X SO BIG K1 UNDERFLOWS', 1, 1)
      IF (X.GT.XMAX) RETURN
C
      DBESK1 = EXP(-X) * DBSK1E(X)
C
      RETURN
      END
*DECK DBSK1E
      DOUBLE PRECISION FUNCTION DBSK1E (X)
C***BEGIN PROLOGUE  DBSK1E
C***PURPOSE  Compute the exponentially scaled modified (hyperbolic)
C            Bessel function of the third kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (BESK1E-S, DBSK1E-D)
C***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBSK1E(S) computes the double precision exponentially scaled
C modified (hyperbolic) Bessel function of the third kind of order
C one for positive double precision argument X.
C
C Series for BK1        on the interval  0.          to  4.00000E+00
C                                        with weighted error   9.16E-32
C                                         log weighted error  31.04
C                               significant figures required  30.61
C                                    decimal places required  31.64
C
C Series for AK1        on the interval  1.25000E-01 to  5.00000E-01
C                                        with weighted error   3.07E-32
C                                         log weighted error  31.51
C                               significant figures required  30.71
C                                    decimal places required  32.30
C
C Series for AK12       on the interval  0.          to  1.25000E-01
C                                        with weighted error   2.41E-32
C                                         log weighted error  31.62
C                               significant figures required  30.25
C                                    decimal places required  32.38
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DBESI1, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBSK1E
      DOUBLE PRECISION X, BK1CS(16), AK1CS(38), AK12CS(33), XMIN,
     1  XSML, Y, D1MACH, DCSEVL, DBESI1
      LOGICAL FIRST
      SAVE BK1CS, AK1CS, AK12CS, NTK1, NTAK1, NTAK12, XMIN, XSML,
     1  FIRST
      DATA BK1CS(  1) / +.2530022733 8947770532 5311208685 33 D-1     /
      DATA BK1CS(  2) / -.3531559607 7654487566 7238316918 01 D+0     /
      DATA BK1CS(  3) / -.1226111808 2265714823 4790679300 42 D+0     /
      DATA BK1CS(  4) / -.6975723859 6398643501 8129202960 83 D-2     /
      DATA BK1CS(  5) / -.1730288957 5130520630 1765073689 79 D-3     /
      DATA BK1CS(  6) / -.2433406141 5659682349 6007350301 64 D-5     /
      DATA BK1CS(  7) / -.2213387630 7347258558 3152525451 26 D-7     /
      DATA BK1CS(  8) / -.1411488392 6335277610 9583302126 08 D-9     /
      DATA BK1CS(  9) / -.6666901694 1993290060 8537512643 73 D-12    /
      DATA BK1CS( 10) / -.2427449850 5193659339 2631968648 53 D-14    /
      DATA BK1CS( 11) / -.7023863479 3862875971 7837971200 00 D-17    /
      DATA BK1CS( 12) / -.1654327515 5100994675 4910293333 33 D-19    /
      DATA BK1CS( 13) / -.3233834745 9944491991 8933333333 33 D-22    /
      DATA BK1CS( 14) / -.5331275052 9265274999 4666666666 66 D-25    /
      DATA BK1CS( 15) / -.7513040716 2157226666 6666666666 66 D-28    /
      DATA BK1CS( 16) / -.9155085717 6541866666 6666666666 66 D-31    /
      DATA AK1CS(  1) / +.2744313406 9738829695 2576662272 66 D+0     /
      DATA AK1CS(  2) / +.7571989953 1993678170 8923781492 90 D-1     /
      DATA AK1CS(  3) / -.1441051556 4754061229 8531161756 25 D-2     /
      DATA AK1CS(  4) / +.6650116955 1257479394 2513854770 36 D-4     /
      DATA AK1CS(  5) / -.4369984709 5201407660 5808450891 67 D-5     /
      DATA AK1CS(  6) / +.3540277499 7630526799 4171390085 34 D-6     /
      DATA AK1CS(  7) / -.3311163779 2932920208 9826882457 04 D-7     /
      DATA AK1CS(  8) / +.3445977581 9010534532 3114997709 92 D-8     /
      DATA AK1CS(  9) / -.3898932347 4754271048 9819374927 58 D-9     /
      DATA AK1CS( 10) / +.4720819750 4658356400 9474493390 05 D-10    /
      DATA AK1CS( 11) / -.6047835662 8753562345 3735915628 90 D-11    /
      DATA AK1CS( 12) / +.8128494874 8658747888 1938379856 63 D-12    /
      DATA AK1CS( 13) / -.1138694574 7147891428 9239159510 42 D-12    /
      DATA AK1CS( 14) / +.1654035840 8462282325 9729482050 90 D-13    /
      DATA AK1CS( 15) / -.2480902567 7068848221 5160104405 33 D-14    /
      DATA AK1CS( 16) / +.3829237890 7024096948 4292272991 57 D-15    /
      DATA AK1CS( 17) / -.6064734104 0012418187 7682103773 86 D-16    /
      DATA AK1CS( 18) / +.9832425623 2648616038 1940046506 66 D-17    /
      DATA AK1CS( 19) / -.1628416873 8284380035 6666201156 26 D-17    /
      DATA AK1CS( 20) / +.2750153649 6752623718 2841203370 66 D-18    /
      DATA AK1CS( 21) / -.4728966646 3953250924 2810695680 00 D-19    /
      DATA AK1CS( 22) / +.8268150002 8109932722 3920503466 66 D-20    /
      DATA AK1CS( 23) / -.1468140513 6624956337 1939648853 33 D-20    /
      DATA AK1CS( 24) / +.2644763926 9208245978 0858948266 66 D-21    /
      DATA AK1CS( 25) / -.4829015756 4856387897 9698688000 00 D-22    /
      DATA AK1CS( 26) / +.8929302074 3610130180 6563327999 99 D-23    /
      DATA AK1CS( 27) / -.1670839716 8972517176 9977514666 66 D-23    /
      DATA AK1CS( 28) / +.3161645603 4040694931 3686186666 66 D-24    /
      DATA AK1CS( 29) / -.6046205531 2274989106 5064106666 66 D-25    /
      DATA AK1CS( 30) / +.1167879894 2042732700 7184213333 33 D-25    /
      DATA AK1CS( 31) / -.2277374158 2653996232 8678400000 00 D-26    /
      DATA AK1CS( 32) / +.4481109730 0773675795 3058133333 33 D-27    /
      DATA AK1CS( 33) / -.8893288476 9020194062 3360000000 00 D-28    /
      DATA AK1CS( 34) / +.1779468001 8850275131 3920000000 00 D-28    /
      DATA AK1CS( 35) / -.3588455596 7329095821 9946666666 66 D-29    /
      DATA AK1CS( 36) / +.7290629049 2694257991 6799999999 99 D-30    /
      DATA AK1CS( 37) / -.1491844984 5546227073 0240000000 00 D-30    /
      DATA AK1CS( 38) / +.3073657387 2934276300 7999999999 99 D-31    /
      DATA AK12CS(  1) / +.6379308343 7390010366 0048853410 2 D-1      /
      DATA AK12CS(  2) / +.2832887813 0497209358 3503028470 8 D-1      /
      DATA AK12CS(  3) / -.2475370673 9052503454 1454556673 2 D-3      /
      DATA AK12CS(  4) / +.5771972451 6072488204 7097662576 3 D-5      /
      DATA AK12CS(  5) / -.2068939219 5365483027 4553319655 2 D-6      /
      DATA AK12CS(  6) / +.9739983441 3818041803 0921309788 7 D-8      /
      DATA AK12CS(  7) / -.5585336140 3806249846 8889551112 9 D-9      /
      DATA AK12CS(  8) / +.3732996634 0461852402 2121285473 1 D-10     /
      DATA AK12CS(  9) / -.2825051961 0232254451 3506575492 8 D-11     /
      DATA AK12CS( 10) / +.2372019002 4841441736 4349695548 6 D-12     /
      DATA AK12CS( 11) / -.2176677387 9917539792 6830166793 8 D-13     /
      DATA AK12CS( 12) / +.2157914161 6160324539 3956268970 6 D-14     /
      DATA AK12CS( 13) / -.2290196930 7182692759 9155133815 4 D-15     /
      DATA AK12CS( 14) / +.2582885729 8232749619 1993956522 6 D-16     /
      DATA AK12CS( 15) / -.3076752641 2684631876 2109817344 0 D-17     /
      DATA AK12CS( 16) / +.3851487721 2804915970 9489684479 9 D-18     /
      DATA AK12CS( 17) / -.5044794897 6415289771 1728250880 0 D-19     /
      DATA AK12CS( 18) / +.6888673850 4185442370 1829222399 9 D-20     /
      DATA AK12CS( 19) / -.9775041541 9501183030 0213248000 0 D-21     /
      DATA AK12CS( 20) / +.1437416218 5238364610 0165973333 3 D-21     /
      DATA AK12CS( 21) / -.2185059497 3443473734 9973333333 3 D-22     /
      DATA AK12CS( 22) / +.3426245621 8092206316 4538880000 0 D-23     /
      DATA AK12CS( 23) / -.5531064394 2464082325 0124800000 0 D-24     /
      DATA AK12CS( 24) / +.9176601505 6859954037 8282666666 6 D-25     /
      DATA AK12CS( 25) / -.1562287203 6180249114 4874666666 6 D-25     /
      DATA AK12CS( 26) / +.2725419375 4843331323 4943999999 9 D-26     /
      DATA AK12CS( 27) / -.4865674910 0748279923 7802666666 6 D-27     /
      DATA AK12CS( 28) / +.8879388552 7235025873 5786666666 6 D-28     /
      DATA AK12CS( 29) / -.1654585918 0392575489 3653333333 3 D-28     /
      DATA AK12CS( 30) / +.3145111321 3578486743 0399999999 9 D-29     /
      DATA AK12CS( 31) / -.6092998312 1931276124 1600000000 0 D-30     /
      DATA AK12CS( 32) / +.1202021939 3698158346 2399999999 9 D-30     /
      DATA AK12CS( 33) / -.2412930801 4594088413 8666666666 6 D-31     /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBSK1E
      IF (FIRST) THEN
         ETA = 0.1*REAL(D1MACH(3))
         NTK1 = INITDS (BK1CS, 16, ETA)
         NTAK1 = INITDS (AK1CS, 38, ETA)
         NTAK12 = INITDS (AK12CS, 33, ETA)
C
         XMIN = EXP (MAX(LOG(D1MACH(1)), -LOG(D1MACH(2))) + 0.01D0)
         XSML = SQRT(4.0D0*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0.D0) CALL XERMSG ('SLATEC', 'DBSK1E',
     +   'X IS ZERO OR NEGATIVE', 2, 2)
      IF (X.GT.2.0D0) GO TO 20
C
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'DBSK1E',
     +   'X SO SMALL K1 OVERFLOWS', 3, 2)
      Y = 0.D0
      IF (X.GT.XSML) Y = X*X
      DBSK1E = EXP(X)*(LOG(0.5D0*X)*DBESI1(X) + (0.75D0 +
     1  DCSEVL (0.5D0*Y-1.D0, BK1CS, NTK1))/X )
      RETURN
C
 20   IF (X.LE.8.D0) DBSK1E = (1.25D0 + DCSEVL ((16.D0/X-5.D0)/3.D0,
     1  AK1CS, NTAK1))/SQRT(X)
      IF (X.GT.8.D0) DBSK1E = (1.25D0 +
     1  DCSEVL (16.D0/X-1.D0, AK12CS, NTAK12))/SQRT(X)
C
      RETURN
      END

*DECK DBESI1
      DOUBLE PRECISION FUNCTION DBESI1 (X)
C***BEGIN PROLOGUE  DBESI1
C***PURPOSE  Compute the modified (hyperbolic) Bessel function of the
C            first kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (BESI1-S, DBESI1-D)
C***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESI1(X) calculates the double precision modified (hyperbolic)
C Bessel function of the first kind of order one and double precision
C argument X.
C
C Series for BI1        on the interval  0.          to  9.00000E+00
C                                        with weighted error   1.44E-32
C                                         log weighted error  31.84
C                               significant figures required  31.45
C                                    decimal places required  32.46
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DBSI1E, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBESI1
      DOUBLE PRECISION X, BI1CS(17), XMAX, XMIN, XSML, Y, D1MACH,
     1  DCSEVL, DBSI1E
      LOGICAL FIRST
      SAVE BI1CS, NTI1, XMIN, XSML, XMAX, FIRST
      DATA BI1CS(  1) / -.1971713261 0998597316 1385032181 49 D-2     /
      DATA BI1CS(  2) / +.4073488766 7546480608 1553936520 14 D+0     /
      DATA BI1CS(  3) / +.3483899429 9959455866 2450377837 87 D-1     /
      DATA BI1CS(  4) / +.1545394556 3001236038 5984010584 89 D-2     /
      DATA BI1CS(  5) / +.4188852109 8377784129 4588320041 20 D-4     /
      DATA BI1CS(  6) / +.7649026764 8362114741 9597039660 69 D-6     /
      DATA BI1CS(  7) / +.1004249392 4741178689 1798080372 38 D-7     /
      DATA BI1CS(  8) / +.9932207791 9238106481 3712980548 63 D-10    /
      DATA BI1CS(  9) / +.7663801791 8447637275 2001716813 49 D-12    /
      DATA BI1CS( 10) / +.4741418923 8167394980 3880919481 60 D-14    /
      DATA BI1CS( 11) / +.2404114404 0745181799 8631720320 00 D-16    /
      DATA BI1CS( 12) / +.1017150500 7093713649 1211007999 99 D-18    /
      DATA BI1CS( 13) / +.3645093565 7866949458 4917333333 33 D-21    /
      DATA BI1CS( 14) / +.1120574950 2562039344 8106666666 66 D-23    /
      DATA BI1CS( 15) / +.2987544193 4468088832 0000000000 00 D-26    /
      DATA BI1CS( 16) / +.6973231093 9194709333 3333333333 33 D-29    /
      DATA BI1CS( 17) / +.1436794822 0620800000 0000000000 00 D-31    /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESI1
      IF (FIRST) THEN
         NTI1 = INITDS (BI1CS, 17, 0.1*REAL(D1MACH(3)))
         XMIN = 2.0D0*D1MACH(1)
         XSML = SQRT(4.5D0*D1MACH(3))
         XMAX = LOG (D1MACH(2))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.3.0D0) GO TO 20
C
      DBESI1 = 0.D0
      IF (Y.EQ.0.D0)  RETURN
C
      IF (Y .LE. XMIN) CALL XERMSG ('SLATEC', 'DBESI1',
     +   'ABS(X) SO SMALL I1 UNDERFLOWS', 1, 1)
      IF (Y.GT.XMIN) DBESI1 = 0.5D0*X
      IF (Y.GT.XSML) DBESI1 = X*(0.875D0 + DCSEVL (Y*Y/4.5D0-1.D0,
     1  BI1CS, NTI1))
      RETURN
C
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'DBESI1',
     +   'ABS(X) SO BIG I1 OVERFLOWS', 2, 2)
C
      DBESI1 = EXP(Y) * DBSI1E(X)
C
      RETURN
      END


*DECK DBSI1E
      DOUBLE PRECISION FUNCTION DBSI1E (X)
C***BEGIN PROLOGUE  DBSI1E
C***PURPOSE  Compute the exponentially scaled modified (hyperbolic)
C            Bessel function of the first kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (BESI1E-S, DBSI1E-D)
C***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
C             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
C             ORDER ONE, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBSI1E(X) calculates the double precision exponentially scaled
C modified (hyperbolic) Bessel function of the first kind of order
C one for double precision argument X.  The result is I1(X)
C multiplied by EXP(-ABS(X)).
C
C Series for BI1        on the interval  0.          to  9.00000E+00
C                                        with weighted error   1.44E-32
C                                         log weighted error  31.84
C                               significant figures required  31.45
C                                    decimal places required  32.46
C
C Series for AI1        on the interval  1.25000E-01 to  3.33333E-01
C                                        with weighted error   2.81E-32
C                                         log weighted error  31.55
C                               significant figures required  29.93
C                                    decimal places required  32.38
C
C Series for AI12       on the interval  0.          to  1.25000E-01
C                                        with weighted error   1.83E-32
C                                         log weighted error  31.74
C                               significant figures required  29.97
C                                    decimal places required  32.66
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBSI1E
      DOUBLE PRECISION X, BI1CS(17), AI1CS(46), AI12CS(69), XMIN,
     1  XSML, Y, D1MACH, DCSEVL
      LOGICAL FIRST
      SAVE BI1CS, AI1CS, AI12CS, NTI1, NTAI1, NTAI12, XMIN, XSML,
     1  FIRST
      DATA BI1CS(  1) / -.1971713261 0998597316 1385032181 49 D-2     /
      DATA BI1CS(  2) / +.4073488766 7546480608 1553936520 14 D+0     /
      DATA BI1CS(  3) / +.3483899429 9959455866 2450377837 87 D-1     /
      DATA BI1CS(  4) / +.1545394556 3001236038 5984010584 89 D-2     /
      DATA BI1CS(  5) / +.4188852109 8377784129 4588320041 20 D-4     /
      DATA BI1CS(  6) / +.7649026764 8362114741 9597039660 69 D-6     /
      DATA BI1CS(  7) / +.1004249392 4741178689 1798080372 38 D-7     /
      DATA BI1CS(  8) / +.9932207791 9238106481 3712980548 63 D-10    /
      DATA BI1CS(  9) / +.7663801791 8447637275 2001716813 49 D-12    /
      DATA BI1CS( 10) / +.4741418923 8167394980 3880919481 60 D-14    /
      DATA BI1CS( 11) / +.2404114404 0745181799 8631720320 00 D-16    /
      DATA BI1CS( 12) / +.1017150500 7093713649 1211007999 99 D-18    /
      DATA BI1CS( 13) / +.3645093565 7866949458 4917333333 33 D-21    /
      DATA BI1CS( 14) / +.1120574950 2562039344 8106666666 66 D-23    /
      DATA BI1CS( 15) / +.2987544193 4468088832 0000000000 00 D-26    /
      DATA BI1CS( 16) / +.6973231093 9194709333 3333333333 33 D-29    /
      DATA BI1CS( 17) / +.1436794822 0620800000 0000000000 00 D-31    /
      DATA AI1CS(  1) / -.2846744181 8814786741 0037246830 7 D-1      /
      DATA AI1CS(  2) / -.1922953231 4432206510 4444877497 9 D-1      /
      DATA AI1CS(  3) / -.6115185857 9437889822 5624991778 5 D-3      /
      DATA AI1CS(  4) / -.2069971253 3502277088 8282377797 9 D-4      /
      DATA AI1CS(  5) / +.8585619145 8107255655 3694467313 8 D-5      /
      DATA AI1CS(  6) / +.1049498246 7115908625 1745399786 0 D-5      /
      DATA AI1CS(  7) / -.2918338918 4479022020 9343232669 7 D-6      /
      DATA AI1CS(  8) / -.1559378146 6317390001 6068096907 7 D-7      /
      DATA AI1CS(  9) / +.1318012367 1449447055 2530287390 9 D-7      /
      DATA AI1CS( 10) / -.1448423418 1830783176 3913446781 5 D-8      /
      DATA AI1CS( 11) / -.2908512243 9931420948 2504099301 0 D-9      /
      DATA AI1CS( 12) / +.1266388917 8753823873 1115969040 3 D-9      /
      DATA AI1CS( 13) / -.1664947772 9192206706 2417839858 0 D-10     /
      DATA AI1CS( 14) / -.1666653644 6094329760 9593715499 9 D-11     /
      DATA AI1CS( 15) / +.1242602414 2907682652 3216847201 7 D-11     /
      DATA AI1CS( 16) / -.2731549379 6724323972 5146142863 3 D-12     /
      DATA AI1CS( 17) / +.2023947881 6458037807 0026268898 1 D-13     /
      DATA AI1CS( 18) / +.7307950018 1168836361 9869812612 3 D-14     /
      DATA AI1CS( 19) / -.3332905634 4046749438 1377861713 3 D-14     /
      DATA AI1CS( 20) / +.7175346558 5129537435 4225466567 0 D-15     /
      DATA AI1CS( 21) / -.6982530324 7962563558 5062922365 6 D-16     /
      DATA AI1CS( 22) / -.1299944201 5627607600 6044608058 7 D-16     /
      DATA AI1CS( 23) / +.8120942864 2427988920 5467834286 0 D-17     /
      DATA AI1CS( 24) / -.2194016207 4107368981 5626664378 3 D-17     /
      DATA AI1CS( 25) / +.3630516170 0296548482 7986093233 4 D-18     /
      DATA AI1CS( 26) / -.1695139772 4391041663 0686679039 9 D-19     /
      DATA AI1CS( 27) / -.1288184829 8979078071 1688253822 2 D-19     /
      DATA AI1CS( 28) / +.5694428604 9670527801 0999107310 9 D-20     /
      DATA AI1CS( 29) / -.1459597009 0904800565 4550990028 7 D-20     /
      DATA AI1CS( 30) / +.2514546010 6757173140 8469133448 5 D-21     /
      DATA AI1CS( 31) / -.1844758883 1391248181 6040002901 3 D-22     /
      DATA AI1CS( 32) / -.6339760596 2279486419 2860979199 9 D-23     /
      DATA AI1CS( 33) / +.3461441102 0310111111 0814662656 0 D-23     /
      DATA AI1CS( 34) / -.1017062335 3713935475 9654102357 3 D-23     /
      DATA AI1CS( 35) / +.2149877147 0904314459 6250077866 6 D-24     /
      DATA AI1CS( 36) / -.3045252425 2386764017 4620617386 6 D-25     /
      DATA AI1CS( 37) / +.5238082144 7212859821 7763498666 6 D-27     /
      DATA AI1CS( 38) / +.1443583107 0893824464 1678950399 9 D-26     /
      DATA AI1CS( 39) / -.6121302074 8900427332 0067071999 9 D-27     /
      DATA AI1CS( 40) / +.1700011117 4678184183 4918980266 6 D-27     /
      DATA AI1CS( 41) / -.3596589107 9842441585 3521578666 6 D-28     /
      DATA AI1CS( 42) / +.5448178578 9484185766 5051306666 6 D-29     /
      DATA AI1CS( 43) / -.2731831789 6890849891 6256426666 6 D-30     /
      DATA AI1CS( 44) / -.1858905021 7086007157 7190399999 9 D-30     /
      DATA AI1CS( 45) / +.9212682974 5139334411 2776533333 3 D-31     /
      DATA AI1CS( 46) / -.2813835155 6535611063 7083306666 6 D-31     /
      DATA AI12CS(  1) / +.2857623501 8280120474 4984594846 9 D-1      /
      DATA AI12CS(  2) / -.9761097491 3614684077 6516445730 2 D-2      /
      DATA AI12CS(  3) / -.1105889387 6262371629 1256921277 5 D-3      /
      DATA AI12CS(  4) / -.3882564808 8776903934 5654477627 4 D-5      /
      DATA AI12CS(  5) / -.2512236237 8702089252 9452002212 1 D-6      /
      DATA AI12CS(  6) / -.2631468846 8895195068 3705236523 2 D-7      /
      DATA AI12CS(  7) / -.3835380385 9642370220 4500678796 8 D-8      /
      DATA AI12CS(  8) / -.5589743462 1965838068 6811252222 9 D-9      /
      DATA AI12CS(  9) / -.1897495812 3505412344 9892503323 8 D-10     /
      DATA AI12CS( 10) / +.3252603583 0154882385 5508067994 9 D-10     /
      DATA AI12CS( 11) / +.1412580743 6613781331 6336633284 6 D-10     /
      DATA AI12CS( 12) / +.2035628544 1470895072 2452613684 0 D-11     /
      DATA AI12CS( 13) / -.7198551776 2459085120 9258989044 6 D-12     /
      DATA AI12CS( 14) / -.4083551111 0921973182 2849963969 1 D-12     /
      DATA AI12CS( 15) / -.2101541842 7726643130 1984572746 2 D-13     /
      DATA AI12CS( 16) / +.4272440016 7119513542 9778833699 7 D-13     /
      DATA AI12CS( 17) / +.1042027698 4128802764 1741449994 8 D-13     /
      DATA AI12CS( 18) / -.3814403072 4370078047 6707253539 6 D-14     /
      DATA AI12CS( 19) / -.1880354775 5107824485 1273453396 3 D-14     /
      DATA AI12CS( 20) / +.3308202310 9209282827 3190335240 5 D-15     /
      DATA AI12CS( 21) / +.2962628997 6459501390 6854654205 2 D-15     /
      DATA AI12CS( 22) / -.3209525921 9934239587 7837353288 7 D-16     /
      DATA AI12CS( 23) / -.4650305368 4893583255 7128281897 9 D-16     /
      DATA AI12CS( 24) / +.4414348323 0717079499 4611375964 1 D-17     /
      DATA AI12CS( 25) / +.7517296310 8421048054 2545808029 5 D-17     /
      DATA AI12CS( 26) / -.9314178867 3268833756 8484784515 7 D-18     /
      DATA AI12CS( 27) / -.1242193275 1948909561 1678448869 7 D-17     /
      DATA AI12CS( 28) / +.2414276719 4548484690 0515390217 6 D-18     /
      DATA AI12CS( 29) / +.2026944384 0532851789 7192286069 2 D-18     /
      DATA AI12CS( 30) / -.6394267188 2690977870 4391988681 1 D-19     /
      DATA AI12CS( 31) / -.3049812452 3730958960 8488450357 1 D-19     /
      DATA AI12CS( 32) / +.1612841851 6514802251 3462230769 1 D-19     /
      DATA AI12CS( 33) / +.3560913964 3099250545 1027090462 0 D-20     /
      DATA AI12CS( 34) / -.3752017947 9364390796 6682800324 6 D-20     /
      DATA AI12CS( 35) / -.5787037427 0747993459 5198231074 1 D-22     /
      DATA AI12CS( 36) / +.7759997511 6481619619 8236963209 2 D-21     /
      DATA AI12CS( 37) / -.1452790897 2022333940 6445987408 5 D-21     /
      DATA AI12CS( 38) / -.1318225286 7390367021 2192275337 4 D-21     /
      DATA AI12CS( 39) / +.6116654862 9030707018 7999133171 7 D-22     /
      DATA AI12CS( 40) / +.1376279762 4271264277 3024338363 4 D-22     /
      DATA AI12CS( 41) / -.1690837689 9593478849 1983938230 6 D-22     /
      DATA AI12CS( 42) / +.1430596088 5954331539 8720108538 5 D-23     /
      DATA AI12CS( 43) / +.3409557828 0905940204 0536772990 2 D-23     /
      DATA AI12CS( 44) / -.1309457666 2707602278 4573872642 4 D-23     /
      DATA AI12CS( 45) / -.3940706411 2402574360 9352141755 7 D-24     /
      DATA AI12CS( 46) / +.4277137426 9808765808 0616679735 2 D-24     /
      DATA AI12CS( 47) / -.4424634830 9826068819 0028312302 9 D-25     /
      DATA AI12CS( 48) / -.8734113196 2307149721 1530978874 7 D-25     /
      DATA AI12CS( 49) / +.4045401335 6835333921 4340414242 8 D-25     /
      DATA AI12CS( 50) / +.7067100658 0946894656 5160771780 6 D-26     /
      DATA AI12CS( 51) / -.1249463344 5651052230 0286451860 5 D-25     /
      DATA AI12CS( 52) / +.2867392244 4034370329 7948339142 6 D-26     /
      DATA AI12CS( 53) / +.2044292892 5042926702 8177957421 0 D-26     /
      DATA AI12CS( 54) / -.1518636633 8204625683 7134680291 1 D-26     /
      DATA AI12CS( 55) / +.8110181098 1875758861 3227910703 7 D-28     /
      DATA AI12CS( 56) / +.3580379354 7735860911 2717370327 0 D-27     /
      DATA AI12CS( 57) / -.1692929018 9279025095 9305717544 8 D-27     /
      DATA AI12CS( 58) / -.2222902499 7024276390 6775852777 4 D-28     /
      DATA AI12CS( 59) / +.5424535127 1459696550 4860040112 8 D-28     /
      DATA AI12CS( 60) / -.1787068401 5780186887 6491299330 4 D-28     /
      DATA AI12CS( 61) / -.6565479068 7228149388 2392943788 0 D-29     /
      DATA AI12CS( 62) / +.7807013165 0611452809 2206770683 9 D-29     /
      DATA AI12CS( 63) / -.1816595260 6689797173 7933315222 1 D-29     /
      DATA AI12CS( 64) / -.1287704952 6600848203 7687559895 9 D-29     /
      DATA AI12CS( 65) / +.1114548172 9881645474 1370927369 4 D-29     /
      DATA AI12CS( 66) / -.1808343145 0393369391 5936887668 7 D-30     /
      DATA AI12CS( 67) / -.2231677718 2037719522 3244822893 9 D-30     /
      DATA AI12CS( 68) / +.1619029596 0803415106 1790980361 4 D-30     /
      DATA AI12CS( 69) / -.1834079908 8049414139 0130843921 0 D-31     /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBSI1E
      IF (FIRST) THEN
         ETA = 0.1*REAL(D1MACH(3))
         NTI1 = INITDS (BI1CS, 17, ETA)
         NTAI1 = INITDS (AI1CS, 46, ETA)
         NTAI12 = INITDS (AI12CS, 69, ETA)
C
         XMIN = 2.0D0*D1MACH(1)
         XSML = SQRT(4.5D0*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.3.0D0) GO TO 20
C
      DBSI1E = 0.0D0
      IF (Y.EQ.0.D0)  RETURN
C
      IF (Y .LE. XMIN) CALL XERMSG ('SLATEC', 'DBSI1E',
     +   'ABS(X) SO SMALL I1 UNDERFLOWS', 1, 1)
      IF (Y.GT.XMIN) DBSI1E = 0.5D0*X
      IF (Y.GT.XSML) DBSI1E = X*(0.875D0 + DCSEVL (Y*Y/4.5D0-1.D0,
     1  BI1CS, NTI1) )
      DBSI1E = EXP(-Y) * DBSI1E
      RETURN
C
 20   IF (Y.LE.8.D0) DBSI1E = (0.375D0 + DCSEVL ((48.D0/Y-11.D0)/5.D0,
     1  AI1CS, NTAI1))/SQRT(Y)
      IF (Y.GT.8.D0) DBSI1E = (0.375D0 + DCSEVL (16.D0/Y-1.D0, AI12CS,
     1  NTAI12))/SQRT(Y)
      DBSI1E = SIGN (DBSI1E, X)
C
      RETURN
      END
